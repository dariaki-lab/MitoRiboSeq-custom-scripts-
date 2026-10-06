import importlib.util
import json
import zipfile
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / {"run_pipeline": "01_process_reads.py", "summarize_logs": "summarize_processing_logs.py"}[name])
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


pipeline = module("run_pipeline")
summary = module("summarize_logs")


class PipelineTests(unittest.TestCase):
    def test_each_profile_selects_only_its_twelve_cell_line_samples(self):
        for profile, cell in (("hek_processed", "HEK"), ("hek_unprocessed", "HEK"),
                              ("n2a_processed", "N2a")):
            flow = pipeline.Workflow({"hek_processed": "HEK_settings.json", "hek_unprocessed": "HEK_unprocessed_settings.json", "n2a_processed": "N2a_settings.json"}[profile])
            self.assertEqual(len(flow.samples), 12)
            self.assertEqual({row["cell_line"] for row in flow.samples}, {cell})
            self.assertEqual(len({row["sample_label"] for row in flow.samples}), 12)
            for condition in ("WT", "MTRF1 KO", "MTRF1A KO", "MTRF1/MTRF1A dKO"):
                self.assertEqual({row["replicate"] for row in flow.samples
                                  if row["condition"] == condition}, {"1", "2", "3"})

    def test_cross_cell_line_selection_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown samples"):
            pipeline.Workflow("HEK_settings.json", samples=["N2a_WT1"])

    def test_display_slashes_never_enter_output_file_names(self):
        flow = pipeline.Workflow("N2a_settings.json",
                                 samples=["N2a_MTRF1_MTRF1A_dKO3"])
        row = flow.samples[0]
        self.assertEqual(row["sample_label"], "MTRF1/MTRF1A dKO3")
        self.assertEqual(flow.bam(row).name, "N2a_MTRF1_MTRF1A_dKO3_mapped.sorted.bam")

    def test_processed_and_unprocessed_references_are_separate(self):
        mature = pipeline.Workflow("HEK_settings.json")
        genome = pipeline.Workflow("HEK_unprocessed_settings.json")
        self.assertNotEqual(mature.path("bowtie2_index"), genome.path("bowtie2_index"))
        self.assertEqual(mature.trimmed, genome.trimmed)
        self.assertEqual(len(mature.samples), 12)
        self.assertEqual(len(pipeline.read_samples(mature.sheet_path)), 24)

    def test_alignment_preserves_local_mode_and_mapped_flag_filter(self):
        flow = pipeline.Workflow("N2a_settings.json")
        commands = flow.commands("align", flow.samples[0])
        self.assertIn("--local", commands[0])
        self.assertIn("-U", commands[0])
        filter_cmd = commands[2]
        self.assertEqual(filter_cmd[filter_cmd.index("-F") + 1], "4")
        self.assertFalse(any("--very-sensitive" in command for command in commands))

    def test_profile_read_lengths_and_distinct_cox1_offsets(self):
        flow = pipeline.Workflow("HEK_settings.json")
        self.assertEqual(flow.lengths("phase"), ["--min_length", "30", "--max_length", "45"])
        self.assertEqual(flow.lengths("phase_cox1"), ["--min_length", "25", "--max_length", "45"])
        self.assertNotEqual(flow.offset(flow.samples[0]), flow.offset(flow.samples[0], "phase-cox1"))
        genome = pipeline.Workflow("HEK_unprocessed_settings.json")
        self.assertEqual(genome.lengths("phase"), ["--min_length", "25", "--max_length", "45"])
        self.assertEqual(genome.lengths("vectors"), ["--min_length", "30", "--max_length", "45"])

    def test_incomplete_offsets_fail_instead_of_defaulting_silently(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            path = Path(tmp) / "offsets.txt"
            path.write_text("30\t19\n31\t18\n")
            with self.assertRaisesRegex(ValueError, "offsets missing"):
                pipeline.validate_offsets(path, 30, 45)
            path.write_text("default\t19\n")
            pipeline.validate_offsets(path, 30, 45)
            path.write_text("30\t30\n")
            with self.assertRaisesRegex(ValueError, "invalid offset"):
                pipeline.validate_offsets(path, 30, 30)

    def test_seven_recovered_tables_are_complete(self):
        with zipfile.ZipFile(ROOT / "original_workflows.zip") as archive:
            files = [name for name in archive.namelist() if name.startswith("N2a_recovered_offset_")]
            self.assertEqual(len(files), 7)
            with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
                for name in files:
                    file = Path(tmp) / name
                    file.write_bytes(archive.read(name))
                    pipeline.validate_offsets(file, 30, 45)
            wt2 = archive.read("N2a_recovered_offset_N2a_WT_2.txt").decode()
            self.assertIn("32\t20\n", wt2)
            self.assertIn("33\t18\n", wt2)

    def test_unconfirmed_settings_block_execution(self):
        flow = pipeline.Workflow("N2a_settings.json")
        with self.assertRaisesRegex(ValueError, "author-confirmed"):
            flow.preflight("align")

    def test_bad_sample_name_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            path = Path(tmp) / "samples.tsv"
            path.write_text("sample_id\tfastq_path\tcondition\tbatch\toffset_file\n"
                            "../outside\treads.fastq.gz\tWT\tDK\t\n")
            with self.assertRaisesRegex(ValueError, "Sample identifiers"):
                pipeline.read_samples(path)

    def test_unknown_sample_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown samples"):
            pipeline.Workflow("HEK_settings.json", samples=["wrong"])

    def test_failed_bowtie2_prevents_bam_processing_and_final_outputs(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            flow = pipeline.Workflow("N2a_settings.json", samples=["N2a_WT1"])
            flow.output = Path(tmp) / "results"
            flow.trimmed = Path(tmp) / "trimmed"
            # Isolate command ordering; this does not emulate alignment or scientific output.
            with patch.object(flow, "preflight"), patch.object(flow, "write_provenance"):
                (flow.output / "logs").mkdir(parents=True)
                with patch.object(pipeline.subprocess, "run",
                                  side_effect=subprocess.CalledProcessError(2, ["bowtie2"])) as run:
                    with self.assertRaises(subprocess.CalledProcessError):
                        flow.run("align")
                    self.assertEqual(run.call_count, 1)
                self.assertFalse(flow.bam(flow.samples[0]).exists())
                self.assertFalse(Path(str(flow.bam(flow.samples[0])) + ".bai").exists())

    def test_dry_run_does_not_create_results(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmp:
            flow = pipeline.Workflow("HEK_settings.json")
            flow.output = Path(tmp) / "unused"
            with patch("sys.stdout"):
                flow.run("align", dry_run=True)
            self.assertFalse(flow.output.exists())

    def test_log_summaries_handle_commas_and_missing_entries(self):
        cut = summary.parse_cutadapt(
            "Total reads processed: 35,489,635\nReads with adapters: 34,973,949 (98.5%)\n"
            "Reads written (passing filters): 23,020,016 (64.9%)\n")
        self.assertEqual(cut["Total"], 35489635)
        self.assertEqual(cut["Written"], 23020016)
        bow = summary.parse_bowtie2(
            "23020016 (100.00%) were unpaired; of these:\n"
            "  21042294 (91.41%) aligned 0 times\n"
            "  1977715 (8.59%) aligned exactly 1 time\n"
            "  7 (0.00%) aligned >1 times\n")
        self.assertEqual(sum(bow[k] for k in ("Unaligned","Aligned_1x","Aligned_multi")), bow["Total"])
        self.assertIsNone(summary.parse_cutadapt("partial log")["Written"])


if __name__ == "__main__":
    unittest.main()

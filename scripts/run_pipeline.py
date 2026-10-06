#!/usr/bin/env python3
"""Run the supplied mitoriboseq processing commands with explicit inputs.

This is a refactoring of the attached Cutadapt/Bowtie2/plastid workflows.
It does not choose biological settings or infer missing references/offsets.
Python 3.9+; only standard-library modules are used by this runner.
"""
import argparse
import csv
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGES = ("index", "trim", "align", "psite", "phase", "phase-cox1", "counts", "vectors")
TOOLS = {
    "index": ["bowtie2-build"],
    "trim": ["cutadapt"],
    "align": ["bowtie2", "samtools"],
    "psite": ["psite"],
    "phase": ["phase_by_size"],
    "phase-cox1": ["phase_by_size"],
    "counts": ["counts_in_region"],
    "vectors": ["get_count_vectors"],
}


def resolve(root, value):
    p = Path(value)
    return p if p.is_absolute() else root / p


def require_file(p):
    if not p.is_file() or p.stat().st_size == 0:
        raise ValueError(f"Missing or empty input: {p}")


def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validate_offsets(path, lower, upper):
    """Require an offset for each retained length, or an explicit default."""
    require_file(path)
    offsets = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(f"{path}:{number}: expected length and offset")
        length, value = fields
        try:
            key = length if length == "default" else int(length)
            offset = int(value)
        except ValueError as error:
            raise ValueError(f"{path}:{number}: non-integer offset entry") from error
        if key in offsets or offset < 0 or (key != "default" and offset >= key):
            raise ValueError(f"{path}:{number}: duplicate or invalid offset")
        offsets[key] = offset
    missing = [n for n in range(lower, upper + 1) if n not in offsets and "default" not in offsets]
    if missing:
        raise ValueError(f"{path}: offsets missing for lengths {missing}")
    for n in range(lower, upper + 1):
        value = offsets.get(n, offsets.get("default"))
        if value is not None and value >= n:
            raise ValueError(f"{path}: offset {value} is outside read length {n}")


def read_samples(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="," if path.suffix.lower() == ".csv" else "\t")
        required = {"sample_id", "fastq_path", "condition", "batch", "offset_file"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing sample-sheet columns {required}")
        rows = list(reader)
    identifiers = [row["sample_id"] for row in rows]
    if not rows or len(set(identifiers)) != len(identifiers):
        raise ValueError("Sample identifiers must be nonempty and unique.")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", x) for x in identifiers):
        raise ValueError("Sample identifiers must contain only letters, digits, _, . or -.")
    return rows


class Workflow:
    def __init__(self, config_path, root=ROOT, samples=None):
        self.root = Path(root).resolve()
        self.config_path = resolve(self.root, str(config_path))
        self.config = json.loads(self.config_path.read_text())
        self.sheet_path = resolve(self.root, self.config["sample_sheet"])
        self.samples = read_samples(self.sheet_path)
        if self.config.get("cell_line"):
            self.samples = [row for row in self.samples
                            if row.get("cell_line") == self.config["cell_line"]]
        if samples:
            requested = set(samples)
            unknown = requested - {x["sample_id"] for x in self.samples}
            if unknown:
                raise ValueError(f"Unknown samples: {sorted(unknown)}")
            self.samples = [row for row in self.samples if row["sample_id"] in requested]
        elif self.config.get("processing_batches"):
            self.samples = [row for row in self.samples if row["batch"] in self.config["processing_batches"]]
        if not self.samples:
            raise ValueError("No samples selected for processing.")
        self.output = resolve(self.root, self.config["output_dir"])
        self.trimmed = resolve(self.root, self.config["trimmed_dir"])
        if not isinstance(self.config["threads"], int) or self.config["threads"] < 1:
            raise ValueError("threads must be a positive integer")
        for key in ("trim", "psite", "phase", "counts", "vectors", "phase_cox1"):
            lower, upper = self.config[key]["min_length"], self.config[key]["max_length"]
            if not (isinstance(lower, int) and isinstance(upper, int) and 1 <= lower <= upper):
                raise ValueError(f"Invalid read-length range: {key}")

    def path(self, key):
        return resolve(self.root, self.config[key])

    def bam(self, row):
        return self.output / "bam" / (row["sample_id"] + "_mapped.sorted.bam")

    def offset(self, row, stage=None):
        value = self.config.get("cox1_offset_file") if stage == "phase-cox1" else None
        value = value or row["offset_file"] or self.config["offset_file"]
        return resolve(self.root, value.format(sample_id=row["sample_id"], analysis_id=self.config["analysis_id"]))

    def trimmed_fastq(self, row):
        return self.trimmed / (row["sample_id"] + "_cut.fastq.gz")

    def preflight(self, stage, force=False):
        cfg = self.config
        if cfg["settings_confirmed"] is not True:
            raise ValueError(
                "Settings are recorded but not author-confirmed. Review docs/RELEASE_CHECKLIST.md "
                "and set settings_confirmed=true in your working config before execution."
            )
        for executable in TOOLS[stage]:
            if shutil.which(executable) is None:
                raise ValueError(f"Required program is not on PATH: {executable}")
        if stage == "index":
            require_file(self.path("reference_fasta"))
            existing = list(self.path("bowtie2_index").parent.glob(self.path("bowtie2_index").name + ".*.bt2*"))
            if existing and not force:
                raise ValueError("Index files already exist. Use --force after reviewing the reference.")
            return
        if stage == "align":
            prefix = self.path("bowtie2_index")
            for suffix in (".bt2", ".bt2l"):
                parts = [Path(str(prefix) + "." + part + suffix) for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
                if all(p.is_file() and p.stat().st_size for p in parts):
                    break
            else:
                raise ValueError(f"No complete Bowtie2 index at {prefix}; run the index stage.")
        if stage in ("psite", "phase", "phase-cox1", "counts", "vectors"):
            if stage != "psite" and cfg["mapping_site"] not in ("P", "A", "custom"):
                raise ValueError("Set mapping_site to P, A, or custom after checking the approved offsets.")
            key = {"psite": cfg.get("offset_roi_key", "roi_stop"), "phase": "roi_start", "phase-cox1": "roi_cox1",
                   "counts": "annotation_gtf", "vectors": "annotation_bed"}[stage]
            require_file(self.path(key))
        for row in self.samples:
            if stage == "trim":
                if not row["fastq_path"]:
                    raise ValueError(f"Add the confirmed fastq_path for {row['sample_id']} in {self.sheet_path}")
                require_file(resolve(self.root, row["fastq_path"]))
            elif stage == "align":
                require_file(self.trimmed_fastq(row))
            else:
                bam = self.bam(row)
                require_file(bam)
                require_file(Path(str(bam) + ".bai"))
                if stage != "psite":
                    options = cfg["phase_cox1" if stage == "phase-cox1" else stage]
                    validate_offsets(self.offset(row, stage), options["min_length"], options["max_length"])
            # A failed previous run leaves its log; require an explicit rerun.
            log = self.output / "logs" / (row["sample_id"] + "_" + stage + ".log")
            if log.exists() and not force:
                raise ValueError(f"Existing run log: {log}. Use --force for an intentional rerun.")
            primary = {"trim": self.trimmed_fastq(row), "align": self.bam(row),
                       "psite": self.output / "offset_estimates" / (row["sample_id"] + "_p_offsets.txt"),
                       "phase": self.output / "phasing" / (row["sample_id"] + "_phasing.txt"),
                       "phase-cox1": self.output / "phasing_COX1" / (row["sample_id"] + "_phasing.txt"),
                       "counts": self.output / "counts" / ("Counts_" + row["sample_id"] + ".txt"),
                       "vectors": self.output / "vectors" / row["sample_id"]}[stage]
            if primary.exists() and not force:
                raise ValueError(f"Existing output: {primary}. Use --force for an intentional rerun.")

    def lengths(self, key):
        return ["--min_length", str(self.config[key]["min_length"]),
                "--max_length", str(self.config[key]["max_length"])]

    def commands(self, stage, row=None, work=None):
        """Return lists of argv; shell expansion and evaluation are never used."""
        cfg = self.config
        if stage == "index":
            return [["bowtie2-build", str(self.path("reference_fasta")), str(self.path("bowtie2_index"))]]
        sample = row["sample_id"]
        work = Path(work) if work else self.output / "_work" / sample
        if stage == "trim":
            source = str(resolve(self.root, row["fastq_path"])) if row["fastq_path"] else f"<ADD_FASTQ_PATH_FOR_{sample}>"
            cmd = ["cutadapt", "-a", cfg["adapter"], "-o", str(work / "trimmed.fastq.gz"), source,
                   "-m", str(cfg["trim"]["min_length"]), "-M", str(cfg["trim"]["max_length"]),
                   "--cores", str(cfg["threads"])]
            if cfg["max_n"] is not None:
                cmd += ["--max-n", str(cfg["max_n"])]
            return [cmd]
        if stage == "align":
            return [
                ["bowtie2", "--local", "-x", str(self.path("bowtie2_index")),
                 "-U", str(self.trimmed_fastq(row)), "-p", str(cfg["threads"]), "-S", str(work / "alignment.sam")],
                ["samtools", "view", "-bS", "-o", str(work / "all.bam"), str(work / "alignment.sam")],
                ["samtools", "view", "-b", "-F", "4", "-o", str(work / "mapped.bam"), str(work / "all.bam")],
                ["samtools", "sort", str(work / "mapped.bam"), "-o", str(work / "sorted.bam")],
                ["samtools", "index", str(work / "sorted.bam")],
                ["samtools", "quickcheck", str(work / "sorted.bam")],
            ]
        if stage == "psite":
            return [["psite", str(self.path(cfg.get("offset_roi_key", "roi_stop"))), str(self.output / "offset_estimates" / sample),
                     *self.lengths("psite"), "--count_files", str(self.bam(row))]]
        common = ["--count_files", str(self.bam(row)), "--fiveprime_variable", "--offset", str(self.offset(row, stage))]
        if stage in ("phase", "phase-cox1"):
            cox1 = stage == "phase-cox1"
            destination = self.output / ("phasing_COX1" if cox1 else "phasing") / sample
            return [["phase_by_size", str(self.path("roi_cox1" if cox1 else "roi_start")),
                     str(destination), *common, "--codon_buffer", str(cfg["codon_buffer"]),
                     *self.lengths("phase_cox1" if cox1 else "phase")]]
        if stage == "counts":
            return [["counts_in_region", str(self.output / "counts" / ("Counts_" + sample + ".txt")),
                     *common, "--annotation_files", str(self.path("annotation_gtf")), *self.lengths("counts")]]
        return [["get_count_vectors", "--annotation_files", str(self.path("annotation_bed")),
                 "--annotation_format", "BED", *common, *self.lengths("vectors"),
                 str(self.output / "vectors" / sample)]]

    def write_provenance(self, stage, force):
        self.output.mkdir(parents=True, exist_ok=True)
        logdir = self.output / "logs"
        logdir.mkdir(exist_ok=True)
        timestamp = datetime.now(timezone.utc).isoformat()
        record = {"time_utc": timestamp, "stage": stage, "force": force,
                  "config": self.config, "config_sha256": digest(self.config_path),
                  "sample_sheet_sha256": digest(self.sheet_path),
                  "samples": [row["sample_id"] for row in self.samples], "runner_sha256": digest(Path(__file__)),
                  "references": {}}
        for key in ("reference_fasta", "annotation_gtf", "annotation_bed", "roi_start", "roi_stop", "roi_cox1"):
            p = self.path(key)
            if p.is_file():
                record["references"][key] = {"path": str(p), "sha256": digest(p)}
        if stage in ("phase", "phase-cox1", "counts", "vectors"):
            record["approved_offsets"] = {row["sample_id"]: {"path": str(self.offset(row, stage)),
                                                         "sha256": digest(self.offset(row, stage))}
                                          for row in self.samples}
        with (logdir / "runs.jsonl").open("a") as handle:
            handle.write(json.dumps(record) + "\n")
        with (logdir / "versions.txt").open("a") as handle:
            handle.write(f"\nRun {timestamp}; stage {stage}\n")
            for executable in TOOLS[stage]:
                try:
                    result = subprocess.run([executable, "--version"], text=True, capture_output=True, timeout=20)
                    handle.write(f"{executable}; exit={result.returncode}\n{result.stdout}{result.stderr}\n")
                except (OSError, subprocess.TimeoutExpired) as error:
                    handle.write(f"{executable}: version query unavailable: {error}\n")

    def run(self, stage, dry_run=False, force=False):
        if dry_run:
            if not self.config["settings_confirmed"]:
                print("# Preview only: settings still need author confirmation.")
            for row in ([None] if stage == "index" else self.samples):
                for cmd in self.commands(stage, row):
                    print(shlex.join(cmd))
            return
        self.preflight(stage, force)
        self.write_provenance(stage, force)
        if stage == "index":
            self.path("bowtie2_index").parent.mkdir(parents=True, exist_ok=True)
            self.execute(self.commands(stage)[0], self.output / "logs" / "index.log")
            return
        for row in self.samples:
            sample = row["sample_id"]
            for name in ("bam", "_work", "offset_estimates", "phasing", "phasing_COX1", "counts", "vectors"):
                (self.output / name).mkdir(parents=True, exist_ok=True)
            self.trimmed.mkdir(parents=True, exist_ok=True)
            log = self.output / "logs" / (sample + "_" + stage + ".log")
            if force and log.exists():
                log.unlink()
            with tempfile.TemporaryDirectory(prefix=sample + "_", dir=self.output / "_work") as tmp:
                work = Path(tmp)
                for cmd in self.commands(stage, row, work):
                    self.execute(cmd, log)
                if stage == "trim":
                    require_file(work / "trimmed.fastq.gz")
                    (work / "trimmed.fastq.gz").replace(self.trimmed_fastq(row))
                elif stage == "align":
                    require_file(work / "sorted.bam")
                    require_file(work / "sorted.bam.bai")
                    (work / "sorted.bam").replace(self.bam(row))
                    (work / "sorted.bam.bai").replace(Path(str(self.bam(row)) + ".bai"))
                    if self.config["keep_intermediate_bams"]:
                        for name, suffix in (("all.bam", ".bam"), ("mapped.bam", "_mapped.bam")):
                            (work / name).replace(self.output / "bam" / (sample + suffix))
            print(f"Completed {sample}: {stage}")

    def execute(self, command, log):
        print(shlex.join(command), flush=True)
        with log.open("a") as handle:
            handle.write("\nCOMMAND " + shlex.join(command) + "\n")
            handle.flush()
            with (self.output / "logs" / "commands.jsonl").open("a") as trace:
                trace.write(json.dumps({"argv": command, "log": str(log)}) + "\n")
            # Separate commands replace the source pipelines so upstream errors are observed.
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True, cwd=self.root)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON profile; relative paths use the repository root.")
    parser.add_argument("--stage", required=True, choices=STAGES)
    parser.add_argument("--root", type=Path, default=ROOT, help="Override repository root.")
    parser.add_argument("--samples", nargs="+", help="Run only these sample_id values.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running or writing files.")
    parser.add_argument("--force", action="store_true", help="Intentionally replace outputs from this stage.")
    args = parser.parse_args(argv)
    try:
        Workflow(args.config, args.root, args.samples).run(args.stage, args.dry_run, args.force)
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

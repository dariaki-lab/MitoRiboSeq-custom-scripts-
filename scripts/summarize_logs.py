#!/usr/bin/env python3
"""Recover Cutadapt/Bowtie2 summaries without assuming numeric field positions."""
import argparse
import csv
import json
import re
from pathlib import Path


def extract(text, label):
    match = re.search(label, text, flags=re.MULTILINE)
    return int(match.group(1).replace(",", "")) if match else None


def parse_cutadapt(text):
    return {
        "Total": extract(text, r"Total reads processed:\s*([\d,]+)"),
        "WithAdapters": extract(text, r"Reads with adapters:\s*([\d,]+)"),
        "Written": extract(text, r"Reads written \(passing filters\):\s*([\d,]+)"),
    }


def parse_bowtie2(text):
    return {
        "Total": extract(text, r"^\s*([\d,]+)\s+\([^)]+\) were unpaired"),
        "Unaligned": extract(text, r"^\s*([\d,]+)\s+\([^)]+\) aligned 0 times"),
        "Aligned_1x": extract(text, r"^\s*([\d,]+)\s+\([^)]+\) aligned exactly 1 time"),
        "Aligned_multi": extract(text, r"^\s*([\d,]+)\s+\([^)]+\) aligned >1 times"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    cfgpath = Path(args.config)
    cfg = json.loads((cfgpath if cfgpath.is_absolute() else root / cfgpath).read_text())
    sheet = Path(cfg["sample_sheet"])
    with (sheet if sheet.is_absolute() else root / sheet).open(newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="," if sheet.suffix.lower() == ".csv" else "\t"))
    if cfg.get("cell_line"):
        samples = [row for row in samples if row.get("cell_line") == cfg["cell_line"]]
    if cfg.get("processing_batches"):
        samples = [row for row in samples if row["batch"] in cfg["processing_batches"]]
    output = Path(cfg["output_dir"])
    output = output if output.is_absolute() else root / output
    output.mkdir(parents=True, exist_ok=True)
    for stage, name, parse, columns in (
        ("trim", "cutadapt_summary.tsv", parse_cutadapt, ["Total", "WithAdapters", "Written"]),
        ("align", "bowtie2_summary.tsv", parse_bowtie2, ["Total", "Unaligned", "Aligned_1x", "Aligned_multi"]),
    ):
        with (output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Sample", *columns, "Status"], delimiter="\t")
            writer.writeheader()
            for sample in samples:
                path = output / "logs" / (sample["sample_id"] + "_" + stage + ".log")
                parsed = parse(path.read_text()) if path.exists() else dict.fromkeys(columns)
                status = "complete" if all(v is not None for v in parsed.values()) else ("missing_log" if not path.exists() else "incomplete_log")
                writer.writerow({"Sample": sample["sample_id"], **parsed, "Status": status})
        print(output / name)


if __name__ == "__main__":
    main()

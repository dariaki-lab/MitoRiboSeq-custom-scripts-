#!/usr/bin/env python3
"""Check the general repository template without sequencing software or inputs."""
import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONDITIONS = {
    "WT": "WT",
    "MTRF1_KO": "MTRF1 KO",
    "MTRF1A_KO": "MTRF1A KO",
    "MTRF1_MTRF1A_dKO": "MTRF1/MTRF1A dKO",
}


def rows(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta(path):
    sequences = {}
    key = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            if key in sequences:
                raise ValueError(f"Duplicate FASTA identifier: {key}")
            sequences[key] = ""
        elif line.strip():
            if key is None:
                raise ValueError(f"Sequence before FASTA header: {path}")
            sequences[key] += line.strip()
    return sequences


def validate(root=ROOT, check_checksums=True):
    root = Path(root)
    errors = []
    samples = rows(root / "metadata/sample_info.csv", ",")
    expected = {f"{cell}_{code}{rep}": (cell, label, str(rep), f"{label}{rep}")
                for cell in ("HEK", "N2a") for code, label in CONDITIONS.items()
                for rep in range(1, 4)}
    if len(samples) != 24 or {s["sample_id"] for s in samples} != set(expected):
        errors.append("The general template must contain exactly the 24 requested sample IDs.")
    for sample in samples:
        if (sample["cell_line"], sample["condition"], sample["replicate"], sample["sample_label"]) != expected.get(sample["sample_id"]):
            errors.append(f"Inconsistent sample label/condition: {sample['sample_id']}")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample["sample_id"]):
            errors.append(f"Unsafe file ID: {sample['sample_id']}")
    profiles = []
    for cfgfile in sorted((root / "config").glob("*.json")):
        cfg = json.loads(cfgfile.read_text())
        profiles.append(cfg["analysis_id"])
        selected = [s for s in samples if s["cell_line"] == cfg["cell_line"]]
        if len(selected) != 12:
            errors.append(f"{cfgfile.name}: expected 12 samples for its cell line")
        if cfg["sample_sheet"] != "metadata/sample_info.csv":
            errors.append(f"{cfgfile.name}: not using the common sample table")
        mapping = rows(root / cfg["r_analysis"]["vector_mapping"])
        if len(mapping) != 13 or len({m["Gene"] for m in mapping}) != 13:
            errors.append(f"{cfgfile.name}: expected 13 distinct transcript mappings")
        sequence_path = cfg["r_analysis"].get("sequence_fasta")
        if sequence_path:
            seqs = read_fasta(root / sequence_path)
            if set(seqs) != {m["sequence_gene"] for m in mapping}:
                errors.append(f"{cfgfile.name}: sequence/mapping identifiers disagree")
            if any(not re.fullmatch(r"[ACGTN]+", seq) for seq in seqs.values()):
                errors.append(f"{cfgfile.name}: invalid sequence characters")
        if cfg["cell_line"] == "N2a":
            for key in ("annotation_gtf", "annotation_bed", "roi_start", "roi_stop", "roi_cox1"):
                path = root / cfg[key]
                if path.is_file() and "Human_mtDNA" in path.read_text():
                    errors.append(f"Human coordinates in N2a input: {path.name}")
    for ref in rows(root / "metadata/reference_manifest.tsv"):
        if ref["status"].startswith("INCLUDED"):
            path = root / ref["repository_path"]
            if not path.is_file() or sha256(path) != ref["sha256"]:
                errors.append(f"Included reference changed or missing: {ref['repository_path']}")
    for coordinate in rows(root / "metadata/sequence_coordinate_review.tsv"):
        seqs = read_fasta(root / "references" / coordinate["profile"] / "analysis_sequences.fasta")
        seq = seqs[coordinate["Gene"]]
        if len(seq) != int(coordinate["sequence_length"]) or hashlib.sha256(seq.encode()).hexdigest() != coordinate["sequence_sha256"]:
            errors.append(f"Recovered sequence differs: {coordinate['profile']} {coordinate['Gene']}")
    # Active documentation links only; historical archive links are source records.
    for doc in root.rglob("*.md"):
        if "archive" in doc.relative_to(root).parts:
            continue
        for target in re.findall(r"!?\[[^\]]*\]\(([^\s)]+)\)", doc.read_text()):
            if re.match(r"[A-Za-z][A-Za-z0-9+.-]*:", target) or target.startswith("#"):
                continue
            target = target.split("#", 1)[0]
            if target and not (doc.parent / target).exists():
                errors.append(f"Broken link in {doc.relative_to(root)}: {target}")
    # Public copies should not retain the original user's machine paths.
    for path in root.rglob("*"):
        if path.is_file() and path.suffix in (".py", ".R", ".Rmd", ".txt", ".md", ".gtf2", ".bed"):
            if re.search(r"/Users/[A-Za-z0-9]", path.read_text()):
                errors.append(f"Personal absolute path: {path.relative_to(root)}")
    checksums = root / "metadata/package_files.sha256"
    if check_checksums and checksums.is_file():
        for line in checksums.read_text().splitlines():
            expected_hash, relative = line.split("  ", 1)
            path = root / relative
            if not path.is_file() or sha256(path) != expected_hash:
                errors.append(f"Package checksum mismatch: {relative}")
    missing = [r["repository_path"] for r in rows(root / "metadata/reference_manifest.tsv")
               if not r["status"].startswith("INCLUDED")]
    return {"ok": not errors, "sample_slots": len(samples), "profiles": profiles,
            "errors": errors, "experimental_inputs_still_to_supply": missing,
            "scientific_results_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--skip-checksums", action="store_true",
                        help="After adapting files, validate structure without original package hashes.")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    try:
        report = validate(args.root, not args.skip_checksums)
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"{'PASS' if report['ok'] else 'FAIL'}: {report['sample_slots']} sample slots; "
              f"{len(report['profiles'])} analysis profiles.")
        for error in report["errors"]:
            print("ERROR:", error)
        print(f"Experimental inputs to supply: {len(report['experimental_inputs_still_to_supply'])} manifest entries.")
        print("This checks the general template and included files, not manuscript results.")
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

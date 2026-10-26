#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 samples.csv output_directory" >&2
    echo "CSV header: sample,hap1,hap2 (relative FASTA paths resolve from the CSV directory)" >&2
    exit 1
fi

mkdir -p "$2"
outdir=$(cd "$2" && pwd)
scratch=$(mktemp -d "$outdir/.pansn.XXXXXXXX")
trap 'rm -f "$scratch/manifest" "$scratch/output.fa.gz"; rmdir "$scratch"' EXIT

# Parse quoted CSV fields without treating commas or spaces as shell separators.
python3 - "$1" "$outdir" > "$scratch/manifest" <<'PY'
import csv
from pathlib import Path
import re
import sys

csv_path = Path(sys.argv[1]).resolve()
outdir = Path(sys.argv[2])
seen = set()
with csv_path.open(newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)
    if not {"sample", "hap1", "hap2"}.issubset(reader.fieldnames or []):
        sys.exit("CSV must have columns: sample,hap1,hap2")
    for row in reader:
        sample = (row["sample"] or "").strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            sys.exit(f"Invalid sample name at CSV line {reader.line_num}: {sample!r}")
        if sample in seen:
            sys.exit(f"Duplicate sample: {sample}")
        seen.add(sample)
        for hap in (1, 2):
            source = (row[f"hap{hap}"] or "").strip()
            if not source or "\0" in source:
                sys.exit(f"Invalid hap{hap} path for {sample}")
            path = Path(source)
            if not path.is_absolute():
                path = csv_path.parent / path
            path = path.resolve()
            if not path.is_file():
                sys.exit(f"FASTA not found: {path}")
            target = outdir / f"{sample}_hap{hap}.fa.gz"
            if target.exists() or target.is_symlink():
                sys.exit(f"Refusing to overwrite: {target}")
            for field in (sample, str(hap), str(path)):
                sys.stdout.buffer.write(field.encode() + b"\0")
if not seen:
    sys.exit("CSV contains no samples")
PY

while IFS= read -r -d '' sample &&
      IFS= read -r -d '' hap &&
      IFS= read -r -d '' source; do
    target="$outdir/${sample}_hap${hap}.fa.gz"
    gzip -cdf -- "$source" |
        awk -v prefix="$sample#$hap#" '
            /^>/ {
                if (length($1) == 1 || index($1, "#")) {
                    print "Empty contig name or contig name containing # is not supported" > "/dev/stderr"
                    invalid = 1
                    exit 1
                }
                sub(/^>/, ">" prefix)
                headers++
            }
            { print }
            END {
                if (!headers && !invalid) {
                    print "Input contains no FASTA headers" > "/dev/stderr"
                    exit 1
                }
            }
        ' | gzip -c > "$scratch/output.fa.gz"
    mv "$scratch/output.fa.gz" "$target"
    printf '%s\n' "$target"
done < "$scratch/manifest"

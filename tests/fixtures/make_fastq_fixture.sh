#!/usr/bin/env bash
# make_fastq_fixture.sh -- rebuild the phenotyping test FASTQs (tests/fixtures/fastq) and their
# samplesheet/metadata (tests/fixtures/phenotyping) from snuh_chr22_fastq.tar.
#
# Source: 8 placentas of the SNUH cohort (BioProject PRJNA820329, public). The archive is written by
# tests/fixtures/extract_snuh_chr22.sh from the atlas STAR BAMs (reads already Trim Galore-trimmed):
# read pairs with both mates aligned to chr22, about 170k-275k pairs per sample.
#
# Reduction for the repository: base qualities are replaced by a constant '?' (Phred 30), which is
# how ENA/SRA distribute this project's FASTQs (SRA Lite). Sequences, read names and read order are
# unchanged. Output is deterministic (gzip level 9, mtime 0), so a rebuild gives identical files;
# md5s are recorded in tests/fixtures/fastq/MD5SUMS.
#
# Metadata: sex and gestational age (days / 7, one decimal) from each run's BioSample record;
# genotype_id is the BioSample isolate name (no genotypes are used by the phenotyping test).
#
# Usage: tests/fixtures/make_fastq_fixture.sh snuh_chr22_fastq.tar
set -euo pipefail

TAR=${1:?usage: make_fastq_fixture.sh snuh_chr22_fastq.tar}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
tar -C "$TMP" -xf "$TAR"

python3 - "$TMP/snuh_chr22" "$HERE" <<'EOF'
import csv, gzip, hashlib, os, sys

src, here = sys.argv[1], sys.argv[2]
# run accession -> (BioSample, isolate, sex, gestational days); from https://www.ebi.ac.uk/biosamples
BIOSAMPLE = {
    "SRR18553497": ("SAMN26993873", "OGF212", "male", 301),
    "SRR18553498": ("SAMN26993872", "OGF209", "male", 296),
    "SRR18553501": ("SAMN26993869", "OGF206", "male", 296),
    "SRR18553511": ("SAMN26993860", "OGF172", "male", 275),
    "SRR18553512": ("SAMN26993859", "OGF171", "male", 299),
    "SRR18553513": ("SAMN26993787", "OGF056", "male", 294),
    "SRR18553520": ("SAMN26993852", "OGF157", "male", 310),
    "SRR18553521": ("SAMN26993851", "OGF155", "male", 281),
}

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = list(csv.DictReader(open(os.path.join(src, "manifest.tsv")), delimiter="\t"))
fqdir = os.path.join(here, "fastq")
phdir = os.path.join(here, "phenotyping")
os.makedirs(fqdir, exist_ok=True)
os.makedirs(phdir, exist_ok=True)
sums, rows = [], []
for m in manifest:
    acc = m["sample"]
    if acc not in BIOSAMPLE:
        sys.exit(f"{acc}: no BioSample record in this script")
    n_pairs = None
    for r in ("R1", "R2"):
        inp = os.path.join(src, f"{acc}_{r}.fastq.gz")
        if md5(inp) != m[f"{r}_md5"]:
            sys.exit(f"{inp}: md5 differs from manifest.tsv")
        out = os.path.join(fqdir, f"{acc}_{r}.fastq.gz")
        n = 0
        with gzip.open(inp, "rt") as f, open(out, "wb") as raw, \
                gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as g:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline().rstrip("\n")
                f.readline()
                f.readline()
                g.write(f"{h.rstrip()}\n{s}\n+\n{'?' * len(s)}\n".encode())
                n += 1
        if n != int(m["pairs_written"]):
            sys.exit(f"{inp}: {n} reads, manifest says {m['pairs_written']}")
        n_pairs = n
        sums.append(f"{md5(out)}  {acc}_{r}.fastq.gz")
    bs, iso, sex, days = BIOSAMPLE[acc]
    rows.append((acc, iso, sex, round(days / 7, 1), bs, n_pairs))
    print(f"{acc}: {n_pairs} pairs", flush=True)

with open(os.path.join(fqdir, "MD5SUMS"), "w") as f:
    f.write("\n".join(sums) + "\n")
with open(os.path.join(phdir, "samplesheet.csv"), "w") as f:
    f.write("sample_id,fastq_1,fastq_2\n")
    for acc, *_ in rows:
        f.write(f"{acc},../fastq/{acc}_R1.fastq.gz,../fastq/{acc}_R2.fastq.gz\n")
with open(os.path.join(phdir, "metadata.tsv"), "w") as f:
    f.write("rnaseq_id\tgenotype_id\tsex\tGA\tbiosample\tchr22_pairs\n")
    for acc, iso, sex, ga, bs, n in rows:
        f.write(f"{acc}\t{iso}\t{sex}\t{ga}\t{bs}\t{n}\n")
EOF
du -sh "$HERE/fastq"

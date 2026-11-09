#!/usr/bin/env python3
"""check_imputation_prep.py lifted|swapped OUTDIR GRCh38_VCF  -- checks for test_imputation_prep.sh.

Both: chr22.vcf.gz (+ .tbi), imputation_prep_report.tsv, imputation_prep_samples.tsv and
TOPMED_CHECKLIST.md exist; no other chromosome VCF; every output REF/ALT pair is a biallelic SNP.
lifted:  >= 200 output variants share a position with the GRCh38 fixture and >= 99 % of those have
         the same REF and ALT; liftover_summary.tsv reports lifted variants.
swapped: the report counts REF/ALT swaps; at shared positions the output REF/ALT equal the GRCh38
         fixture's, and each sample's ALT allele count equals that in the fixture (>= 99.9 %).
Pure Python (gzip), so it runs on a CI host without htslib.
"""
import gzip
import os
import sys


def read_vcf(path, want_gt=False):
    samples, rows = [], {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            gt = None
            if want_gt:
                gt = [sum(a == "1" for a in g.split(":")[0].replace("|", "/").split("/")) if "." not in g.split(":")[0] else None
                      for g in f[9:]]
            rows[int(f[1])] = (f[3], f[4], gt)
    return samples, rows


def main():
    kind, out, ref_vcf = sys.argv[1:4]
    fails = 0

    def check(ok, msg):
        nonlocal fails
        print(("PASS  " if ok else "FAIL  ") + msg)
        fails += 0 if ok else 1

    files = sorted(os.listdir(out))
    need = ["chr22.vcf.gz", "chr22.vcf.gz.tbi", "imputation_prep_report.tsv", "imputation_prep_samples.tsv", "TOPMED_CHECKLIST.md"]
    check(all(f in files for f in need), f"outputs present: {files}")
    check(not [f for f in files if f.startswith("chr") and f.endswith(".vcf.gz") and f != "chr22.vcf.gz"], "only chr22 VCF written")
    samples, got = read_vcf(os.path.join(out, "chr22.vcf.gz"), want_gt=(kind == "swapped"))
    check(len(got) > 500 and all(len(r) == 1 and len(a) == 1 for r, a, _ in got.values()),
          f"{len(got)} biallelic SNPs in chr22.vcf.gz")
    rsamples, ref = read_vcf(ref_vcf, want_gt=(kind == "swapped"))
    shared = sorted(set(got) & set(ref))
    same = [p for p in shared if got[p][:2] == ref[p][:2]]
    frac = len(same) / len(shared) if shared else 0
    check(len(shared) >= 200 and frac >= 0.99, f"{len(shared)} positions shared with GRCh38 genotypes; {100 * frac:.2f}% same REF/ALT")
    report = open(os.path.join(out, "imputation_prep_report.tsv")).read().splitlines()
    hdr = report[0].split("\t")
    total = dict(zip(hdr, report[-1].split("\t")))
    if kind == "lifted":
        lift = os.path.join(out, "liftover_summary.tsv")
        n = int(open(lift).read().splitlines()[1].split("\t")[1]) if os.path.exists(lift) else 0
        check(n > 1000, f"liftover_summary.tsv: {n} variants lifted")
    else:
        check(int(total["ref_alt_swapped"]) > 100, f"report: {total['ref_alt_swapped']} REF/ALT swaps")
        idx = [rsamples.index(s) for s in samples]
        n = agree = 0
        for p in same:
            g, r = got[p][2], ref[p][2]
            for i, j in enumerate(idx):
                if g[i] is not None and r[j] is not None:
                    n += 1
                    agree += g[i] == r[j]
        check(n > 0 and agree / n >= 0.999, f"ALT allele counts agree with the source genotypes: {agree}/{n}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

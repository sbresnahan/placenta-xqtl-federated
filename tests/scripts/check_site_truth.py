#!/usr/bin/env python3
"""check_site_truth.py EXPORT_COHORT_DIR TRUTH.tsv ANC...  -- M3 acceptance checks on a test export.

For every planted effect in TRUTH.tsv and every listed stratum:
  * the phenotype is significant (qval <= 0.05) in the permutation table B of its arm
    (expression -> 'expression'; splicing -> 'splicing_ungrouped', the per-intron arm);
  * at the planted variant, the nominal table A has p <= 1e-3 and beta with the planted sign;
  * the B lead variant is reported (it may be a variant in LD with the planted one).
Also: no stratum outside ANC... is exported, and MANIFEST.tsv lists every file with its sha256.
Exit 1 on any failure.
"""
import glob
import hashlib
import os
import sys

import pandas as pd

ARM = {"expression": "expression", "splicing": "splicing_ungrouped"}


def read(path_noext):
    for ext, f in ((".parquet", pd.read_parquet), (".tsv.gz", lambda p: pd.read_csv(p, sep="\t"))):
        if os.path.exists(path_noext + ext):
            return f(path_noext + ext)
    raise FileNotFoundError(path_noext)


def main():
    root, truth_path, ancs = sys.argv[1], sys.argv[2], sys.argv[3:]
    cohort = os.path.basename(root.rstrip("/"))
    truth = pd.read_csv(truth_path, sep="\t")
    fails = 0

    def check(ok, msg):
        nonlocal fails
        print(("PASS  " if ok else "FAIL  ") + msg)
        fails += 0 if ok else 1

    exported = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    check(exported == sorted(ancs), f"exported strata {exported} == expected {sorted(ancs)}")
    for anc in ancs:
        for t in truth.itertuples():
            arm = ARM[t.modality]
            b = read(os.path.join(root, anc, "B_mapping", f"{cohort}_{anc}_B_permutation_{arm}"))
            row = b[b["phenotype_id"] == t.phenotype_id]
            sig = len(row) == 1 and float(row["qval"].iloc[0]) <= 0.05
            lead = row["variant_id"].iloc[0] if len(row) else "NA"
            q = float(row["qval"].iloc[0]) if len(row) else float("nan")
            check(sig, f"{anc} {arm} {t.phenotype_id}: q={q:.2e} lead={lead} (planted {t.variant_id})")
            chrom = t.variant_id.split(":")[0]
            a = read(os.path.join(root, anc, "A_nominal", f"{cohort}_{anc}_A_nominal_{t.modality}_{chrom}"))
            hit = a[(a["phenotype_id"] == t.phenotype_id) & (a["variant_id"] == t.variant_id)]
            ok = len(hit) == 1 and hit["pval"].iloc[0] <= 1e-3 and (hit["beta"].iloc[0] > 0) == (t.effect_per_alt > 0)
            desc = f"beta={hit['beta'].iloc[0]:.2f} p={hit['pval'].iloc[0]:.1e}" if len(hit) else "pair missing"
            check(ok, f"{anc} A nominal at planted variant {t.variant_id}: {desc} (planted {t.effect_per_alt:+})")

    man = pd.read_csv(os.path.join(root, "MANIFEST.tsv"), sep="\t")
    files = sorted(os.path.relpath(p, os.path.dirname(root.rstrip("/")))
                   for p in glob.glob(os.path.join(root, "**", "*"), recursive=True) if os.path.isfile(p))
    files = [f for f in files if not f.endswith("MANIFEST.tsv")]
    check(sorted(man["path"]) == files, f"MANIFEST lists all {len(files)} files")
    bad = [r.path for r in man.head(20).itertuples()
           if hashlib.sha256(open(os.path.join(os.path.dirname(root.rstrip("/")), r.path), "rb").read()).hexdigest() != r.sha256]
    check(not bad, "MANIFEST sha256 matches (first 20 files)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

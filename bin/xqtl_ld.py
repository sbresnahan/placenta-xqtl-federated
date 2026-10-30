#!/usr/bin/env python3
"""xqtl_ld.py -- in-sample LD reference for one stratum (export C).

  xqtl_ld.py --geno EUR --variant-qc EUR.variant_qc.tsv --out EUR_ld \
             [--maf-min 0.05 --rsq-min 0.8 --window-mb 3 --step-mb 1]

Variants with stratum MAF >= --maf-min and imputation R2 >= --rsq-min are kept. For every
chromosome, windows of --window-mb start every --step-mb; each window with >= 2 variants gives
  OUT/<OUT>_<chrom>_<start>_<end>.npz      key 'R': signed Pearson r of ALT dosages, float16, square
  OUT/<OUT>_<chrom>_<start>_<end>.variants.tsv.gz   variant_id, chrom, pos, ref, alt, eaf (row order of R)
plus OUT/<OUT>.windows.tsv (one row per window). No sample-level values are written.
"""
import argparse
import os

import numpy as np
import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geno", required=True)
    ap.add_argument("--variant-qc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--rsq-min", type=float, default=0.8)
    ap.add_argument("--window-mb", type=float, default=3)
    ap.add_argument("--step-mb", type=float, default=1)
    a = ap.parse_args(argv)

    from tensorqtl import pgen
    vq = pd.read_csv(a.variant_qc, sep="\t", dtype={"variant_id": str, "chrom": str})
    keep = vq[(vq["maf"] >= a.maf_min) & (vq["rsq"].fillna(1.0) >= a.rsq_min)]
    os.makedirs(a.out, exist_ok=True)
    base = os.path.basename(a.out)
    pgr = pgen.PgenReader(a.geno)
    pv = pgr.pvar_df
    pv_idx = pd.Series(np.arange(len(pv)), index=pv["id"].values)
    W, S = int(a.window_mb * 1e6), int(a.step_mb * 1e6)
    rows = []
    for chrom, sub in keep.groupby("chrom", sort=False):
        sub = sub.sort_values("pos")
        ids = sub["variant_id"].tolist()
        n_miss = int(pv_idx.reindex(ids).isna().sum())
        if n_miss:
            raise SystemExit(f"xqtl_ld.py: {n_miss} variant_qc variants missing from {a.geno}.pvar")
        D = np.asarray(pgr.read_dosages_list(ids), dtype=np.float32)   # variants x samples
        D[D < 0] = np.nan
        mu = np.nanmean(D, axis=1, keepdims=True)
        D = np.where(np.isnan(D), mu, D) - mu
        sd = np.sqrt((D ** 2).sum(axis=1, keepdims=True))
        sd[sd == 0] = 1.0
        Z = D / sd
        pos = sub["pos"].to_numpy()
        start = (pos.min() // S) * S
        while start <= pos.max():
            end = start + W
            m = (pos >= start) & (pos < end)
            if m.sum() >= 2:
                R = np.clip(Z[m] @ Z[m].T, -1, 1).astype(np.float16)
                stem = os.path.join(a.out, f"{base}_{chrom}_{start}_{end}")
                np.savez_compressed(stem + ".npz", R=R)
                sub.loc[m, ["variant_id", "chrom", "pos", "ref", "alt", "af"]].rename(columns={"af": "eaf"}).to_csv(
                    stem + ".variants.tsv.gz", sep="\t", index=False, float_format="%.6g")
                rows.append({"chrom": chrom, "start": int(start), "end": int(end), "n_variants": int(m.sum()),
                             "file": os.path.basename(stem) + ".npz"})
            start += S
    pd.DataFrame(rows, columns=["chrom", "start", "end", "n_variants", "file"]).to_csv(
        os.path.join(a.out, f"{base}.windows.tsv"), sep="\t", index=False)
    print(f"{a.out}: {len(rows)} LD windows from {len(keep)} variants (MAF >= {a.maf_min}, R2 >= {a.rsq_min})")


if __name__ == "__main__":
    main()

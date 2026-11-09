#!/usr/bin/env python3
"""Synthetic GWAS for the central acceptance test.

Simulates z-scores for one trait over the variants of one exported LD window (file C, one ancestry) with two
causal variants:
  * the planted expression lead of ENSGFIX0000011 (chr22:19730097:A:G) -> should colocalize with that eQTL;
  * a variant in weak LD (|r| < 0.05) with the planted ENSGFIX0000006 lead, near that gene -> should NOT
    colocalize with the ENSGFIX0000006 eQTL (H3).
z ~ MVN(R @ lambda, R), lambda = non-centrality at the causal variants. beta/se are on a standardized-genotype
scale for N individuals. Deterministic for a given seed.

  make_gwas.py --export-dir <site export>/<COHORT>/EUR --out gwas_EUR.tsv.gz [--n 20000 --seed 7]
               [--also-export-dir <other export>/<COHORT2>/EUR ...]

--also-export-dir restricts the distinct (H3) causal variant to variants that are also in the LD files (C)
of the other exports, so it is in the fine-mapping variant set shared by all contributions.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

COLOC_LEAD = "chr22:19730097:A:G"
H3_LEAD = "chr22:19504472:T:A"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--lambda-coloc", type=float, default=9.0)
    ap.add_argument("--lambda-h3", type=float, default=8.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--also-export-dir", action="append", default=[])
    a = ap.parse_args(argv)
    win = pd.read_csv(glob.glob(os.path.join(a.export_dir, "C_ld", "*_C_ld.windows.tsv"))[0], sep="\t")
    pos = [int(COLOC_LEAD.split(":")[1]), int(H3_LEAD.split(":")[1])]
    w = win[(win["start"] <= min(pos)) & (win["end"] >= max(pos))].sort_values("n_variants").iloc[-1]
    stem = glob.glob(os.path.join(a.export_dir, "C_ld", f"*_C_ld_{w['chrom']}_{w['start']}_{w['end']}.npz"))[0]
    var = pd.read_csv(stem[:-4] + ".variants.tsv.gz", sep="\t")
    R = np.load(stem)["R"].astype(np.float64)
    np.fill_diagonal(R, 1.0)
    ids = list(var["variant_id"])
    c1, h = ids.index(COLOC_LEAD), ids.index(H3_LEAD)
    vpos = var["pos"].to_numpy()
    maf = np.minimum(var["eaf"], 1 - var["eaf"]).to_numpy()
    shared = np.ones(len(ids), dtype=bool)
    for d in a.also_export_dir:
        other = set()
        for f in glob.glob(os.path.join(d, "C_ld", "*_C_ld_*.variants.tsv*")):
            other |= set(pd.read_csv(f, sep="\t", usecols=["variant_id"])["variant_id"])
        shared &= var["variant_id"].isin(other).to_numpy()
    cand = np.where((np.abs(R[h]) < 0.05) & (np.abs(R[c1]) < 0.05) & (np.abs(vpos - vpos[h]) < 150000)
                    & (maf > 0.1) & shared)[0]
    c2 = int(cand[np.argmin(np.abs(vpos[cand] - vpos[h]))])
    lam = np.zeros(len(ids))
    lam[c1], lam[c2] = a.lambda_coloc, a.lambda_h3
    mu = R @ lam
    ev, evec = np.linalg.eigh(R)
    L = evec * np.sqrt(np.clip(ev, 0, None))
    rng = np.random.default_rng(a.seed)
    z = mu + L @ rng.standard_normal(len(ids))
    se = np.full(len(ids), 1 / np.sqrt(a.n))
    out = pd.DataFrame({"variant_id": ids, "beta": z * se, "se": se, "eaf": var["eaf"], "n": a.n})
    out.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
    print(f"causal_coloc\t{ids[c1]}\ncausal_h3\t{ids[c2]}\tr_with_{H3_LEAD}={R[h, c2]:.3f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""xqtl_covariates.py -- assemble the tensorQTL covariate table for one stratum x modality x k.

  xqtl_covariates.py --bed EUR_expression.bed.gz --samples EUR.samples.tsv --metadata metadata.tsv
      --geno-pcs EUR.geno_pcs.tsv [--cell-props cell_proportions.tsv] [--hcp EUR_expression.hcp_k10.tsv]
      [--exclude Maternal] --out EUR_expression.covariates_k10.tsv

Atlas 25_build_covariates.py at a single site (no cohort indicators):
  tiers (pruning priority)  sex < genotype PCs < cell types < HCP factors
  sex                       M -> 0, F -> 1 (replicates must agree), missing -> mode
  cell types                replicate mean; the dominant type is dropped as the compositional
                            reference; --exclude names (bare or ct_) are removed; arcsinh;
                            mean-centred; prefixed ct_
  pruning                   drop sd < 1e-8, then iteratively drop the lower-priority member of the
                            max-|r| pair while |r| > 0.9
Columns follow the phenotype BED (genotype IDs). Writes OUT and OUT.pruning.tsv.
"""
import argparse
import sys

import numpy as np
import pandas as pd

TIER = {"demographic": 0, "geno_pc": 1, "celltype": 2, "hcp": 3}


def prune_correlated(cov, priority, threshold):
    dropped = []
    cur = cov.copy()
    while cur.shape[0] > 1:
        c = np.abs(np.corrcoef(cur.values.astype(float)))
        np.fill_diagonal(c, 0.0)
        c[~np.isfinite(c)] = 0.0
        if c.max() <= threshold:
            break
        i, j = np.unravel_index(c.argmax(), c.shape)
        a, b = cur.index[i], cur.index[j]
        drop, keep = (b, a) if priority[a] <= priority[b] else (a, b)
        dropped.append({"covariate": drop, "reason": "correlated", "correlated_with": keep,
                        "r_value": round(float(c[i, j]), 4)})
        cur = cur.drop(index=drop)
    return cur, dropped


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--geno-pcs", required=True)
    p.add_argument("--cell-props")
    p.add_argument("--hcp")
    p.add_argument("--exclude", nargs="*", default=[])
    p.add_argument("--cor-threshold", type=float, default=0.9)
    p.add_argument("--min-sd", type=float, default=1e-8)
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    ids = list(pd.read_csv(a.bed, sep="\t", nrows=0).columns[4:])
    smp = pd.read_csv(a.samples, sep="\t", dtype=str)
    rmap = dict(zip(smp.rnaseq_id, smp.genotype_id))
    blocks, priority, dropped = [], {}, []

    meta = pd.read_csv(a.metadata, sep="\t", dtype={"rnaseq_id": str, "genotype_id": str})
    meta = meta[meta.rnaseq_id.isin(rmap)].copy()
    meta["gid"] = meta.rnaseq_id.map(rmap)
    sex = meta["sex"].astype(str).str.strip().map({"M": 0.0, "F": 1.0, "m": 0.0, "f": 1.0,
                                                   "Male": 0.0, "Female": 1.0})
    sex = sex.groupby(meta.gid).mean().reindex(ids)
    bad = sex[(sex > 0) & (sex < 1)]
    if len(bad):
        sys.exit(f"xqtl_covariates.py: replicate RNA-seq samples disagree on sex for {len(bad)} individual(s); fix metadata")
    if sex.isna().all():
        sex[:] = 0.0
    sex = sex.fillna(sex.mode().iloc[0])
    blocks.append(pd.DataFrame({"sex": sex}).T)
    priority["sex"] = (TIER["demographic"], 0)

    pcs = pd.read_csv(a.geno_pcs, sep="\t", dtype={"genotype_id": str}).set_index("genotype_id")
    miss = [i for i in ids if i not in pcs.index]
    if miss:
        sys.exit(f"xqtl_covariates.py: {len(miss)} phenotype columns lack genotype PCs")
    pcs = pcs.loc[ids].T
    for n in pcs.index:
        priority[n] = (TIER["geno_pc"], int(n[2:]))
    blocks.append(pcs)

    if a.cell_props:
        ct = pd.read_csv(a.cell_props, sep="\t", dtype={"sample": str})
        ct = ct[ct["sample"].isin(rmap)].copy()
        ct["gid"] = ct["sample"].map(rmap)
        ct = ct.drop(columns="sample").groupby("gid").mean().reindex(ids)
        if ct.isna().any().any():
            ct = ct.fillna(ct.median())
        dominant = ct.mean().idxmax()
        dropped.append({"covariate": "ct_" + dominant, "reason": "reference_cell_type", "correlated_with": "", "r_value": ""})
        ct = ct.drop(columns=[dominant])
        for e in a.exclude:
            name = e[3:] if e.startswith("ct_") else e
            if name in ct.columns:
                ct = ct.drop(columns=[name])
                dropped.append({"covariate": "ct_" + name, "reason": "excluded", "correlated_with": "", "r_value": ""})
        ct = np.arcsinh(ct)
        ct = (ct - ct.mean()).T
        ct.index = ["ct_" + str(c) for c in ct.index]
        for i, n in enumerate(ct.index):
            priority[n] = (TIER["celltype"], i)
        blocks.append(ct)

    if a.hcp:
        h = pd.read_csv(a.hcp, sep="\t", index_col=0)
        if h.shape[0]:
            h.columns = h.columns.astype(str)
            h = h[ids]
            for n in h.index:
                priority[n] = (TIER["hcp"], int(str(n).split("_")[1]))
            blocks.append(h)

    cov = pd.concat(blocks)[ids].astype(float)
    sd = cov.std(axis=1)
    for n in sd.index[sd < a.min_sd]:
        dropped.append({"covariate": n, "reason": "zero_variance", "correlated_with": "", "r_value": ""})
    cov = cov.loc[sd >= a.min_sd]
    cov, d2 = prune_correlated(cov, priority, a.cor_threshold)
    dropped += d2
    cov.index.name = "covariate"
    cov.to_csv(a.out, sep="\t", float_format="%.8g")
    pd.DataFrame(dropped, columns=["covariate", "reason", "correlated_with", "r_value"]).to_csv(
        a.out.replace(".tsv", "") + ".pruning.tsv", sep="\t", index=False)
    print(f"{a.out}: {cov.shape[0]} covariates x {cov.shape[1]} samples ({len(dropped)} dropped)")


if __name__ == "__main__":
    main()

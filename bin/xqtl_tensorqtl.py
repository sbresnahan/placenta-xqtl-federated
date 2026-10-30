#!/usr/bin/env python3
"""xqtl_tensorqtl.py -- tensorQTL cis mapping for one stratum x modality arm (atlas 27_run_tensorqtl.py).

Subcommands (common: --bed --geno PREFIX --cov COVARIATES.tsv [--groups FILE] --window --maf --seed):
  cis          permutation scan (map_cis, nperm), Storey q on pval_beta via compute_qvalues.R
               (BH fallback, recorded in the 'qvalue_method' column) -> OUT.cis.parquet + OUT.cis_summary.tsv
  eval         HCP-k evaluation: permutation scan on --eval-chrom when it has >= --eval-min phenotypes,
               otherwise genome-wide; writes OUT.eval.tsv (k, scope, n_tested, n_egenes)
  independent  conditional scan (map_independent) for chunk --chunk-index of the q <= fdr units;
               the FULL cis table is passed so results equal a monolithic run -> OUT.parquet
  nominal      all cis pairs (map_nominal, ungrouped) -> OUT.nominal/<OUT>.cis_qtl_pairs.<chr>.parquet

Genotypes are the stratum pgen (imputed dosages). Samples are the intersection of phenotype,
genotype and covariate columns, in phenotype order.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


def log(msg):
    print(msg, flush=True)


def load_inputs(a):
    import tensorqtl
    from tensorqtl import genotypeio
    ph, ph_pos = tensorqtl.read_phenotype_bed(a.bed)
    gdf, vdf = genotypeio.load_genotypes(a.geno, dosages=True)
    cov = pd.read_csv(a.cov, sep="\t", index_col=0, dtype={0: str})
    cov.columns = cov.columns.astype(str)
    cov = cov.T
    common = [s for s in ph.columns if s in set(gdf.columns) and s in set(cov.index)]
    if len(common) < 10:
        sys.exit(f"xqtl_tensorqtl.py: only {len(common)} samples shared by phenotypes, genotypes and covariates")
    ph = ph[common]
    gdf = gdf[common]
    cov = cov.loc[common] if cov.shape[1] else None
    group_s = None
    if a.groups and os.path.exists(a.groups) and os.path.getsize(a.groups) > 0:
        g = pd.read_csv(a.groups, sep="\t", header=None, names=["phenotype_id", "group_id"], dtype=str)
        group_s = pd.Series(g.group_id.values, index=g.phenotype_id.values)
        group_s = group_s[~group_s.index.duplicated()]
        miss = set(ph.index) - set(group_s.index)
        if miss:
            sys.exit(f"xqtl_tensorqtl.py: {len(miss)} phenotypes lack a group")
        group_s = group_s.loc[ph.index]
        pcols = [c for c in ph_pos.columns if c != "chr"]
        span = ph_pos.assign(g=group_s).groupby("g")[pcols].nunique()
        if (span > 1).any().any():
            sys.exit("xqtl_tensorqtl.py: phenotypes of a group must share the gene TSS coordinates "
                     f"(BED start/end); {int((span > 1).any(axis=1).sum())} group(s) do not")
    log(f"{ph.shape[0]} phenotypes, {gdf.shape[0]} variants, {len(common)} samples, "
        f"{0 if cov is None else cov.shape[1]} covariates, grouped={group_s is not None}")
    return ph, ph_pos, gdf, vdf, cov, group_s


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def add_qvalues(res):
    p = pd.to_numeric(res["pval_beta"], errors="coerce").to_numpy(dtype=float)
    bad = ~np.isfinite(p)
    pin = np.where(bad, 1.0, p)
    bridge = shutil.which("compute_qvalues.R") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "compute_qvalues.R")
    method = "storey"
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "p.tsv"), os.path.join(td, "q.tsv")
        pd.DataFrame({"pval_beta": pin}).to_csv(fi, sep="\t", index=False)
        pr = subprocess.run(["Rscript", bridge, fi, fo], capture_output=True, text=True)
        if pr.returncode == 0 and os.path.exists(fo):
            q = pd.read_csv(fo, sep="\t")["qval"].to_numpy(dtype=float)
            log(pr.stdout.strip())
        else:
            log("WARNING: Storey q-values failed (" + pr.stderr.strip().splitlines()[-1:][0] if pr.stderr.strip() else "WARNING: Storey q-values failed")
            log("         using Benjamini-Hochberg on pval_beta")
            q, method = bh(pin), "bh"
    q[bad] = 1.0
    res["qval"] = q
    res["qvalue_method"] = method
    return res


def run_cis(a, ph, ph_pos, gdf, vdf, cov, group_s):
    from tensorqtl import cis
    res = cis.map_cis(gdf, vdf, ph, ph_pos, covariates_df=cov, group_s=group_s, nperm=a.nperm,
                      window=a.window, maf_threshold=a.maf, seed=a.seed, verbose=False)
    res = res.copy()
    # ungrouped: index = phenotype_id; grouped: index = top phenotype_id, plus group_id/group_size
    res.index.name = "phenotype_id"
    res = res.reset_index()
    return add_qvalues(res)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["cis", "eval", "independent", "nominal"])
    ap.add_argument("--bed", required=True)
    ap.add_argument("--geno", required=True)
    ap.add_argument("--cov", required=True)
    ap.add_argument("--groups")
    ap.add_argument("--window", type=int, default=1_000_000)
    ap.add_argument("--maf", type=float, default=0.01)
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261029)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--k", default="NA")
    ap.add_argument("--eval-chrom", default="chr1")
    ap.add_argument("--eval-min", type=int, default=300)
    ap.add_argument("--cis-table")
    ap.add_argument("--chunk-index", type=int, default=1)
    ap.add_argument("--chunk-size", type=int, default=100)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    import torch
    torch.set_num_threads(max(1, a.threads))
    ph, ph_pos, gdf, vdf, cov, group_s = load_inputs(a)

    if a.command == "cis":
        res = run_cis(a, ph, ph_pos, gdf, vdf, cov, group_s)
        res.to_parquet(a.out + ".cis.parquet", index=False)
        n_sig = int((res["qval"] <= a.fdr).sum())
        pd.DataFrame([{"n_tested": len(res), "n_significant": n_sig, "fdr": a.fdr, "nperm": a.nperm,
                       "grouped": group_s is not None, "qvalue_method": res["qvalue_method"].iloc[0],
                       "n_samples": ph.shape[1],
                       "n_covariates": 0 if cov is None else cov.shape[1]}]).to_csv(
            a.out + ".cis_summary.tsv", sep="\t", index=False)
        log(f"{a.out}: {n_sig}/{len(res)} units at q <= {a.fdr}")

    elif a.command == "eval":
        chrom_ph = ph_pos.index[ph_pos["chr"].astype(str) == a.eval_chrom]
        scope = "genome"
        if len(chrom_ph) >= a.eval_min:
            ph, ph_pos = ph.loc[chrom_ph], ph_pos.loc[chrom_ph]
            if group_s is not None:
                group_s = group_s.loc[chrom_ph]
            scope = a.eval_chrom
        res = run_cis(a, ph, ph_pos, gdf, vdf, cov, group_s)
        unit = "group_id" if group_s is not None else "phenotype_id"
        n_e = res.loc[res["qval"] <= a.fdr, unit].nunique()
        pd.DataFrame([{"k": a.k, "scope": scope, "n_tested": len(res), "n_egenes": int(n_e),
                       "n_covariates": 0 if cov is None else cov.shape[1]}]).to_csv(a.out + ".eval.tsv", sep="\t", index=False)
        log(f"{a.out}: k={a.k} scope={scope} eGenes={n_e}/{len(res)}")

    elif a.command == "independent":
        from tensorqtl import cis
        full = pd.read_parquet(a.cis_table)
        unit = "group_id" if group_s is not None else "phenotype_id"
        sig = full[full["qval"] <= a.fdr]
        chunk = sig.iloc[(a.chunk_index - 1) * a.chunk_size:a.chunk_index * a.chunk_size]
        members = set(chunk[unit])
        ids = [p for p in ph.index if (group_s.get(p) if group_s is not None else p) in members]
        out = pd.DataFrame()
        if ids:
            cis_df = full.set_index("phenotype_id")  # map_cis layout: top phenotype as index, group_id column
            try:
                out = cis.map_independent(gdf, vdf, cis_df, ph.loc[ids], ph_pos.loc[ids], covariates_df=cov,
                                          group_s=group_s, maf_threshold=a.maf, fdr=a.fdr, fdr_col="qval",
                                          nperm=a.nperm, window=a.window, seed=a.seed, verbose=False)
            except ValueError as e:
                if "No significant phenotypes" not in str(e):
                    raise
                log(f"WARNING: {e}")
        out = out.reset_index(drop=True) if len(out) else out
        out.to_parquet(a.out + ".parquet", index=False)
        log(f"{a.out}: chunk {a.chunk_index} of {math.ceil(len(sig) / a.chunk_size) if len(sig) else 0}: "
            f"{len(members)} units -> {len(out)} independent signals")

    else:  # nominal
        from tensorqtl import cis
        od = a.out + ".nominal"
        os.makedirs(od, exist_ok=True)
        cis.map_nominal(gdf, vdf, ph, ph_pos, os.path.basename(a.out), covariates_df=cov,
                        maf_threshold=a.maf, window=a.window, output_dir=od, write_top=False, verbose=False)
        log(f"{a.out}: nominal pairs in {od}: {sorted(os.listdir(od))}")


if __name__ == "__main__":
    main()

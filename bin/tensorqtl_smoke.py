#!/usr/bin/env python3
"""tensorqtl_smoke.py -- infrastructure smoke test for tensorQTL inside a container.

Generates a small SYNTHETIC software-test fixture in memory (random genotype
dosages, phenotypes with one planted cis effect, random covariates), runs
tensorQTL cis.map_cis (permutations) and cis.map_nominal, and writes a one-line
PASS/FAIL result. These values are not data and carry no biological meaning;
they only check that tensorQTL, torch, and the requested device work.

Usage: tensorqtl_smoke.py --device cpu|cuda --threads N --out result.tsv
Also writes cis_perm.tsv (phenotype_id, pval_beta) for the q-value bridge test (T9q).

The tensorQTL import warning "'rfunc' cannot be imported" is expected: like the
atlas, we never use tensorQTL's rpy2 path (calculate_qvalues / logp=True); Storey
q-values come from R qvalue through compute_qvalues.R instead.
"""
import argparse
import os
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", choices=["cpu", "cuda"], default="cpu")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--test-id", default="T9_tensorqtl_cpu")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    def write(status, detail):
        with open(args.out, "w") as fh:
            fh.write(f"{args.test_id}\t{status}\t{detail}\n")

    try:
        import numpy as np
        import pandas as pd
        import torch
        import tensorqtl
        from tensorqtl import cis
    except Exception as exc:  # import failure is itself the finding
        write("FAIL", f"import error: {exc!r}")
        return 0

    torch.set_num_threads(max(1, args.threads))
    if args.device == "cuda" and not torch.cuda.is_available():
        write("FAIL", f"torch {torch.__version__}: CUDA not available (check --nv and GPU queue)")
        return 0
    device = torch.device(args.device)

    rng = np.random.default_rng(1)
    n, n_var, n_phen = 120, 3000, 40
    samples = [f"s{i}" for i in range(n)]
    pos = np.sort(rng.choice(np.arange(1, 5_000_000), n_var, replace=False))
    var_ids = [f"chr1:{p}:A:G" for p in pos]
    maf = rng.uniform(0.05, 0.5, n_var)
    G = rng.binomial(2, maf[:, None], size=(n_var, n)).astype(np.float32)
    genotype_df = pd.DataFrame(G, index=var_ids, columns=samples)
    variant_df = pd.DataFrame({"chrom": "chr1", "pos": pos}, index=var_ids)

    tss = np.sort(rng.choice(np.arange(1_000_000, 4_000_000), n_phen, replace=False))
    Y = rng.normal(size=(n_phen, n))
    # plant one cis effect so the pipeline must find something
    nearest = int(np.argmin(np.abs(pos - tss[0])))
    Y[0] += 1.5 * (G[nearest] - G[nearest].mean())
    phen_ids = [f"feat{i}" for i in range(n_phen)]
    phenotype_df = pd.DataFrame(Y, index=phen_ids, columns=samples)
    phenotype_pos_df = pd.DataFrame({"chr": "chr1", "pos": tss}, index=phen_ids)
    covariates_df = pd.DataFrame(rng.normal(size=(n, 3)), index=samples,
                                 columns=["c1", "c2", "c3"])

    t0 = time.time()
    try:
        res = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                          covariates_df=covariates_df, nperm=1000, seed=1,
                          window=1_000_000, verbose=False)
        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
                        prefix="smoke", covariates_df=covariates_df,
                        output_dir=".", window=1_000_000, verbose=False)
    except Exception as exc:
        write("FAIL", f"tensorQTL error on {args.device}: {exc!r}")
        return 0
    dt = time.time() - t0
    # pval_beta table for the Storey q-value bridge (T9q, compute_qvalues.R in the R image),
    # exactly as 27_run_tensorqtl.py hands it over in the atlas
    res[["pval_beta"]].rename_axis("phenotype_id").to_csv("cis_perm.tsv", sep="\t")
    top_p = float(res.loc["feat0", "pval_beta"])
    ok = top_p < 1e-3
    detail = (f"tensorqtl={getattr(tensorqtl, '__version__', 'unknown')} torch={torch.__version__} "
              f"device={device} threads={args.threads} runtime_s={dt:.1f} "
              f"planted_feature_pval_beta={top_p:.2e} host={os.uname().nodename}")
    write("PASS" if ok else "FAIL", detail)
    return 0


if __name__ == "__main__":
    sys.exit(main())

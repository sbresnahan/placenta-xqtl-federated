#!/usr/bin/env python3
"""xqtl_select_k.py -- choose the number of HCP factors (k*) from the eval tables and write settings.

  xqtl_select_k.py --anc EUR --modality expression --arm expression [--arm expression_ungrouped]
                   --norm-summary EUR_expression.norm_summary.tsv --nperm 10000 --window 1000000
                   --maf 0.01 --fdr 0.05 --seed 20261029 [--grouped] --out-k k.txt EVAL_TSV...

k* = the k with the most eGenes (unique phenotypes or groups at q <= fdr in the eval scope); ties go
to the smaller k (atlas rule). The same k* is used for the grouped and ungrouped arm of a modality.
Writes {ANC}_{arm}.settings.json for each --arm and the bare k* to --out-k.
"""
import argparse
import json

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anc", required=True)
    ap.add_argument("--modality", required=True)
    ap.add_argument("--arm", action="append", required=True)
    ap.add_argument("--norm-summary")
    ap.add_argument("--n-samples", type=int)
    ap.add_argument("--nperm", type=int, required=True)
    ap.add_argument("--window", type=int, required=True)
    ap.add_argument("--maf", type=float, required=True)
    ap.add_argument("--fdr", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--grouped", action="store_true", help="the first --arm is the grouped arm")
    ap.add_argument("--out-k", required=True)
    ap.add_argument("evals", nargs="+")
    a = ap.parse_args(argv)

    ev = pd.concat([pd.read_csv(p, sep="\t") for p in a.evals], ignore_index=True)
    ev["k"] = ev["k"].astype(int)
    ev = ev.sort_values("k").drop_duplicates("k")
    best = ev.sort_values(["n_egenes", "k"], ascending=[False, True]).iloc[0]
    k = int(best["k"])

    method, n = None, a.n_samples
    if a.norm_summary:
        ns = pd.read_csv(a.norm_summary, sep="\t").iloc[0]
        method, n = str(ns["method"]), int(ns["individuals"])
    grid = [{"k": int(r.k), "n_egenes": int(r.n_egenes), "n_tested": int(r.n_tested), "scope": str(r.scope),
             "n_covariates": int(r.n_covariates)} for r in ev.itertuples()]
    for i, arm in enumerate(a.arm):
        st = {"modality": a.modality, "arm": arm, "grouped": bool(a.grouped and i == 0),
              "normalization": method, "hcp_grid": grid, "k_selected": k, "k_rule": "argmax eGenes, ties to smaller k",
              "n_samples": n, "nperm": a.nperm, "cis_window": a.window, "maf_threshold": a.maf,
              "fdr": a.fdr, "seed": a.seed}
        with open(f"{a.anc}_{arm}.settings.json", "w") as fh:
            json.dump(st, fh, indent=1)
    with open(a.out_k, "w") as fh:
        fh.write(f"{k}\n")
    print(f"{a.anc} {a.modality}: k* = {k} ({int(best['n_egenes'])} eGenes)")


if __name__ == "__main__":
    main()

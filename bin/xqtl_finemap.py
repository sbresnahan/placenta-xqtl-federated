#!/usr/bin/env sushie-python
"""Multi-ancestry fine-mapping of central xQTL loci with SuShiE in summary-statistics mode.

Runs in the SuShiE environment (`sushie-python`: jax, numpy 2, polars; no pandas). Each input .npz is
written by `xqtl_central.py finemap-prep` and holds, for one feature, the shared variant list, one
row per ancestry of within-ancestry meta z-scores, the N-weighted LD matrix (file C) and N.

  xqtl_finemap.py --loci L1.npz L2.npz ... --L 10 --out-prefix shard1
writes <prefix>.pips.tsv (variant PIPs, credible-set membership), <prefix>.cs.tsv (credible sets),
<prefix>.weights.tsv (posterior effect per ancestry = SuShiE prediction weights, standardized
genotype and phenotype scale), <prefix>.corr.tsv (estimated cross-ancestry effect correlation per
single effect) and <prefix>.status.tsv (one row per locus).
"""
import argparse
import csv
import os
import sys
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")


def fit(npz, L, purity, max_iter, seed):
    from sushie.infer_ss import infer_sushie_ss
    z = np.load(npz)
    zs = [np.asarray(x, dtype=np.float64) for x in z["zs"]]
    lds = [np.asarray(x, dtype=np.float64) for x in z["lds"]]
    ns = np.asarray(z["ns"], dtype=np.float64)
    p = len(z["variants"])
    res = infer_sushie_ss(lds=lds, zs=zs, ns=ns, L=min(L, p), purity=purity, max_iter=max_iter,
                          min_snps=min(100, p), max_select=min(250, p), seed=seed)
    return z, res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", nargs="+", required=True)
    ap.add_argument("--L", type=int, default=10)
    ap.add_argument("--purity", type=float, default=0.5)
    ap.add_argument("--max-iter", type=int, default=500)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    out = {k: open(f"{a.out_prefix}.{k}.tsv", "w", newline="") for k in ("pips", "cs", "weights", "corr", "status")}
    w = {k: csv.writer(fh, delimiter="\t", lineterminator="\n") for k, fh in out.items()}
    w["pips"].writerow(["locus_id", "variant_id", "pip", "pip_cs", "cs_index"])
    w["cs"].writerow(["locus_id", "cs_index", "variant_id", "alpha", "c_alpha"])
    w["weights"].writerow(["locus_id", "variant_id", "ancestry", "weight"])
    w["corr"].writerow(["locus_id", "effect_index", "ancestry_1", "ancestry_2", "rho"])
    w["status"].writerow(["locus_id", "status", "n_variants", "ancestries", "n_cs", "elbo_iterations", "seconds", "detail"])
    for npz in a.loci:
        lid = os.path.basename(npz)[:-len(".npz")]
        t0 = time.time()
        try:
            z, r = fit(npz, a.L, a.purity, a.max_iter, a.seed)
        except Exception as e:  # noqa: BLE001  (recorded per locus; other loci continue)
            w["status"].writerow([lid, "error", "", "", "", "", f"{time.time() - t0:.1f}", str(e)[:300]])
            print(f"{lid}: {e}", file=sys.stderr)
            continue
        vids = [str(v) for v in z["variants"]]
        ancs = [str(x) for x in z["ancestries"]]
        pip_all = np.asarray(r.pip_all)
        pip_cs = np.asarray(r.pip_cs)
        cs = r.cs  # polars DataFrame: CSIndex, SNPIndex, alpha, c_alpha, pip_all, pip_cs
        cs_of = {}
        rows = cs.rows(named=True) if hasattr(cs, "rows") else []
        for row in rows:
            ci, si = int(row["CSIndex"]), int(row["SNPIndex"])
            cs_of.setdefault(si, ci)
            w["cs"].writerow([lid, ci, vids[si], f"{row['alpha']:.6g}", f"{row['c_alpha']:.6g}"])
        for i, v in enumerate(vids):
            w["pips"].writerow([lid, v, f"{pip_all[i]:.6g}", f"{pip_cs[i]:.6g}", cs_of.get(i, "")])
        wt = np.asarray(r.posteriors.post_mean).sum(axis=0)       # p x k
        for i, v in enumerate(vids):
            for k in range(len(ancs)):
                w["weights"].writerow([lid, v, ancs[k], f"{wt[i, k]:.6g}"])
        cov = np.asarray(r.priors.effect_covar)                     # L x k x k
        n_cs = len({int(x["CSIndex"]) for x in rows})
        for li in range(cov.shape[0]):
            for k1 in range(len(ancs)):
                for k2 in range(k1 + 1, len(ancs)):
                    den = np.sqrt(cov[li, k1, k1] * cov[li, k2, k2])
                    w["corr"].writerow([lid, li + 1, ancs[k1], ancs[k2], f"{cov[li, k1, k2] / den:.4f}" if den > 0 else "NA"])
        w["status"].writerow([lid, "ok" if n_cs else "no_cs", len(vids), ";".join(ancs), n_cs,
                              len(np.atleast_1d(np.asarray(r.elbo))),
                              f"{time.time() - t0:.1f}", ""])
        print(f"{lid}: {len(vids)} variants, {len(ancs)} ancestries, {n_cs} credible sets", file=sys.stderr)
    for fh in out.values():
        fh.close()


if __name__ == "__main__":
    main()

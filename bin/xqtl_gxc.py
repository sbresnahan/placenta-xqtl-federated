#!/usr/bin/env python3
"""xqtl_gxc.py -- genotype x cell-type interaction at lead xQTL signals (atlas 49 port).

  xqtl_gxc.py --anc EUR --modality expression --pairs EUR_expression.gx_pairs.tsv \
      --bed EUR_expression.bed.gz --cov EUR_expression.covariates.tsv --pfile EUR \
      --samples EUR.samples.tsv --cell-props cell_proportions.tsv [--exclude Maternal] \
      [--min-n 20] [--maternal-threshold 0.10] [--fdr 0.05] [--workers 1] --out EUR_expression

Model per pair x target cell type (atlas 06/49 _fit_gxc), OLS with a t-test:
    Y ~ 1 + g0 + c0 + g0:c0 + Z
  Y   the stratum's normalized phenotype (BED columns = genotype IDs)
  g0  ALT dosage from the stratum pgen, centred over the analysed samples
  c0  arcsinh(cell proportion), centred over the stratum, then re-centred over the analysed samples
  Z   the modality's final mapping covariates without the target's ct_ covariate; covariates
      constant in the analysed subset are dropped; rows with any non-finite value are removed
Cell proportions are per RNA-seq sample and are averaged over technical replicates to genotype
IDs. Targets: every cell type except --exclude (default Maternal). Maternal sensitivity: the same
fit after setting Y to missing for samples with maternal fraction > --maternal-threshold.
Status: ok | no_common_samples | insufficient_n | monomorphic | no_celltype_variation |
rank_deficient | insufficient_df | invalid_se (atlas). BH q within stratum x modality over ok tests.

Writes OUT.gxc.parquet (one row per pair x cell type; untested pairs one row with cell_type empty)
and OUT.gxc_annotation.tsv.gz (one row per pair: best cell type by (q, p); primary_cell_type is
the best cell type if q <= --fdr, else 'unassigned').
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

STATE = None


def norm_name(v):
    return "".join(ch.lower() for ch in str(v) if ch.isalnum())


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def fit_gxc(y, g, cell, cov, min_n=20):
    """Atlas _fit_gxc: y, g, cell are Series on sample IDs; cov a samples x covariates frame."""
    from scipy import stats
    common = y.index.intersection(g.index).intersection(cell.index)
    if cov is not None and len(cov.columns):
        common = common.intersection(cov.index)
    if len(common) == 0:
        return {"status": "no_common_samples"}
    yv = pd.to_numeric(y.loc[common], errors="coerce").to_numpy(dtype=float)
    gv = pd.to_numeric(g.loc[common], errors="coerce").to_numpy(dtype=float)
    cv = pd.to_numeric(cell.loc[common], errors="coerce").to_numpy(dtype=float)
    if cov is not None and len(cov.columns):
        zdf = cov.loc[common].apply(pd.to_numeric, errors="coerce")
        keep = [c for c in zdf.columns if np.nanstd(zdf[c].to_numpy(dtype=float)) > 1e-12]
        z = zdf[keep].to_numpy(dtype=float) if keep else np.empty((len(common), 0))
    else:
        z = np.empty((len(common), 0))
    finite = np.isfinite(yv) & np.isfinite(gv) & np.isfinite(cv)
    if z.shape[1]:
        finite &= np.isfinite(z).all(axis=1)
    yv, gv, cv, z = yv[finite], gv[finite], cv[finite], z[finite]
    n = len(yv)
    if n < min_n:
        return {"status": "insufficient_n", "n": n}
    if np.std(gv) <= 1e-12:
        return {"status": "monomorphic", "n": n}
    if np.std(cv) <= 1e-12:
        return {"status": "no_celltype_variation", "n": n}
    af = float(np.mean(gv) / 2.0)
    maf = min(af, 1.0 - af)
    g0 = gv - gv.mean()
    c0 = cv - cv.mean()
    X = np.column_stack([np.ones(n), g0, c0, g0 * c0, z])
    p = X.shape[1]
    beta, _, rank, _ = np.linalg.lstsq(X, yv, rcond=None)
    rank = int(rank)
    if rank < p:
        return {"status": "rank_deficient", "n": n, "af": af, "maf": maf, "design_rank": rank, "design_cols": p}
    df = n - p
    if df <= 2:
        return {"status": "insufficient_df", "n": n, "af": af, "maf": maf, "design_rank": rank, "design_cols": p}
    resid = yv - X @ beta
    sigma2 = float(resid @ resid / df)
    try:
        var_b = sigma2 * np.linalg.inv(X.T @ X)[3, 3]
    except np.linalg.LinAlgError:       # numerically singular although lstsq reports full rank
        return {"status": "rank_deficient", "n": n, "af": af, "maf": maf, "design_rank": rank, "design_cols": p}
    if not np.isfinite(var_b) or var_b <= 0:
        return {"status": "invalid_se", "n": n, "af": af, "maf": maf}
    se = float(np.sqrt(var_b))
    t = float(beta[3] / se)
    q25, q50, q75 = np.quantile(cv, [0.25, 0.5, 0.75])
    cm = float(cv.mean())
    return {"status": "ok", "n": n, "af": af, "maf": maf, "design_rank": rank, "design_cols": p,
            "beta_G": float(beta[1]), "beta_C": float(beta[2]), "beta_GxC": float(beta[3]),
            "se_GxC": se, "t_GxC": t, "p_GxC": float(2 * stats.t.sf(abs(t), df)), "df": df,
            "cell_transformed_q25": float(q25), "cell_transformed_q50": float(q50),
            "cell_transformed_q75": float(q75), "cell_transformed_iqr": float(q75 - q25),
            "beta_GxC_iqr": float(beta[3] * (q75 - q25)),
            "beta_G_at_cell_q25": float(beta[1] + beta[3] * (q25 - cm)),
            "beta_G_at_cell_q50": float(beta[1] + beta[3] * (q50 - cm)),
            "beta_G_at_cell_q75": float(beta[1] + beta[3] * (q75 - cm))}


def drop_target(cov, cell_type):
    target = norm_name("ct_" + str(cell_type))
    return cov.drop(columns=[c for c in cov.columns if norm_name(c) == target], errors="ignore")


def fit_pair(pair):
    s = STATE
    pid, vid = pair["phenotype_id"], pair["variant_id"]
    y = s["phen"].loc[pid]
    g = s["dos"][vid]
    y_s = y.where(s["maternal_keep"].reindex(y.index).fillna(False)) if s["maternal_keep"] is not None else None
    rows = []
    for ct in s["targets"]:
        fit = fit_gxc(y, g, s["deconv_t"][ct], s["cov_by_ct"][ct], s["min_n"])
        q25, q50, q75 = s["raw_q"][ct]
        rec = {**pair, "cell_type": ct, "cell_prop_raw_q25": q25, "cell_prop_raw_q50": q50,
               "cell_prop_raw_q75": q75, **fit,
               "maternal_sensitivity_n": np.nan, "maternal_sensitivity_beta_GxC": np.nan,
               "maternal_sensitivity_se_GxC": np.nan, "maternal_sensitivity_p_GxC": np.nan,
               "maternal_sensitivity_status": "no_maternal_fraction"}
        if y_s is not None:
            sens = fit_gxc(y_s, g, s["deconv_t"][ct], s["cov_by_ct"][ct], s["min_n"])
            rec.update({"maternal_sensitivity_status": sens.get("status"),
                        "maternal_sensitivity_n": sens.get("n", np.nan),
                        "maternal_sensitivity_beta_GxC": sens.get("beta_GxC", np.nan),
                        "maternal_sensitivity_se_GxC": sens.get("se_GxC", np.nan),
                        "maternal_sensitivity_p_GxC": sens.get("p_GxC", np.nan)})
        rows.append(rec)
    return rows


def run(pairs, phen, dos, cell_raw, cov, exclude, min_n=20, maternal_threshold=0.10, fdr=0.05, workers=1):
    """pairs: gx_pairs frame; phen: phenotypes x samples; dos: samples x variants; cell_raw: samples x
    cell types (raw proportions, stratum samples); cov: samples x covariates. Returns (tests, annotation)."""
    global STATE
    excl = {norm_name(e[3:] if str(e).startswith("ct_") else e) for e in exclude}
    targets = [c for c in cell_raw.columns if norm_name(c) not in excl]
    mat = next((c for c in cell_raw.columns if "maternal" in norm_name(c)), None)
    deconv_t = np.arcsinh(cell_raw[targets])
    deconv_t = deconv_t - deconv_t.mean(axis=0)
    STATE = {"phen": phen, "dos": dos, "targets": targets, "deconv_t": deconv_t, "min_n": min_n,
             "cov_by_ct": {ct: drop_target(cov, ct) for ct in targets},
             "raw_q": {ct: tuple(float(x) for x in cell_raw[ct].quantile([0.25, 0.5, 0.75])) for ct in targets},
             "maternal_keep": (cell_raw[mat] <= maternal_threshold) if mat else None}
    pairs = pairs.copy()
    pairs["untested_reason"] = pairs["untested_reason"].fillna("").astype(str)
    miss_v = (pairs["untested_reason"] == "") & ~pairs["variant_id"].isin(dos.columns)
    pairs.loc[miss_v, "untested_reason"] = "variant_absent"
    miss_p = (pairs["untested_reason"] == "") & ~pairs["phenotype_id"].isin(phen.index)
    pairs.loc[miss_p, "untested_reason"] = "phenotype_absent"
    todo = pairs[pairs["untested_reason"] == ""].to_dict("records")
    rows = []
    if workers > 1 and len(todo) > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(processes=min(workers, len(todo))) as pool:
            for r in pool.imap(fit_pair, todo, chunksize=max(1, len(todo) // (workers * 8))):
                rows.extend(r)
    else:
        for pr in todo:
            rows.extend(fit_pair(pr))
    STATE = None
    for pr in pairs[pairs["untested_reason"] != ""].to_dict("records"):
        rows.append({**pr, "cell_type": "", "status": pr["untested_reason"]})
    tests = pd.DataFrame(rows)
    for c in ("n", "beta_GxC", "se_GxC", "p_GxC"):
        if c not in tests.columns:
            tests[c] = np.nan
    if len(tests) == 0:
        tests = pd.DataFrame(columns=list(pairs.columns) + ["cell_type", "status", "p_GxC"])
    tests["untested_reason"] = np.where(tests["status"] == "ok", "", tests["status"].astype(str))
    tests["q_GxC"] = np.nan
    ok = (tests["status"] == "ok") & tests["p_GxC"].notna()
    if ok.any():
        tests.loc[ok, "q_GxC"] = bh(tests.loc[ok, "p_GxC"].to_numpy(float))
    tests = tests.sort_values(["phenotype_id", "variant_id", "cell_type"]).reset_index(drop=True)

    ann = []
    keys = ["ancestry", "modality", "phenotype_id", "variant_id", "source"]
    for key, sub in tests.groupby(keys, sort=True):
        base = dict(zip(keys, key))
        valid = sub[(sub["status"] == "ok") & sub["p_GxC"].notna()].sort_values(["q_GxC", "p_GxC"])
        if valid.empty:
            ann.append({**base, "primary_cell_type": "unassigned", "best_cell_type": "", "n_celltypes_tested": 0,
                        "untested_reason": sub["untested_reason"].iloc[0] or "no_ok_model"})
            continue
        b = valid.iloc[0]
        q = b["q_GxC"]
        ann.append({**base, "primary_cell_type": b["cell_type"] if q <= fdr else "unassigned",
                    "best_cell_type": b["cell_type"], "primary_beta_GxC": b["beta_GxC"],
                    "primary_se_GxC": b["se_GxC"], "primary_p_GxC": b["p_GxC"], "primary_q_GxC": q,
                    "n_celltypes_tested": len(valid),
                    "maternal_sensitivity_beta_GxC": b["maternal_sensitivity_beta_GxC"],
                    "maternal_sensitivity_p_GxC": b["maternal_sensitivity_p_GxC"], "untested_reason": ""})
    return tests, pd.DataFrame(ann)


def load_cell_props(path, samples_tsv, ids):
    smp = pd.read_csv(samples_tsv, sep="\t", dtype=str)
    rmap = dict(zip(smp.rnaseq_id, smp.genotype_id))
    ct = pd.read_csv(path, sep="\t", dtype={"sample": str})
    ct = ct[ct["sample"].isin(rmap)].copy()
    ct["gid"] = ct["sample"].map(rmap)
    num = ct.drop(columns="sample").set_index("gid").apply(pd.to_numeric, errors="coerce")
    return num.groupby(level=0).mean().reindex(ids)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--anc", "--modality", "--pairs", "--bed", "--cov", "--pfile", "--samples", "--cell-props", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--exclude", nargs="*", default=["Maternal"])
    p.add_argument("--min-n", type=int, default=20)
    p.add_argument("--maternal-threshold", type=float, default=0.10)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--workers", type=int, default=1)
    a = p.parse_args(argv)

    pairs = pd.read_csv(a.pairs, sep="\t", dtype={"untested_reason": str}, keep_default_na=False,
                        na_values={"site_qval": [""]})
    want = set(pairs.loc[pairs["untested_reason"] == "", "phenotype_id"])
    bed = pd.read_csv(a.bed, sep="\t", dtype={3: str})
    ids = [str(c) for c in bed.columns[4:]]
    phen = bed[bed.iloc[:, 3].isin(want)].set_index(bed.columns[3]).iloc[:, 3:]
    phen.columns = ids
    cov = pd.read_csv(a.cov, sep="\t", index_col=0)
    cov.columns = cov.columns.astype(str)
    cov = cov.T.apply(pd.to_numeric, errors="coerce")
    cell = load_cell_props(a.cell_props, a.samples, ids)

    variants = sorted(set(pairs.loc[pairs["untested_reason"] == "", "variant_id"]))
    if variants:
        from tensorqtl.pgen import PgenReader
        pgr = PgenReader(a.pfile, select_samples=ids)
        present = [v for v in variants if v in set(pgr.pvar_df["id"])]
        dos = pd.DataFrame(pgr.read_dosages_list(present).T if present else np.empty((len(ids), 0)),
                           index=ids, columns=present) if present else pd.DataFrame(index=ids)
        dos = dos.replace(-9, np.nan)
    else:
        dos = pd.DataFrame(index=ids)

    tests, ann = run(pairs, phen, dos, cell, cov, a.exclude, a.min_n, a.maternal_threshold, a.fdr, a.workers)
    tests.to_parquet(a.out + ".gxc.parquet", index=False)
    ann.to_csv(a.out + ".gxc_annotation.tsv.gz", sep="\t", index=False, compression="gzip")
    n_ok = int((tests["status"] == "ok").sum())
    n_ass = int((ann["primary_cell_type"] != "unassigned").sum()) if len(ann) else 0
    print(f"{a.out}: {len(pairs)} pairs, {n_ok} ok G x C tests, {n_ass} pairs with a primary cell type")
    if len(pairs) and n_ok == 0 and (pairs["untested_reason"] == "").any():
        print("xqtl_gxc.py: WARNING no G x C model could be fitted", file=sys.stderr)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    main()

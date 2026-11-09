"""Central meta-analysis (bin/xqtl_central.py) and its validator checks."""
import json

import numpy as np
import pandas as pd
import pytest
from conftest import load
from scipy import stats

cen = load("xqtl_central")
val = load("xqtl_validate")


def hand_ivw(b, s):
    """Textbook fixed-effect IVW and DerSimonian-Laird random effects for one group."""
    b, s = np.asarray(b, float), np.asarray(s, float)
    w = 1 / s ** 2
    bf = np.sum(w * b) / np.sum(w)
    sf = np.sqrt(1 / np.sum(w))
    q = np.sum(w * (b - bf) ** 2)
    k = len(b)
    tau2 = max(0.0, (q - (k - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    wr = 1 / (s ** 2 + tau2)
    br = np.sum(wr * b) / np.sum(wr)
    sr = np.sqrt(1 / np.sum(wr))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return dict(beta_fe=bf, se_fe=sf, q_het=q, tau2=tau2, beta_re=br, se_re=sr, i2=i2,
                p_het=stats.chi2.sf(q, k - 1), p_fe=2 * stats.norm.sf(abs(bf / sf)))


def test_ivw_matches_hand_calculation():
    rng = np.random.default_rng(1)
    groups = [(rng.normal(0.2, 0.15, k), rng.uniform(0.05, 0.2, k)) for k in (2, 3, 5)]
    gid = np.concatenate([np.full(len(b), i) for i, (b, _) in enumerate(groups)])
    beta = np.concatenate([b for b, _ in groups])
    se = np.concatenate([s for _, s in groups])
    res = cen.ivw(gid, beta, se, len(groups))
    for i, (b, s) in enumerate(groups):
        h = hand_ivw(b, s)
        for key, v in h.items():
            assert res[key][i] == pytest.approx(v, rel=1e-9, abs=1e-12), key


def test_ivw_single_contribution():
    res = cen.ivw(np.array([0]), np.array([0.3]), np.array([0.1]), 1)
    assert res["beta_fe"][0] == pytest.approx(0.3) and res["se_fe"][0] == pytest.approx(0.1)
    assert res["tau2"][0] == 0 and np.isnan(res["i2"][0]) and np.isnan(res["p_het"][0])
    assert res["beta_re"][0] == pytest.approx(0.3)


def test_meta_frame_mask_n_eaf():
    df = pd.DataFrame({"pid": ["a", "a", "b"], "beta": [0.1, 0.3, 0.5], "se": [0.1, 0.1, 0.2],
                       "n": [100, 300, 50], "eaf": [0.2, 0.4, 0.5], "contrib": [0, 2, 1],
                       "freq_outlier": [False, True, False]})
    m = cen.meta_frame(df, ["pid"]).set_index("pid")
    assert m.loc["a", "contrib_mask"] == 0b101 and m.loc["b", "contrib_mask"] == 0b010
    assert m.loc["a", "n"] == 400 and m.loc["a", "eaf"] == pytest.approx((100 * 0.2 + 300 * 0.4) / 400)
    assert m.loc["a", "n_freq_outlier"] == 1
    assert m.loc["a", "beta_fe"] == pytest.approx(0.2)


def test_acat():
    assert cen.acat([0.03]) == 0.03
    assert np.isnan(cen.acat([np.nan]))
    p = np.array([0.01, 0.2, 0.7])
    t = np.mean(np.tan((0.5 - p) * np.pi))
    assert cen.acat(p) == pytest.approx(0.5 - np.arctan(t) / np.pi)
    assert cen.acat([1e-20, 0.5]) == pytest.approx(2e-20, rel=1e-6)   # small-p approximation


def test_bh():
    q = cen.bh([0.01, 0.04, 0.03, np.nan])
    assert np.allclose(q[:3], [0.03, 0.04, 0.04]) and np.isnan(q[3])


def _d(rows):
    return pd.DataFrame(rows, columns=["contrib", "ancestry", "variant_id", "chrom", "pos", "ref", "alt", "eaf", "n"])


def test_harmonize_swap_flip_palindromic():
    rows = []
    for c, (r1, a1, e1), (r2, a2, e2), (r3, a3, e3), (r4, a4, e4), (r5, a5, e5) in [
        (0, ("A", "G", 0.3), ("A", "C", 0.2), ("A", "T", 0.45), ("A", "T", 0.45), ("C", "G", 0.1)),
        (1, ("G", "A", 0.7), ("T", "G", 0.2), ("T", "A", 0.55), ("A", "T", 0.45), ("G", "C", 0.9)),
    ]:
        for pos, (r, a, e) in zip((100, 200, 300, 400, 500), ((r1, a1, e1), (r2, a2, e2), (r3, a3, e3),
                                                               (r4, a4, e4), (r5, a5, e5))):
            rows.append([c, "EUR", f"chr1:{pos}:{r}:{a}", "chr1", pos, r, a, e, 100])
    h = cen.harmonize(_d(rows), 0.4, 0.2).set_index(["contrib", "variant_id"])
    get = lambda c, v: h.loc[(c, v)]  # noqa: E731
    assert get(1, "chr1:100:G:A")["action"] == "swap" and get(1, "chr1:100:G:A")["sign"] == -1
    assert get(1, "chr1:100:G:A")["canonical_id"] == "chr1:100:A:G"
    assert get(1, "chr1:100:G:A")["eaf_h"] == pytest.approx(0.3)
    assert get(1, "chr1:200:T:G")["action"] == "flip" and get(1, "chr1:200:T:G")["sign"] == 1
    # palindromic, orientations differ, MAF > 0.4: dropped everywhere
    assert get(0, "chr1:300:A:T")["action"] == "drop_palindromic"
    assert get(1, "chr1:300:T:A")["action"] == "drop_palindromic"
    # palindromic, same orientation everywhere: kept (forward-strand GRCh38 IDs)
    assert get(1, "chr1:400:A:T")["action"] == "keep"
    # palindromic, orientations differ, low MAF: read literally as a swap
    assert get(1, "chr1:500:G:C")["action"] == "swap" and get(1, "chr1:500:G:C")["sign"] == -1
    assert not h["freq_outlier"].any()


def test_harmonize_frequency_outlier():
    rows = [[c, "EUR", "chr1:100:A:G", "chr1", 100, "A", "G", e, 100] for c, e in ((0, 0.2), (1, 0.21), (2, 0.22), (3, 0.6))]
    h = cen.harmonize(_d(rows), 0.4, 0.2).set_index("contrib")
    assert bool(h.loc[3, "freq_outlier"]) and not h.loc[[0, 1, 2], "freq_outlier"].any()


def test_exposure_pooled_sd_and_rescale(tmp_path):
    rng = np.random.default_rng(3)
    xs = [rng.normal(39, 1.2, 50), rng.normal(38, 1.6, 80)]
    con = []
    for c, x in enumerate(xs):
        pre = tmp_path / f"s{c}" / "EUR" / f"C{c}_EUR"
        pre.parent.mkdir(parents=True)
        site_sd = float(np.std(x, ddof=1))
        pd.DataFrame([
            {"exposure": "GA", "ancestry": "EUR", "n": len(x), "mean_raw": x.mean(), "sd_raw": np.std(x),
             "scale_mean": x.mean(), "scale_sd": site_sd, "type": "continuous", "transform": "none",
             "status": "available"},
            {"exposure": "GA", "ancestry": "pooled", "n": len(x), "mean_raw": x.mean(), "sd_raw": np.std(x),
             "scale_mean": x.mean(), "scale_sd": site_sd, "type": "continuous", "transform": "none",
             "status": "available"},
        ]).to_csv(str(pre) + "_H_exposures.tsv", sep="\t", index=False)
        con.append({"contrib": c, "prefix": str(pre)})
    scales, tab = cen.exposure_scale(pd.DataFrame(con))
    allx = np.concatenate(xs)
    assert scales["GA"]["pooled_sd"] == pytest.approx(np.std(allx), rel=1e-12)
    assert scales["GA"]["pooled_mean"] == pytest.approx(allx.mean())
    d = pd.DataFrame({"b_gi": [0.5], "b_gi_se": [0.1]})
    d2, f = cen.rescale(d.copy(), 1, "GA", scales)
    assert f == pytest.approx(np.std(allx) / np.std(xs[1], ddof=1))
    assert d2["b_gi"][0] == pytest.approx(0.5 * f) and d2["b_gi_se"][0] == pytest.approx(0.1 * f)
    # transformed exposures are not rescaled
    scales["GA"]["rescale"] = False
    assert cen.rescale(d.copy(), 1, "GA", scales)[1] == 1.0


def test_concat_skips_empty(tmp_path):
    a, b, e = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "e.tsv"
    a.write_text("x\ty\n1\tNA\n")
    b.write_text("x\ty\n2\t3\n")
    e.write_text("")
    out = tmp_path / "o.tsv"
    cen.main(["concat", "--out", str(out), str(a), str(b), str(e)])
    assert out.read_text() == "x\ty\n1\tNA\n2\t3\n"


def _study(tmp_path, name, manifest=True, strata=True):
    d = tmp_path / name
    d.mkdir()
    if manifest:
        (d / "MANIFEST.tsv").write_text("path\tsha256\n")
    if strata:
        (d / f"{name}_cohort_I_strata.tsv").write_text("cohort\tancestry\n")
    return d


def _run_central(tmp_path, central):
    p = tmp_path / "params.json"
    p.write_text(json.dumps({"cohort_id": "MDACC", "central": central}))
    out_j, out_r = tmp_path / "r.json", tmp_path / "r.tsv"
    rc = val.main(["--params", str(p), "--schema", val_schema(), "--mode", "central",
                   "--out-json", str(out_j), "--out-report", str(out_r)])
    return rc, pd.read_csv(out_r, sep="\t")


def val_schema():
    from conftest import SCHEMA
    return SCHEMA


def test_validate_central_ok_and_relative_paths(tmp_path):
    _study(tmp_path, "S1")
    _study(tmp_path, "S2")
    (tmp_path / "studies.tsv").write_text("study_id\texport_dir\nA\tS1\nB\tS2\n")
    (tmp_path / "g.tsv.gz").write_bytes(b"")
    (tmp_path / "gwas.tsv").write_text("trait_id\tpath\tancestry\tn\ttrait_type\tprop_cases\n"
                                       "T1\tg.tsv.gz\tEUR\t1000\tcc\t0.3\n")
    rc, rep = _run_central(tmp_path, {"studies": str(tmp_path / "studies.tsv"),
                                      "gwas_catalog": str(tmp_path / "gwas.tsv")})
    assert rc == 0, rep[rep.status == "FAIL"].to_string()
    assert (rep.check == "central.studies.A").any()


def test_validate_central_failures(tmp_path):
    _study(tmp_path, "S1", manifest=False)
    (tmp_path / "studies.tsv").write_text("study_id\texport_dir\nA\tS1\nA\tmissing\n")
    (tmp_path / "gwas.tsv").write_text("trait_id\tpath\tancestry\tn\ttrait_type\tprop_cases\n"
                                       "T1\tx/g.tsv\tEUR\t1000\tcc\t\nT2\ty/g.tsv\tEUR\t1000\tbinary\t\n")
    rc, rep = _run_central(tmp_path, {"studies": str(tmp_path / "studies.tsv"),
                                      "gwas_catalog": str(tmp_path / "gwas.tsv")})
    assert rc == 1
    fails = " ".join(rep.loc[rep.status == "FAIL", "detail"])
    for msg in ("duplicate study_id", "no MANIFEST.tsv", "export_dir not found", "basenames must be unique",
                "trait_type must be", "needs prop_cases", "missing file"):
        assert msg in fails, msg


def test_validate_central_requires_studies(tmp_path):
    rc, rep = _run_central(tmp_path, {})
    assert rc == 1 and (rep.check == "required.central.studies").any()

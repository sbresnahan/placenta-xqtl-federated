"""Unit tests for bin/xqtl_gxe_core.py and bin/xqtl_gxe.py (G x E, atlas 07_gxe_mapping port).

The interaction math is checked against tensorQTL itself (calculate_interaction_nominal,
filter_maf_interaction), as in the atlas test_gxe_scan.py; those tests need torch + tensorQTL
(present in the image) and are skipped where they are not installed (CI unit job).
"""
import json
import types

import numpy as np
import pandas as pd
import pytest

from conftest import load

X = load("xqtl_gxe")

N_SAMPLES = 200
N_VARIANTS = 400


def simulate_genotypes(n_samples, n_variants=N_VARIANTS, rho=0.85, seed=0):
    """Markov-chain haplotypes -> diploid dosages with decaying LD (samples x variants)."""
    rng = np.random.default_rng(seed)
    mafs = rng.uniform(0.1, 0.4, n_variants)
    haps = np.zeros((n_samples * 2, n_variants), dtype=np.int8)
    haps[:, 0] = rng.binomial(1, mafs[0], n_samples * 2)
    for j in range(1, n_variants):
        flip = rng.binomial(1, 1 - rho, n_samples * 2)
        draw = rng.binomial(1, mafs[j], n_samples * 2)
        haps[:, j] = np.where(flip, draw, haps[:, j - 1])
    return (haps[0::2] + haps[1::2]).astype(np.float64)


@pytest.fixture(scope="module")
def tq():
    torch = pytest.importorskip("torch")
    pytest.importorskip("tensorqtl")
    core = load("xqtl_gxe_core")
    return types.SimpleNamespace(torch=torch, core=core, cpu=torch.device("cpu"))


@pytest.fixture(scope="module")
def inputs():
    rng = np.random.default_rng(7)
    G = simulate_genotypes(N_SAMPLES, seed=3)
    e = rng.normal(size=N_SAMPLES)
    cov = rng.normal(size=(N_SAMPLES, 3))
    return G, e, cov, rng


def core_scan(tq, G, e, cov, P):
    torch, core = tq.torch, tq.core
    Q_t, dof_base = core.make_residualizer(cov, device=tq.cpu)
    dof = dof_base - 2
    g0, e0, ge0 = core.prepare_interaction_terms(torch.tensor(G.T, dtype=torch.float32),
                                                 torch.tensor(e, dtype=torch.float32), Q_t)
    Xt, Xinv, _ = core.build_window_design(g0, e0, ge0)
    out = []
    for j in range(P.shape[0]):
        p0 = core.prepare_phenotype(torch.tensor(P[j], dtype=torch.float32), Q_t)
        out.append([x.numpy() for x in core.window_fit(Xt, Xinv, p0, dof)])
    b, b_se, t = (np.stack([o[i] for o in out]) for i in range(3))
    return b, b_se, t, dof, Xt, Xinv, Q_t


# ---------------------------------------------------------------------------------- tensorQTL parity
def test_interaction_maf_filter_matches_tensorqtl(tq, inputs):
    from tensorqtl.core import filter_maf_interaction
    G, e, _, _ = inputs
    Gv = G.T.astype(np.float32)
    ours = X.interaction_maf_mask(Gv, e, 0.05)
    order = np.argsort(e, kind="mergesort")
    upper = np.ones(len(e), dtype=bool)
    upper[order[:len(e) // 2]] = False
    _, mask_t = filter_maf_interaction(tq.torch.tensor(Gv), interaction_mask_t=tq.torch.tensor(upper),
                                       maf_threshold_interaction=0.05)
    assert np.array_equal(ours, mask_t.numpy().astype(bool))


def test_interaction_maf_filter_rejects_exposure_imbalance():
    e = np.arange(20, dtype=float)
    G = np.zeros((2, 20), dtype=np.float32)
    G[0, :10] = 1.0            # MAF 0.25 overall; absent from the upper exposure half
    G[1, [0, 10]] = 1.0        # MAF 0.05 in both halves
    assert X.interaction_maf_mask(G, e, 0.05).tolist() == [False, True]
    assert X.interaction_maf_mask(G, e, 0.0).tolist() == [True, True]


def test_nominal_matches_tensorqtl(tq, inputs):
    from tensorqtl.core import Residualizer, calculate_interaction_nominal
    G, e, cov, rng = inputs
    P = rng.normal(size=(4, N_SAMPLES))
    P[0] += 0.6 * G[:, 50] + 0.8 * G[:, 50] * e
    P[1] += 0.6 * G[:, 100]
    torch = tq.torch
    res = Residualizer(torch.tensor(cov, dtype=torch.float32))
    t_ref, b_ref, bse_ref, *_ = calculate_interaction_nominal(
        torch.tensor(G.T, dtype=torch.float32), torch.tensor(P, dtype=torch.float32),
        torch.tensor(e.reshape(-1, 1), dtype=torch.float32), residualizer=res)
    b, b_se, t, dof, *_ = core_scan(tq, G, e, cov, P)
    assert dof == res.dof - 2
    for name, mine, ref in [("b", b, b_ref), ("b_se", b_se, bse_ref), ("tstat", t, t_ref)]:
        ref = np.transpose(ref.numpy(), (2, 0, 1))
        assert np.isfinite(mine).all()
        rel = np.max(np.abs(mine - ref) / (np.abs(ref) + 1e-6))
        assert rel < 5e-3, f"{name}: max relative difference {rel:.2e}"


def test_permutation_matches_brute_force(tq, inputs):
    G, e, cov, rng = inputs
    torch, core = tq.torch, tq.core
    y = 0.5 * G[:, 50] + rng.normal(size=N_SAMPLES)
    _, _, t, dof, Xt, Xinv, Q_t = core_scan(tq, G, e, cov, y[None, :])
    r2_max = float(np.max(core.t_to_r2(t[0, :, 2], dof)))
    p_raw = torch.tensor(y, dtype=torch.float32)
    r2_perm, nperm, _ = core.adaptive_permutation_scan(Xt, Xinv, p_raw, Q_t, dof, r2_max,
                                                       np.random.default_rng(99), blocks=(25,), stop_p=2.0)
    rng2 = np.random.default_rng(99)
    bf = []
    for _ in range(25):
        p0 = core.prepare_phenotype(p_raw[torch.from_numpy(rng2.permutation(N_SAMPLES))], Q_t)
        _, _, tp = core.window_fit(Xt, Xinv, p0, dof)
        bf.append(np.nanmax(core.t_to_r2(tp[:, 2].numpy(), dof)))
    assert nperm == 25
    assert np.allclose(r2_perm, bf, rtol=1e-4, atol=1e-6)


def test_interaction_recovery_and_null(tq, inputs):
    G, e, cov, rng = inputs
    y_int = 0.5 * G[:, 50] + 0.9 * G[:, 50] * e + rng.normal(size=N_SAMPLES)
    y_null = 0.8 * G[:, 50] + rng.normal(size=N_SAMPLES)
    b, _, t, *_ = core_scan(tq, G, e, cov, np.stack([y_int, y_null]))
    top = int(np.argmax(np.abs(t[0, :, 2])))
    assert abs(b[0, top, 2] - 0.9) < 0.25 and abs(t[0, top, 2]) > 4
    assert np.max(np.abs(t[1, :, 2])) < 4


def test_adaptive_schedule(tq, inputs):
    G, e, cov, rng = inputs
    core, torch = tq.core, tq.torch
    y_int = 0.5 * G[:, 50] + 1.2 * G[:, 50] * e + rng.normal(size=N_SAMPLES)
    y_null = rng.normal(size=N_SAMPLES)
    _, _, t, dof, Xt, Xinv, Q_t = core_scan(tq, G, e, cov, np.stack([y_int, y_null]))
    res = {}
    for k, y in enumerate([y_int, y_null]):
        r2_max = float(np.max(core.t_to_r2(t[k, :, 2], dof)))
        res[k] = core.adaptive_permutation_scan(Xt, Xinv, torch.tensor(y, dtype=torch.float32), Q_t, dof,
                                                r2_max, np.random.default_rng(1), blocks=(100, 400), stop_p=0.10)
    _, nperm_int, p_int = res[0]
    _, nperm_null, p_null = res[1]
    assert nperm_int == 500 and p_int == pytest.approx(1 / 501)
    assert p_null > p_int
    if nperm_null == 100:
        assert p_null > 0.10


def test_beta_approx_wrapper(tq):
    rng = np.random.default_rng(0)
    out = tq.core.beta_approx_pval(np.sort(rng.beta(0.5, 50, 1000)), r2_nominal=0.9, dof=150)
    assert 0 <= out["pval_beta"] < 0.01
    bad = tq.core.beta_approx_pval(np.full(100, 0.01), 0.5, 150)
    assert not np.isfinite(bad["pval_beta"]) or 0 <= bad["pval_beta"] <= 1


def test_ivw_meta(tq):
    out = tq.core.ivw_meta([0.5, 0.7], [0.1, 0.2])
    assert out["beta"] == pytest.approx((100 * 0.5 + 25 * 0.7) / 125)
    assert out["se"] == pytest.approx(np.sqrt(1 / 125))
    assert out["k"] == 2 and 0 <= out["p_het"] <= 1
    assert tq.core.ivw_meta([np.nan, np.nan], [0.1, 0.2])["k"] == 0


# ---------------------------------------------------------------------------------- exposures
def exposure_inputs():
    smp = {"EUR": pd.DataFrame({"rnaseq_id": ["R1", "R1b", "R2", "R3", "R4"],
                                "genotype_id": ["G1", "G1", "G2", "G3", "G4"]}),
           "EAS": pd.DataFrame({"rnaseq_id": ["R5", "R6", "R7"], "genotype_id": ["G5", "G6", "G7"]})}
    meta = pd.DataFrame({"rnaseq_id": ["R1", "R1b", "R2", "R3", "R4", "R5", "R6", "R7"],
                         "GA": [38.0, 40.0, 37.0, 50.0, 41.0, 36.0, 39.0, 40.0],     # R3 out of range
                         "bmi": [20.0, 20.0, 30.0, 25.0, 40.0, 22.0, 0.0, 28.0],      # R6 not log-able
                         "gdm": [0, 1, 1, 0, 0, 1, 0, 1]})                             # G1 replicates disagree
    return meta, smp


def test_build_exposures_scaling_and_status():
    meta, smp = exposure_inputs()
    cfg = {"min_n": 3, "min_cases": 2,
           "exposures": {"GA": {"enabled": True, "min": 20, "max": 45},
                         "bmi": {"enabled": True, "transform": "log"},
                         "gdm": {"enabled": True, "type": "binary", "min": 0, "max": 1},
                         "ogtt": {"enabled": True},
                         "off": {"enabled": False}}}
    table, summ = X.build_exposures(meta, smp, cfg)
    assert list(table.columns) == ["GA", "bmi", "gdm"]
    # GA: replicate mean (G1 = 39), G3 out of range -> missing; pooled z with ddof = 0
    raw = pd.Series({"G1": 39.0, "G2": 37.0, "G4": 41.0, "G5": 36.0, "G6": 39.0, "G7": 40.0})
    z = (raw - raw.mean()) / raw.std(ddof=0)
    assert np.isnan(table.loc["G3", "GA"])
    assert np.allclose(table.loc[z.index, "GA"], z)
    # log transform on the raw scale; non-positive -> missing
    assert np.isnan(table.loc["G6", "bmi"])
    lb = np.log(pd.Series({"G1": 20.0, "G2": 30.0, "G3": 25.0, "G4": 40.0, "G5": 22.0, "G7": 28.0}))
    assert np.allclose(table.loc[lb.index, "bmi"], (lb - lb.mean()) / lb.std(ddof=0))
    # binary: unscaled 0/1; replicate disagreement (mean 0.5) -> missing
    assert np.isnan(table.loc["G1", "gdm"]) and set(table["gdm"].dropna()) == {0.0, 1.0}
    s = summ.set_index(["ancestry", "exposure"])
    assert s.loc[("EUR", "GA"), "status"] == "available" and s.loc[("EUR", "GA"), "n"] == 3
    assert s.loc[("EAS", "GA"), "status"] == "available"
    assert s.loc[("pooled", "GA"), "n"] == 6
    assert s.loc[("EUR", "gdm"), "status"] == "insufficient_cases"       # 1 case, 2 controls
    assert s.loc[("EUR", "ogtt"), "status"] == "exposure_unavailable"
    assert s.loc[("pooled", "GA"), "scale_sd"] == pytest.approx(raw.std(ddof=0))
    cfg["min_n"] = 4
    _, summ = X.build_exposures(meta, smp, cfg)
    assert summ.set_index(["ancestry", "exposure"]).loc[("EUR", "GA"), "status"] == "insufficient_n"


# ---------------------------------------------------------------------------------- scan (CLI, fake genotypes)
@pytest.fixture()
def scan_inputs(tmp_path, monkeypatch):
    rng = np.random.default_rng(11)
    n, nv = 150, 120
    ids = [f"G{i:03d}" for i in range(n)]
    G = simulate_genotypes(n, n_variants=nv, seed=5)
    pos = 1_000_000 + 1000 * np.arange(nv)
    vids = np.array([f"chr1:{p}:A:G" for p in pos])
    e = rng.normal(size=n)
    cov = pd.DataFrame(rng.normal(size=(3, n)), index=["PC1", "PC2", "sex"], columns=ids)
    cov.index.name = "covariate"
    Y = rng.normal(size=(3, n))
    Y[0] += 0.5 * G[:, 60] + 1.0 * G[:, 60] * e                        # planted G x E
    Y[1] += 0.8 * G[:, 30]                                               # cis only
    bed = pd.DataFrame({"#chr": "chr1", "start": [pos[60], pos[30], pos[90]], "end": [pos[60] + 1, pos[30] + 1, pos[90] + 1],
                        "phenotype_id": ["P_int", "P_cis", "P_null"]})
    bed = pd.concat([bed, pd.DataFrame(Y, columns=ids)], axis=1)
    bed.to_csv(tmp_path / "b.bed.gz", sep="\t", index=False)
    cov.to_csv(tmp_path / "cov.tsv", sep="\t")
    pd.DataFrame({"GA": e}, index=pd.Index(ids, name="genotype_id")).to_csv(tmp_path / "exp.tsv", sep="\t")
    pd.DataFrame([{"ancestry": "EUR", "exposure": "GA", "n": n, "status": "available"}]).to_csv(
        tmp_path / "summ.tsv", sep="\t", index=False)
    gxp = pd.DataFrame([
        {"phenotype_id": "P_int", "site_lead_variant": vids[60], "mdacc_lead_variants": f"{vids[60]},chr1:5:C:T", "untested_reason": ""},
        {"phenotype_id": "P_cis", "site_lead_variant": vids[30], "mdacc_lead_variants": "", "untested_reason": ""},
        {"phenotype_id": "P_gone", "site_lead_variant": "", "mdacc_lead_variants": "chr1:7:A:C", "untested_reason": "phenotype_absent"}])
    gxp.to_csv(tmp_path / "gxp.tsv", sep="\t", index=False)

    def fake_load(pfile, chrom, samples, dosages=False):
        assert samples == ids
        present = set(vids) | {"chr1:5:C:T"}          # chr1:5 exists but fails the MAF filters
        return G.T.astype(np.float32), vids, pos, present
    monkeypatch.setattr(X, "load_genotypes", fake_load)
    common = ["--anc", "EUR", "--modality", "expression", "--exposure", "GA", "--chrom", "chr1",
              "--bed", str(tmp_path / "b.bed.gz"), "--cov", str(tmp_path / "cov.tsv"), "--pfile", "x",
              "--exposures", str(tmp_path / "exp.tsv"), "--summary", str(tmp_path / "summ.tsv"),
              "--window", "20000", "--blocks", "100", "400"]
    return types.SimpleNamespace(tmp=tmp_path, common=common, vids=vids, gxp=tmp_path / "gxp.tsv")


def test_scan_tier1_and_tier2(tq, scan_inputs):
    s = scan_inputs
    t1 = str(s.tmp / "t1")
    X.main(["scan", "--tier", "1", *s.common, "--out", t1])
    r1 = pd.read_parquet(t1 + ".parquet").set_index("phenotype_id")
    assert list(r1.reset_index().columns) == X.T1_COLS
    assert r1.loc["P_int", "variant_id"] == s.vids[60] and r1.loc["P_int", "b_gi"] > 0.7
    assert r1.loc["P_int", "nperm"] == 500 and r1.loc["P_int", "pval_beta"] < 1e-3
    assert r1.loc["P_int", "num_var"] == 41                                 # +-20 kb at 1 kb spacing
    st = pd.read_csv(t1 + ".status.tsv", sep="\t").iloc[0]
    assert st["status"] == "DONE" and st["n"] == 150 and st["dof"] == 150 - 3 - 4
    # the per-phenotype seed depends on the BED row only: same result on a rerun
    X.main(["scan", "--tier", "1", *s.common, "--out", t1 + "b"])
    pd.testing.assert_frame_equal(pd.read_parquet(t1 + ".parquet"), pd.read_parquet(t1 + "b.parquet"))

    t2 = str(s.tmp / "t2")
    X.main(["scan", "--tier", "2", *s.common, "--phenotypes", str(s.gxp), "--out", t2])
    r2 = pd.read_parquet(t2 + ".parquet")
    assert list(r2.columns) == X.T2_COLS
    assert set(r2["phenotype_id"]) == {"P_int", "P_cis"}
    lead = r2[(r2.phenotype_id == "P_int") & (r2.variant_id == s.vids[60])].iloc[0]
    assert lead.site_lead_variant and lead.mdacc_lead_variant
    assert r2.site_lead_variant.sum() == 2 and r2.mdacc_lead_variant.sum() == 1
    unt = pd.read_csv(t2 + ".untested.tsv", sep="\t", keep_default_na=False)
    assert unt.values.tolist() == [["P_int", "chr1:5:C:T", "interaction_maf"]]
    X.main(["scan", "--tier", "2", *s.common, "--phenotypes", str(s.gxp), "--tier2-variants", "pairs", "--out", t2 + "p"])
    rp = pd.read_parquet(t2 + "p.parquet")
    assert sorted(rp["variant_id"]) == sorted([s.vids[30], s.vids[60]])
    pd.testing.assert_frame_equal(rp.reset_index(drop=True),
                                  r2[r2.site_lead_variant | r2.mdacc_lead_variant].reset_index(drop=True))

    X.main(["merge", "--tier", "2", "--anc", "EUR", "--modality", "expression", "--exposure", "GA",
            "--shards", t2 + ".parquet", "--phenotypes", str(s.gxp), "--out", str(s.tmp / "m")])
    m = pd.read_parquet(s.tmp / "m.gxe_tier2.parquet")
    top = m[m.top_variant].set_index("phenotype_id")
    assert top.loc["P_int", "variant_id"] == s.vids[60] and top.loc["P_int", "q_top_bh"] < 1e-6
    mu = pd.read_csv(s.tmp / "m.gxe_tier2_untested.tsv", sep="\t", keep_default_na=False)
    assert ["P_gone", "", "phenotype_absent"] in mu.values.tolist()
    ms = pd.read_csv(s.tmp / "m.gxe_tier2.status.tsv", sep="\t").iloc[0]
    assert ms["status"] == "DONE" and ms["n_tested"] == 2


def test_scan_skips_and_min_n(tq, scan_inputs):
    s = scan_inputs
    out = str(s.tmp / "small")
    X.main(["scan", "--tier", "2", *s.common, "--phenotypes", str(s.gxp), "--min-n", "500", "--out", out])
    assert pd.read_csv(out + ".status.tsv", sep="\t").iloc[0]["reason"] == "insufficient_n"
    assert set(pd.read_csv(out + ".untested.tsv", sep="\t")["untested_reason"]) == {"insufficient_n"}
    pd.DataFrame([{"ancestry": "EUR", "exposure": "GA", "n": 9, "status": "insufficient_n"}]).to_csv(
        s.tmp / "summ.tsv", sep="\t", index=False)
    out = str(s.tmp / "skip")
    X.main(["scan", "--tier", "1", *s.common, "--out", out])
    st = pd.read_csv(out + ".status.tsv", sep="\t").iloc[0]
    assert st["status"] == "SKIPPED" and st["reason"] == "insufficient_n"
    # an empty tier-1 merge never calls the q-value bridge
    X.main(["merge", "--tier", "1", "--anc", "EUR", "--modality", "expression", "--exposure", "GA",
            "--shards", out + ".parquet", "--out", str(s.tmp / "skipm")])
    ms = pd.read_csv(s.tmp / "skipm.gxe_tier1.status.tsv", sep="\t").iloc[0]
    assert ms["status"] == "SKIPPED" and ms["reason"] == "insufficient_n" and ms["n_sig"] == 0
    assert len(pd.read_parquet(s.tmp / "skipm.gxe_tier1.parquet")) == 0


# ---------------------------------------------------------------------------------- sensitivity status
def test_sensitivity_off_and_blocked(tmp_path):
    meta, _ = exposure_inputs()
    meta.to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    args = ["sensitivity", "--anc", "EUR", "--modality", "expression", "--exposure", "GA", "--tier1", "t1.parquet",
            "--bed", "b", "--cov", "c", "--pfile", "p", "--exposures", "e", "--metadata", str(tmp_path / "meta.tsv"),
            "--samples", "s", "--out", str(tmp_path / "o")]
    cfg = {"exposures": {"GA": {"column": "GA"}},
           "sensitivity": {"maternal_age": {"enabled": False}, "GA": {"enabled": True, "column": "GA"}}}
    (tmp_path / "g.json").write_text(json.dumps(cfg))
    X.main(args + ["--config", str(tmp_path / "g.json")])        # GA is never its own sensitivity covariate
    st = pd.read_csv(tmp_path / "o.gxe_sensitivity.status.tsv", sep="\t").iloc[0]
    assert st["status"] == "OFF"
    assert list(pd.read_csv(tmp_path / "o.gxe_sensitivity.tsv", sep="\t").columns) == X.SENS_COLS
    cfg["sensitivity"]["parity"] = {"enabled": True}
    (tmp_path / "g.json").write_text(json.dumps(cfg))
    X.main(args + ["--config", str(tmp_path / "g.json")])
    st = pd.read_csv(tmp_path / "o.gxe_sensitivity.status.tsv", sep="\t").iloc[0]
    assert st["status"] == "BLOCKED" and "parity" in st["reason"]

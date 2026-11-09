"""Unit tests for bin/xqtl_signals.py (signal test sets) and bin/xqtl_gxc.py (G x C)."""
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from conftest import load

S = load("xqtl_signals")
G = load("xqtl_gxc")


def sim(n=80, seed=1):
    rng = np.random.default_rng(seed)
    ids = [f"I{i:03d}" for i in range(n)]
    g = rng.binomial(2, 0.35, n).astype(float)
    cell = np.arcsinh(rng.dirichlet([10, 5, 3], n))
    cell = cell - cell.mean(0)
    z = rng.normal(size=(n, 3))
    y = 0.5 * g + 0.3 * cell[:, 0] + 2.0 * (g - g.mean()) * (cell[:, 0] - cell[:, 0].mean()) \
        + z @ np.array([0.2, -0.1, 0.3]) + rng.normal(0, 1, n)
    return ids, g, cell, z, y


def test_fit_gxc_matches_r_lm(tmp_path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript not available")
    ids, g, cell, z, y = sim()
    cov = pd.DataFrame(z, index=ids, columns=["PC1", "PC2", "ct_B"])
    fit = G.fit_gxc(pd.Series(y, ids), pd.Series(g, ids), pd.Series(cell[:, 0], ids), cov)
    assert fit["status"] == "ok" and fit["n"] == 80
    d = pd.DataFrame({"y": y, "g": g, "c": cell[:, 0], "PC1": z[:, 0], "PC2": z[:, 1], "ctB": z[:, 2]})
    d.to_csv(tmp_path / "d.tsv", sep="\t", index=False)
    code = ("d <- read.delim('%s'); d$g0 <- d$g - mean(d$g); d$c0 <- d$c - mean(d$c);"
            "s <- coef(summary(lm(y ~ g0 * c0 + PC1 + PC2 + ctB, data = d)));"
            "cat(sprintf('%%.10g', c(s['g0', 1], s['c0', 1], s['g0:c0', 1:4])), sep = '\\n')") % (tmp_path / "d.tsv")
    out = subprocess.run([rscript, "-e", code], capture_output=True, text=True, check=True).stdout.split()
    bG, bC, b, se, t, p = map(float, out)
    assert fit["beta_G"] == pytest.approx(bG, rel=1e-6)
    assert fit["beta_C"] == pytest.approx(bC, rel=1e-6)
    assert fit["beta_GxC"] == pytest.approx(b, rel=1e-6)
    assert fit["se_GxC"] == pytest.approx(se, rel=1e-6)
    assert fit["p_GxC"] == pytest.approx(p, rel=1e-5)


def test_fit_gxc_statuses():
    ids, g, cell, z, y = sim(30)
    Y, C = pd.Series(y, ids), pd.Series(cell[:, 0], ids)
    cov = pd.DataFrame(z, index=ids)
    assert G.fit_gxc(Y, pd.Series(g, ids), C, cov, min_n=40)["status"] == "insufficient_n"
    assert G.fit_gxc(Y, pd.Series(np.ones(30), ids), C, cov)["status"] == "monomorphic"
    assert G.fit_gxc(Y, pd.Series(g, ids), pd.Series(np.zeros(30), ids), cov)["status"] == "no_celltype_variation"
    dup = pd.DataFrame({"a": z[:, 0], "b": 2 * z[:, 0]}, index=ids)
    assert G.fit_gxc(Y, pd.Series(g, ids), C, dup)["status"] == "rank_deficient"
    many = pd.DataFrame(np.random.default_rng(0).normal(size=(30, 24)), index=ids)
    assert G.fit_gxc(Y, pd.Series(g, ids), C, many)["status"] == "insufficient_df"
    assert G.fit_gxc(Y, pd.Series(g, [f"X{i}" for i in range(30)]), C, cov)["status"] == "no_common_samples"
    # a covariate that is constant in the subset is dropped rather than making the design singular
    cst = pd.DataFrame({"k": np.ones(30), "a": z[:, 0]}, index=ids)
    assert G.fit_gxc(Y, pd.Series(g, ids), C, cst)["status"] == "ok"


def test_drop_target_covariate():
    cov = pd.DataFrame(np.zeros((3, 3)), columns=["ct_Cytotrophoblast", "ct_Stromal", "PC1"])
    assert list(G.drop_target(cov, "Cytotrophoblast").columns) == ["ct_Stromal", "PC1"]


def test_run_assigns_primary_and_keeps_untested():
    ids, g, cell, z, y = sim(120, seed=3)
    phen = pd.DataFrame([y, np.random.default_rng(9).normal(size=120)], index=["P1", "P2"], columns=ids)
    dos = pd.DataFrame({"chr1:1:A:G": g}, index=ids)
    raw = pd.DataFrame(np.sinh(cell + 0.3), index=ids, columns=["Cyto", "Stromal", "Maternal"])
    cov = pd.DataFrame(z, index=ids, columns=["PC1", "PC2", "ct_Stromal"])
    pairs = pd.DataFrame({"ancestry": "EUR", "modality": "expression",
                          "phenotype_id": ["P1", "P2", "P3", "P1"],
                          "variant_id": ["chr1:1:A:G", "chr1:1:A:G", "chr1:1:A:G", "chr1:9:A:C"],
                          "source": ["both", "mdacc", "mdacc", "mdacc"],
                          "untested_reason": ["", "", "phenotype_absent", ""]})
    tests, ann = G.run(pairs, phen, dos, raw, cov, exclude=["Maternal"])
    assert set(tests.loc[tests.status == "ok", "cell_type"]) == {"Cyto", "Stromal"}   # Maternal excluded
    a = ann.set_index(["phenotype_id", "variant_id"])
    assert a.loc[("P1", "chr1:1:A:G"), "primary_cell_type"] == "Cyto"
    assert a.loc[("P2", "chr1:1:A:G"), "primary_cell_type"] == "unassigned"
    assert a.loc[("P3", "chr1:1:A:G"), "untested_reason"] == "phenotype_absent"
    assert a.loc[("P1", "chr1:9:A:C"), "untested_reason"] == "variant_absent"
    ok = tests[tests.status == "ok"]
    assert ok["q_GxC"].notna().all() and tests.loc[tests.status != "ok", "q_GxC"].isna().all()
    assert (ok["maternal_sensitivity_n"] <= ok["n"]).all()


def test_signal_sets_flags():
    cis = pd.DataFrame({"phenotype_id": ["A", "B", "C"], "variant_id": ["v1", "v2", "v3"], "qval": [0.01, 0.02, 0.5]})
    sig = pd.DataFrame({"phenotype_id": ["A", "B", "D", "A", "E", "F"],
                        "variant_id": ["v1", "v9", "v4", "v5", "v6", "v7"],
                        "modality": "expression", "ancestry": ["EUR", "EAS", "EAS", "EUR", "EUR", "EUR"],
                        "signal_rank": [1, 1, 1, 2, 1, 1]})
    pairs, phen = S.build("EUR", "expression", cis, ["A", "B", "C", "D", "F"], {"v1", "v2", "v3", "v4", "v9"}, sig)
    p = pairs.set_index(["phenotype_id", "variant_id"])
    assert p.loc[("A", "v1"), "source"] == "both"
    assert p.loc[("B", "v2"), "source"] == "site" and p.loc[("B", "v9"), "source"] == "mdacc"
    assert p.loc[("B", "v9"), "mdacc_ancestry"] == "EAS"
    assert p.loc[("D", "v4"), "source"] == "mdacc" and p.loc[("D", "v4"), "untested_reason"] == ""
    assert ("A", "v5") not in p.index                                     # rank 2: not a G x C pair
    assert ("C", "v3") not in p.index                                     # not significant at the site
    assert p.loc[("E", "v6"), "untested_reason"] == "phenotype_absent"
    assert p.loc[("F", "v7"), "untested_reason"] == "variant_absent"
    ph = phen.set_index("phenotype_id")
    assert set(ph.index) == {"A", "B", "D", "E", "F"}
    assert ph.loc["A", "source"] == "both" and ph.loc["A", "mdacc_min_rank"] == 1
    assert ph.loc["E", "untested_reason"] == "phenotype_absent"


def test_ir_containment_mapping():
    site = ["G1__IR_chr1_100_900_+", "G1__IR_chr1_150_400_+", "G1__IR_chr1_100_900_-", "G2__IR_chr1_100_900_+"]
    m = S.map_ir(["G1__IR_G1_chr1_200_350_+_1", "G1__IR_G1_chr1_120_800_+_2", "G1__IR_G1_chr1_50_300_+_3",
                  "G3__IR_G3_chr1_200_300_+_1"], site)
    assert m["G1__IR_G1_chr1_200_350_+_1"] == "G1__IR_chr1_150_400_+"   # shortest containing intron
    assert m["G1__IR_G1_chr1_120_800_+_2"] == "G1__IR_chr1_100_900_+"
    assert m["G1__IR_G1_chr1_50_300_+_3"] is None                       # not contained
    assert m["G3__IR_G3_chr1_200_300_+_1"] is None                       # other gene
    sig = pd.DataFrame({"phenotype_id": ["G1__IR_G1_chr1_200_350_+_1", "G3__IR_G3_chr1_200_300_+_1"],
                        "variant_id": ["v1", "v1"], "modality": "intron_retention", "ancestry": "EUR",
                        "signal_rank": 1})
    empty = pd.DataFrame(columns=["phenotype_id", "variant_id", "qval"])
    pairs, _ = S.build("EUR", "intron_retention", empty, site, {"v1"}, sig)
    p = pairs.set_index("mdacc_phenotype_id")
    assert p.loc["G1__IR_G1_chr1_200_350_+_1", "phenotype_id"] == "G1__IR_chr1_150_400_+"
    assert p.loc["G1__IR_G1_chr1_200_350_+_1", "untested_reason"] == ""
    assert p.loc["G3__IR_G3_chr1_200_300_+_1", "untested_reason"] == "phenotype_absent"


def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.2, 0.5])
    expect = np.array([0.05, 0.0666667, 0.0666667, 0.25, 0.5])
    assert np.allclose(G.bh(p), expect, atol=1e-6)


if __name__ == "__main__":
    pytest.main([os.path.abspath(__file__), "-q"])

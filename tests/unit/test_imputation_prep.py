"""Unit tests for bin/xqtl_impprep.py (sex check, liftover summary, report / build stop, checklist)."""
import gzip
import json

import numpy as np
import pandas as pd
import pytest

from conftest import load

IP = load("xqtl_impprep")


def _gt_lines(n_sites, sexes, rng):
    """chrX GT rows: males hemizygous-as-homozygous, females in HWE."""
    lines = []
    for _ in range(n_sites):
        p = rng.uniform(0.2, 0.5)
        row = []
        for s in sexes:
            if s == "M":
                a = int(rng.random() < p)
                row.append(f"{a}/{a}")
            else:
                row.append(f"{int(rng.random() < p)}/{int(rng.random() < p)}")
        lines.append("\t".join(row) + "\t\n")
    return lines


def test_sex_f_separates_males_and_females():
    rng = np.random.default_rng(1)
    sexes = ["M"] * 10 + ["F"] * 10
    F, sex, sites = IP.sex_f(_gt_lines(400, sexes, rng), len(sexes))
    assert (F[:10] > 0.8).all() and (sex[:10] == "M").all()
    assert (F[10:] < 0.2).all() and (sex[10:] == "F").all()
    assert (sites == 400).all()


def test_sex_f_no_call_below_min_sites_and_missing_gt():
    rng = np.random.default_rng(2)
    lines = _gt_lines(50, ["M", "F", "F"], rng)
    _, sex, sites = IP.sex_f(lines, 3, min_sites=100)
    assert list(sex) == ["no_call"] * 3
    lines = ["./.\t0/1\t0/0\n"] * 5
    _, _, sites = IP.sex_f(lines, 3, min_sites=1)
    assert sites[0] == 0 and sites[1] == 5


def _write(path, text):
    path.write_text(text)
    return str(path)


def test_liftover_summary_counts_reject_reasons(tmp_path):
    rej = tmp_path / "rej.vcf.gz"
    with gzip.open(rej, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        fh.write("22\t1\t.\tA\tG\t.\tNoTarget\t.\n22\t2\t.\tA\tG\t.\tMismatchedRefAllele\t.\n")
        fh.write("22\t3\t.\tA\tG\t.\tNoTarget\t.\n")
    out = tmp_path / "l.tsv"
    IP.main(["liftover", "--rejected", str(rej), "--lifted", "97", "--out", str(out)])
    d = pd.read_csv(out, sep="\t").set_index("reason")["n"].to_dict()
    assert d == {"lifted": 97, "MismatchedRefAllele": 1, "NoTarget": 2}


def _report_inputs(tmp_path, ref_match, pass_hwe=1000, source="array", lift_mm=None):
    counts = tmp_path / "chr22.counts.tsv"
    rows = [("input", 1200), ("biallelic_snps", 1100), ("pass_missing_maf", 1050), ("pass_hwe", pass_hwe),
            ("ref_alt_swapped", 30), ("ref_match", ref_match)]
    counts.write_text("chrom\tstep\tn\n" + "".join(f"chr22\t{s}\t{n}\n" for s, n in rows))
    smiss = _write(tmp_path / "s.smiss", "#IID\tMISSING_CT\tOBS_CT\tF_MISS\nS1\t0\t100\t0.0\nS2\t9\t100\t0.09\n"
                   "S3\t0\t100\t0.0\n")
    sexc = _write(tmp_path / "sex.tsv", "genotype_id\tchrX_sites\tchrX_F\tinferred_sex\nS1\t500\t0.95\tM\n"
                  "S2\t500\t0.01\tF\nS3\t10\t\tno_call\n")
    meta = _write(tmp_path / "m.tsv", "rnaseq_id\tgenotype_id\tsex\tGA\nR1\tS1\tfemale\t39\nR2\tS2\tF\t40\n"
                  "R3\tS3\tM\t38\n")
    cfg = {"cohort_id": "TC", "genotypes": {"raw_build": "GRCh37", "source": source},
           "imputation_prep": {"ref_match_min": 0.95, "sample_missing_max": 0.05}}
    cfgp = _write(tmp_path / "cfg.json", json.dumps(cfg))
    args = ["report", "--counts", str(counts), "--smiss", smiss, "--sexcheck", sexc, "--metadata", meta,
            "--config", cfgp, "--out-dir", str(tmp_path)]
    if lift_mm is not None:
        lp = _write(tmp_path / "lift.tsv", f"reason\tn\nlifted\t1100\nMismatchedRefAllele\t{lift_mm}\n")
        args += ["--liftover", lp]
    return args


def test_report_passes_and_flags_samples(tmp_path, capsys):
    IP.main(_report_inputs(tmp_path, ref_match=1000))
    rep = pd.read_csv(tmp_path / "imputation_prep_report.tsv", sep="\t")
    assert list(rep["chrom"]) == ["chr22", "total"]
    assert rep.loc[0, "ref_mismatch_dropped"] == 0
    smp = pd.read_csv(tmp_path / "imputation_prep_samples.tsv", sep="\t", dtype={"genotype_id": str})
    chk = dict(zip(smp.genotype_id, smp.sex_check))
    assert chk == {"S1": "mismatch", "S2": "ok", "S3": "no_chrX_call"}
    assert list(smp.high_missingness) == [False, True, False]
    err = capsys.readouterr().err
    assert "1 sample(s) have missingness" in err and "differs from the metadata" in err
    md = (tmp_path / "TOPMED_CHECKLIST.md").read_text()
    assert "TOPMed (Version R3 on GRCh38)" in md and "GRCh38/hg38" in md
    assert "7 days" in md and "one-time password" in md
    assert "Input build GRCh37" in md


def test_report_stops_on_wrong_build(tmp_path):
    with pytest.raises(SystemExit) as e:
        IP.main(_report_inputs(tmp_path, ref_match=500))
    assert "Wrong build?" in str(e.value)
    assert (tmp_path / "imputation_prep_report.tsv").exists()


def test_report_rate_counts_liftover_ref_mismatches(tmp_path):
    # 1000 matched SNPs, but 100 more were rejected by the liftover for a REF mismatch -> 90.9% < 95%
    with pytest.raises(SystemExit):
        IP.main(_report_inputs(tmp_path, ref_match=1000, lift_mm=100))
    IP.main(_report_inputs(tmp_path, ref_match=1000, lift_mm=10))


def test_report_wes_warning(tmp_path, capsys):
    IP.main(_report_inputs(tmp_path, ref_match=1000, source="wes"))
    assert "WES input" in capsys.readouterr().err
    assert "WES input" in (tmp_path / "TOPMED_CHECKLIST.md").read_text()

#!/usr/bin/env python3
"""xqtl_impprep.py -- --mode imputation_prep: sex check, liftover summary, report and TOPMed checklist.

  xqtl_impprep.py sexcheck  --gt chrX_gt.tsv --samples chrX_samples.txt --out sexcheck.tsv
  xqtl_impprep.py liftover  --rejected rejected.vcf.gz --lifted N --out liftover.tsv
  xqtl_impprep.py report    --counts chr*.counts.tsv --smiss sample_missing.smiss --sexcheck sexcheck.tsv
                            [--liftover liftover.tsv] [--metadata metadata.tsv] --config cfg.json --out-dir .

sexcheck: X-chromosome inbreeding coefficient per sample over non-PAR SNPs (MAF >= 0.05),
  F = 1 - observed heterozygous calls / expected (sum of 2p(1-p) over the sample's called sites),
  as in PLINK --check-sex: F > 0.8 male, F < 0.2 female, otherwise ambiguous. Haploid calls count as
  homozygous. Fewer than --min-sites sites: no call.
report: per-chromosome filter counts (imputation_prep_report.tsv), per-sample flags
  (imputation_prep_samples.tsv: missingness > sample_missing_max, sex check against metadata sex by
  genotype_id), and TOPMED_CHECKLIST.md. Stops (exit 1) when the reference-match rate is below
  ref_match_min: rate = SNPs whose alleles match GRCh38 (after REF/ALT swaps) / SNPs checked, where
  the liftover's MismatchedRefAllele rejects count as checked and unmatched.
These outputs hold sample IDs and stay at the site; they are never part of the export.
"""
import argparse
import glob
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd


# ------------------------------------------------------------------------------------------------
def sex_f(gt_lines, n_samples, min_sites=100):
    """F per sample from GT strings (one row per site, tab-separated, trailing tab allowed)."""
    obs = np.zeros(n_samples)
    exp = np.zeros(n_samples)
    sites = np.zeros(n_samples, dtype=int)
    for line in gt_lines:
        f = line.rstrip("\n").rstrip("\t").split("\t")
        if len(f) != n_samples:
            continue
        alt = np.zeros(n_samples)
        nal = np.zeros(n_samples)
        het = np.zeros(n_samples, dtype=bool)
        called = np.zeros(n_samples, dtype=bool)
        for i, g in enumerate(f):
            a = g.replace("|", "/").split("/")
            if "." in a:
                continue
            called[i] = True
            nal[i] = len(a)
            alt[i] = sum(x != "0" for x in a)
            het[i] = len(a) == 2 and a[0] != a[1]
        tot = nal[called].sum()
        if tot == 0:
            continue
        p = alt[called].sum() / tot
        if p <= 0 or p >= 1:
            continue
        e = 2 * p * (1 - p)
        exp[called] += e
        obs[het] += 1
        sites[called] += 1
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where((sites >= min_sites) & (exp > 0), 1 - obs / exp, np.nan)
    sex = np.where(np.isnan(F), "no_call", np.where(F > 0.8, "M", np.where(F < 0.2, "F", "ambiguous")))
    return F, sex, sites


def cmd_sexcheck(a):
    with open(a.samples) as fh:
        ids = [x.strip() for x in fh if x.strip()]
    with open(a.gt) as fh:
        F, sex, sites = sex_f(fh, len(ids), a.min_sites)
    pd.DataFrame({"genotype_id": ids, "chrX_sites": sites, "chrX_F": F, "inferred_sex": sex}).to_csv(
        a.out, sep="\t", index=False, float_format="%.4f")


def cmd_liftover(a):
    counts = {}
    if os.path.exists(a.rejected):
        with gzip.open(a.rejected, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                for r in line.split("\t", 7)[6].split(";"):
                    counts[r] = counts.get(r, 0) + 1
    rows = [{"reason": "lifted", "n": int(a.lifted)}] + [{"reason": k, "n": v} for k, v in sorted(counts.items())]
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)


# ------------------------------------------------------------------------------------------------
SEX = {"m": "M", "male": "M", "f": "F", "female": "F"}


def chrom_key(c):
    c = c.replace("chr", "")
    return int(c) if c.isdigit() else 99


def cmd_report(a):
    cfg = json.load(open(a.config))
    ip = cfg.get("imputation_prep", {})
    geno = cfg.get("genotypes", {})
    ref_min = float(ip.get("ref_match_min", 0.95))
    smax = float(ip.get("sample_missing_max", 0.05))
    files = [f for pat in a.counts for f in glob.glob(pat)] or a.counts
    cnt = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    wide = cnt.pivot(index="chrom", columns="step", values="n").fillna(0).astype(int)
    steps = ["input", "biallelic_snps", "pass_missing_maf", "pass_hwe", "ref_alt_swapped", "ref_match"]
    wide = wide.reindex(columns=steps, fill_value=0)
    wide = wide.loc[sorted(wide.index, key=chrom_key)]
    wide["ref_mismatch_dropped"] = wide["pass_hwe"] - wide["ref_match"]
    tot = wide.sum(numeric_only=True)
    lift = pd.read_csv(a.liftover, sep="\t") if a.liftover else None
    lift_mm = int(lift.loc[lift.reason == "MismatchedRefAllele", "n"].sum()) if lift is not None else 0
    checked = int(tot["pass_hwe"]) + lift_mm
    rate = (int(tot["ref_match"]) / checked) if checked else 0.0

    out = wide.reset_index()
    out.loc[len(out)] = ["total"] + [int(tot[c]) for c in out.columns[1:]]
    out.to_csv(os.path.join(a.out_dir, "imputation_prep_report.tsv"), sep="\t", index=False)

    # samples
    sm = pd.read_csv(a.smiss, sep=r"\s+")
    sm = sm.rename(columns={"#IID": "genotype_id", "IID": "genotype_id"})
    smp = sm[["genotype_id", "F_MISS"]].rename(columns={"F_MISS": "missing_rate"})
    smp["genotype_id"] = smp["genotype_id"].astype(str)
    smp["high_missingness"] = smp["missing_rate"] > smax
    sx = pd.read_csv(a.sexcheck, sep="\t", dtype={"genotype_id": str})
    smp = smp.merge(sx, on="genotype_id", how="left")
    smp["metadata_sex"] = ""
    if a.metadata and os.path.exists(a.metadata):
        md = pd.read_csv(a.metadata, sep="\t", dtype=str)
        if {"genotype_id", "sex"} <= set(md.columns):
            ms = md.dropna(subset=["genotype_id"]).drop_duplicates("genotype_id").set_index("genotype_id")["sex"]
            smp["metadata_sex"] = smp["genotype_id"].map(ms).fillna("").map(lambda s: SEX.get(str(s).lower(), ""))

    def chk(r):
        if r["inferred_sex"] in ("no_call", "") or pd.isna(r["inferred_sex"]):
            return "no_chrX_call"
        if not r["metadata_sex"]:
            return "no_metadata_sex"
        if r["inferred_sex"] == "ambiguous":
            return "ambiguous"
        return "ok" if r["inferred_sex"] == r["metadata_sex"] else "mismatch"
    smp["sex_check"] = smp.apply(chk, axis=1)
    smp.to_csv(os.path.join(a.out_dir, "imputation_prep_samples.tsv"), sep="\t", index=False, float_format="%.4f")

    n_hi = int(smp["high_missingness"].sum())
    n_mm = int((smp["sex_check"] == "mismatch").sum())
    chroms = [c for c in wide.index if wide.loc[c, "ref_match"] > 0]
    warnings = []
    if geno.get("source") == "wes":
        warnings.append(f"WES input: {int(tot['ref_match']):,} SNPs pass. Exome variants are sparse and clustered; "
                        "imputation accuracy outside exons will be low. Merge with array data if available.")
    if n_hi:
        warnings.append(f"{n_hi} sample(s) have missingness > {smax} (imputation_prep_samples.tsv). "
                        "Consider removing them before upload.")
    if n_mm:
        warnings.append(f"{n_mm} sample(s) have a chrX-inferred sex that differs from the metadata: "
                        "check for sample swaps before upload.")
    if 0 < int(tot["ref_match"]) < 100_000 and geno.get("source") != "wes":
        warnings.append(f"only {int(tot['ref_match']):,} SNPs pass in total; TOPMed recommends genome-wide array content.")
    write_checklist(a.out_dir, cfg, chroms, tot, rate, lift, warnings, len(smp))
    for w in warnings:
        print("WARNING: " + w, file=sys.stderr)
    print(f"reference-match rate {rate:.4f} ({int(tot['ref_match'])} / {checked}); {len(chroms)} chromosome VCFs")
    if rate < ref_min:
        sys.exit(f"ERROR: only {100 * rate:.1f}% of SNPs match the GRCh38 reference (minimum {100 * ref_min:.0f}%). "
                 f"Wrong build? genotypes.raw_build is '{geno.get('raw_build', 'GRCh38')}'. Data on the Illumina "
                 "TOP/BOT strand must be converted to the forward strand first (bcftools +fixref -m top). "
                 "No VCFs were published; see imputation_prep_report.tsv in the task work directory.")


def write_checklist(out_dir, cfg, chroms, tot, rate, lift, warnings, n_samples):
    g = cfg.get("genotypes", {})
    outdir = cfg.get("outdir", "results")
    lines = [
        "# TOPMed Imputation Server checklist",
        "",
        f"Cohort `{cfg.get('cohort_id', 'NA')}`: {n_samples} samples, {len(chroms)} chromosome VCFs, "
        f"{int(tot['ref_match']):,} SNPs; reference-match rate {100 * rate:.2f}%. "
        f"Input build {g.get('raw_build', 'GRCh38')}{' (lifted to GRCh38 with Picard LiftoverVcf)' if lift is not None else ''}.",
        "",
    ]
    if warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in warnings] + [""]
    lines += [
        "## 1. Upload",
        "",
        "Upload is manual: genotypes leave the site only to the TOPMed Imputation Server, under your institution's approval.",
        "",
        "1. Log in at https://imputation.biodatacatalyst.nhlbi.nih.gov and start a new imputation job.",
        f"2. Upload the files `chr*.vcf.gz` from `{outdir}/imputation_prep/` (one file per chromosome; the `.tbi` files are not needed).",
        "3. Settings:",
        "   - Reference panel: **TOPMed (Version R3 on GRCh38)**.",
        "   - Array build: **GRCh38/hg38** (the files are already on GRCh38, whatever the input build was).",
        "   - rsq filter: **off** (the pipeline applies `qc.rsq_min` itself, as the atlas did).",
        "   - Phasing: **Eagle v2.4** (the files are unphased; choose *No phasing* only if you phased them yourself).",
        "   - Population: **all** / vs. TOPMed panel.",
        "   - Mode: **Quality Control & Imputation**.",
        "   - AES 256 encryption: optional. AES-256 archives need 7-Zip to open.",
        "4. Submit, then read the server's QC report. Many strand flips or allele mismatches mean the build or strand is wrong.",
        "",
        "## 2. Download",
        "",
        "- The server emails a one-time password for the result archives. It is not stored on the server and cannot be re-sent: keep it in your password manager, never in `cohort.yml`, scripts or e-mail threads.",
        "- Results are deleted after 7 days. Download every `chr_*.zip` within that time (the results page gives `wget` links).",
        "- Unzip with the password into one directory on the cluster, for example `genotypes/imputed/`. Keep the `*.dose.vcf.gz` (and `*.info.gz`) files.",
        "",
        "## 3. Run the site analysis",
        "",
        "Set in `cohort.yml`:",
        "",
        "```yaml",
        "genotypes:",
        "  imputed_dir: /path/to/genotypes/imputed",
        "  imputed_pattern: chr{chr}.dose.vcf.gz   # TOPMed file names",
        f"  source: {g.get('source', 'array')}",
        "```",
        "",
        "Then run `--mode site` (see the walkthrough). The pre-imputation files and this checklist stay at the site.",
        "",
    ]
    with open(os.path.join(out_dir, "TOPMED_CHECKLIST.md"), "w") as fh:
        fh.write("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sexcheck")
    s.add_argument("--gt", required=True)
    s.add_argument("--samples", required=True)
    s.add_argument("--min-sites", type=int, default=100)
    s.add_argument("--out", required=True)
    l = sub.add_parser("liftover")
    l.add_argument("--rejected", required=True)
    l.add_argument("--lifted", required=True)
    l.add_argument("--out", required=True)
    r = sub.add_parser("report")
    r.add_argument("--counts", nargs="+", required=True)
    r.add_argument("--smiss", required=True)
    r.add_argument("--sexcheck", required=True)
    r.add_argument("--liftover")
    r.add_argument("--metadata")
    r.add_argument("--config", required=True)
    r.add_argument("--out-dir", default=".")
    a = ap.parse_args(argv)
    {"sexcheck": cmd_sexcheck, "liftover": cmd_liftover, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()

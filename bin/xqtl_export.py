#!/usr/bin/env python3
"""xqtl_export.py -- assemble the shareable summary files (A-F, I per stratum; J per cohort) and scan for leaks.

  xqtl_export.py stratum  --anc EUR --cohort BiSC --in STAGING --qc-metrics qc.tsv [--cell-props ct.tsv]
                          [--format parquet|tsv.gz] [--round-digits N] [--min-mac 0] --out EUR
  xqtl_export.py finalize --cohort BiSC --strata cohort_strata.tsv --ids metadata.tsv samplesheet.csv ...
                          --image-sha256 SHA --restricted-version V --pipeline-version V --out export STRATUM_DIRS...
  xqtl_export.py scan     --dir DIR --ids FILES...  (leak scan only; exit 1 on any finding)

STAGING holds the stratum's pipeline outputs (names as written by the site workflow):
  {ANC}.variant_qc.tsv  {ANC}.geno_qc_counts.tsv  {ANC}.samples.tsv  {ANC}_ld/
  {ANC}_{mod}.features.tsv  {ANC}_{mod}.norm_summary.tsv  {ANC}_{mod}.nominal/
  {ANC}_{arm}.cis.parquet  {ANC}_{arm}.independent.parquet  {ANC}_{arm}.settings.json
  {ANC}_{arm}.covariates.tsv (only the covariate NAMES are exported)
Files A-F and I are summary statistics; nothing per sample is written. The finalize step copies
the stratum folders into OUT/, writes I_strata.tsv, MANIFEST.tsv (J) and DATA_DICTIONARY.md, then
runs the leak scan, which cannot be disabled: any sample identifier (RNA-seq, genotype or VCF
sample ID) in a file name, column name or string value, any npz other than one square matrix
'R', or any unexpected file type fails the export.
"""
import argparse
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd

VARIANT_RE = re.compile(r"^chr[0-9XYM]{1,2}:\d+:[ACGTN]+:[ACGTN]+$")
ALLOWED_EXT = (".parquet", ".tsv.gz", ".tsv", ".npz", ".json", ".md")

DESCRIPTIONS = {
    "A": "nominal cis associations, all variant-phenotype pairs (per modality and chromosome)",
    "B": "permutation (lead per phenotype/group, q-values) and conditionally independent signals",
    "C": "in-sample LD (signed r, float16) for MAF >= 0.05, R2 >= 0.8 variants, 3 Mb windows / 1 Mb step",
    "D": "per-variant QC summary in the stratum",
    "E": "phenotype features: coordinates, group, detection fraction, filter status",
    "F": "model settings: normalization, HCP grid and selected k, covariate names (no values)",
    "I": "aggregate counts, sequencing-metric quantiles, cell-type proportion summaries",
    "J": "manifest and data dictionary",
}
COLUMNS = {
    "A": [("phenotype_id", "feature ID (modality-specific; GENE__feature for grouped modalities)"),
          ("variant_id", "GRCh38 chr:pos:ref:alt"), ("chrom", "chromosome"), ("pos", "GRCh38 position"),
          ("ref", "reference allele"), ("alt", "alternate allele"), ("effect_allele", "allele whose dosage is tested (= alt)"),
          ("other_allele", "= ref"), ("tss_distance", "variant position minus feature TSS (bp)"),
          ("eaf", "effect-allele frequency among mapped samples"), ("maf", "minor-allele frequency"),
          ("ma_count", "minor-allele count"), ("n", "samples in the test"), ("beta", "effect per effect allele (normalized phenotype units)"),
          ("se", "standard error of beta"), ("pval", "nominal p-value")],
    "D": [("variant_id", "GRCh38 chr:pos:ref:alt"), ("chrom", ""), ("pos", ""), ("ref", ""), ("alt", ""),
          ("eaf", "ALT dosage frequency"), ("maf", ""), ("n", "individuals"), ("hwe_p", "Hardy-Weinberg p (hard calls)"),
          ("rsq", "imputation R2 (TOPMed INFO/R2)"), ("typed", "1 if genotyped on the array")],
}


# --------------------------------------------------------------------------------------------
def write_table(df, path_noext, fmt, round_digits=None):
    if round_digits is not None:
        num = df.select_dtypes(include=[np.floating]).columns
        df[num] = df[num].round(round_digits)
    if fmt == "parquet":
        p = path_noext + ".parquet"
        df.to_parquet(p, index=False)
    else:
        p = path_noext + ".tsv.gz"
        df.to_csv(p, sep="\t", index=False, float_format="%.6g", compression="gzip")
    return p


def split_variant(v):
    s = v.str.split(":", expand=True)
    return s[0], s[1].astype(np.int64), s[2], s[3]


def stratum(a):
    S, anc = a.input, a.anc
    out = a.out
    os.makedirs(out, exist_ok=True)
    pre = os.path.join(out, f"{a.cohort}_{anc}")

    # D: variant QC
    vq = pd.read_csv(os.path.join(S, f"{anc}.variant_qc.tsv"), sep="\t")
    vq = vq.rename(columns={"af": "eaf"})
    write_table(vq[["variant_id", "chrom", "pos", "ref", "alt", "eaf", "maf", "n", "hwe_p", "rsq", "typed"]],
                pre + "_D_variant_qc", a.format, a.round_digits)

    # C: LD windows
    ld = os.path.join(S, f"{anc}_ld")
    if os.path.isdir(ld):
        cdir = os.path.join(out, "C_ld")
        os.makedirs(cdir, exist_ok=True)
        for f in sorted(os.listdir(ld)):
            dst = f.replace(f"{anc}_ld", f"{a.cohort}_{anc}_C_ld", 1)
            shutil.copyfile(os.path.join(ld, f), os.path.join(cdir, dst))

    # E + F + A + B per modality / arm
    mods = sorted({os.path.basename(p)[len(anc) + 1:-len(".features.tsv")]
                   for p in glob.glob(os.path.join(S, f"{anc}_*.features.tsv"))})
    feats, norm = [], []
    for m in mods:
        f = pd.read_csv(os.path.join(S, f"{anc}_{m}.features.tsv"), sep="\t")
        f.insert(0, "modality", m)
        feats.append(f)
        ns = os.path.join(S, f"{anc}_{m}.norm_summary.tsv")
        if os.path.exists(ns):
            norm.append(pd.read_csv(ns, sep="\t"))
    if feats:
        write_table(pd.concat(feats, ignore_index=True), pre + "_E_features", a.format)

    n_by_mod = {}
    settings = []
    for js in sorted(glob.glob(os.path.join(S, f"{anc}_*.settings.json"))):
        arm = os.path.basename(js)[len(anc) + 1:-len(".settings.json")]
        with open(js) as fh:
            st = json.load(fh)
        cov = os.path.join(S, f"{anc}_{arm}.covariates.tsv")
        if os.path.exists(cov):
            names = pd.read_csv(cov, sep="\t", usecols=[0]).iloc[:, 0].astype(str).tolist()
            st["covariates"] = names
            st["n_covariates"] = len(names)
        st["arm"] = arm
        settings.append(st)
        n_by_mod[st.get("modality", arm)] = st.get("n_samples")
    fobj = {"cohort": a.cohort, "ancestry": anc, "arms": settings,
            "normalization": norm and pd.concat(norm, ignore_index=True).to_dict(orient="records") or []}
    with open(pre + "_F_model_settings.json", "w") as fh:
        json.dump(fobj, fh, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))

    bdir = os.path.join(out, "B_mapping")
    os.makedirs(bdir, exist_ok=True)
    for kind in ("cis", "independent"):
        for p in sorted(glob.glob(os.path.join(S, f"{anc}_*.{kind}.parquet"))):
            arm = os.path.basename(p)[len(anc) + 1:-len(f".{kind}.parquet")]
            df = pd.read_parquet(p)
            if len(df) == 0:
                df = pd.DataFrame(columns=["phenotype_id", "variant_id"])
            df.insert(0, "arm", arm)
            label = "permutation" if kind == "cis" else "independent"
            write_table(df, os.path.join(bdir, f"{a.cohort}_{anc}_B_{label}_{arm}"), a.format, a.round_digits)

    adir = os.path.join(out, "A_nominal")
    os.makedirs(adir, exist_ok=True)
    for d in sorted(glob.glob(os.path.join(S, f"{anc}_*.nominal"))):
        m = os.path.basename(d)[len(anc) + 1:-len(".nominal")]
        n = n_by_mod.get(m)
        for p in sorted(glob.glob(os.path.join(d, "*.cis_qtl_pairs.*.parquet"))):
            chrom = p.rsplit(".", 2)[-2]
            x = pd.read_parquet(p)
            if a.min_mac > 0:
                x = x[x["ma_count"] >= a.min_mac]
            c, pos, ref, alt = split_variant(x["variant_id"])
            o = pd.DataFrame({"phenotype_id": x["phenotype_id"].values, "variant_id": x["variant_id"].values,
                              "chrom": c.values, "pos": pos.values, "ref": ref.values, "alt": alt.values,
                              "effect_allele": alt.values, "other_allele": ref.values,
                              "tss_distance": x["start_distance"].values, "eaf": x["af"].values,
                              "maf": np.minimum(x["af"].values, 1 - x["af"].values),
                              "ma_count": x["ma_count"].values, "n": n if n is not None else np.nan,
                              "beta": x["slope"].values, "se": x["slope_se"].values, "pval": x["pval_nominal"].values})
            write_table(o, os.path.join(adir, f"{a.cohort}_{anc}_A_nominal_{m}_{chrom}"), a.format, a.round_digits)

    # I: aggregate counts and distributions
    rows = []
    gq = os.path.join(S, f"{anc}.geno_qc_counts.tsv")
    if os.path.exists(gq):
        for r in pd.read_csv(gq, sep="\t").itertuples():
            rows.append({"table": "genotype_qc", "item": r.step, "statistic": "n_variants", "value": r.n_variants})
            rows.append({"table": "genotype_qc", "item": r.step, "statistic": "n_individuals", "value": r.n_samples})
    for r in (pd.concat(norm, ignore_index=True).to_dict(orient="records") if norm else []):
        for k in ("features_in", "features_kept", "samples_in", "connectivity_outliers", "replicates_averaged", "individuals"):
            rows.append({"table": "normalization", "item": r["modality"], "statistic": k, "value": r[k]})
    smp = pd.read_csv(os.path.join(S, f"{anc}.samples.tsv"), sep="\t", dtype=str)
    qs = [0.05, 0.25, 0.5, 0.75, 0.95]

    def quant(df, table):
        for col in df.columns:
            v = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(v) < 10:
                continue
            for q in qs:
                rows.append({"table": table, "item": col, "statistic": f"q{int(q * 100):02d}", "value": float(v.quantile(q))})
            rows.append({"table": table, "item": col, "statistic": "mean", "value": float(v.mean())})
            rows.append({"table": table, "item": col, "statistic": "sd", "value": float(v.std())})
    qm = pd.read_csv(a.qc_metrics, sep="\t", dtype={"sample": str})
    quant(qm[qm["sample"].isin(smp.rnaseq_id)].drop(columns="sample"), "sequencing_metrics")
    if a.cell_props:
        ct = pd.read_csv(a.cell_props, sep="\t", dtype={"sample": str})
        quant(ct[ct["sample"].isin(smp.rnaseq_id)].drop(columns="sample"), "cell_type_proportions")
    write_table(pd.DataFrame(rows, columns=["table", "item", "statistic", "value"]), pre + "_I_summary", a.format)
    print(f"{out}: stratum export written")


# --------------------------------------------------------------------------------------------
def load_ids(paths):
    ids = set()
    for p in paths:
        if not os.path.exists(p):
            continue
        sep = "," if p.endswith(".csv") else "\t"
        if p.endswith((".txt", ".psam")) or os.path.basename(p).endswith("vcf_samples.txt"):
            with open(p) as fh:
                for line in fh:
                    t = line.strip().split("\t")[0].lstrip("#")
                    if t and t not in ("IID", "FID"):
                        ids.add(t)
            continue
        df = pd.read_csv(p, sep=sep, dtype=str)
        for c in ("rnaseq_id", "genotype_id", "sample_id", "sample"):
            if c in df.columns:
                ids.update(df[c].dropna().astype(str).str.strip())
    ids.discard("")
    return ids


class Scanner:
    def __init__(self, ids):
        self.exact = {i for i in ids if i.isdigit() or len(i) < 3}
        tok = sorted((i for i in ids if i not in self.exact), key=len, reverse=True)
        self.rx = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(map(re.escape, tok)) + r")(?![A-Za-z0-9])") if tok else None
        self.seen = set()
        self.findings = []

    def hit(self, s):
        if s in self.exact:
            return s
        if self.rx is not None:
            m = self.rx.search(s)
            if m:
                return m.group(1)
        return None

    def strings(self, values, where):
        for v in values:
            if v in self.seen:
                continue
            self.seen.add(v)
            if VARIANT_RE.match(v):
                continue
            h = self.hit(v)
            if h:
                self.findings.append((where, h))
                return

    def table(self, df, where):
        self.strings([str(c) for c in df.columns], where + " [column names]")
        # Scan every column that is not numeric/boolean/datetime. A whitelist of text dtypes is
        # unsafe: pandas >= 3 infers a new 'str' dtype that an object/string check misses.
        for c in df.columns:
            s = df[c]
            if (pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s)
                    or pd.api.types.is_datetime64_any_dtype(s)):
                continue
            self.strings(pd.unique(s.dropna().astype(str)), f"{where} [column {c}]")

    def file(self, path, rel):
        self.strings([os.path.basename(rel)], rel + " [file name]")
        if path.endswith(".npz"):
            try:
                z = np.load(path, allow_pickle=False)
                keys = list(z.keys())
                R = z["R"] if keys == ["R"] else None
                if R is None or R.ndim != 2 or R.shape[0] != R.shape[1]:
                    self.findings.append((rel, f"npz must hold exactly one square matrix 'R' (keys {keys})"))
            except Exception as e:  # unreadable or pickled content
                self.findings.append((rel, f"unreadable npz: {e}"))
        elif path.endswith(".parquet"):
            import pyarrow.parquet as pq
            pf = pq.ParquetFile(path)
            for i in range(pf.num_row_groups):
                self.table(pf.read_row_group(i).to_pandas(), rel)
        elif path.endswith((".tsv.gz", ".tsv")):
            for chunk in pd.read_csv(path, sep="\t", dtype=str, chunksize=500_000):
                self.table(chunk, rel)
        elif path.endswith((".json", ".md")):
            with open(path) as fh:
                txt = fh.read()
            for tokn in set(re.findall(r"[A-Za-z0-9_.\-]+", txt)):
                h = self.hit(tokn)
                if h:
                    self.findings.append((rel, h))
                    break
        else:
            self.findings.append((rel, "unexpected file type"))


def scan_dir(root, ids):
    sc = Scanner(ids)
    for dp, _, fs in os.walk(root):
        for f in sorted(fs):
            p = os.path.join(dp, f)
            sc.file(p, os.path.relpath(p, root))
    return sc.findings


def count_rows(p):
    if p.endswith(".parquet"):
        import pyarrow.parquet as pq
        return pq.ParquetFile(p).metadata.num_rows
    if p.endswith(".tsv.gz"):
        with gzip.open(p, "rt") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    if p.endswith(".tsv"):
        with open(p) as fh:
            return max(0, sum(1 for _ in fh) - 1)
    if p.endswith(".npz"):
        return int(np.load(p)["R"].shape[0])
    return ""


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def file_class(name):
    m = re.search(r"_([A-J])_", name)
    return m.group(1) if m else ("J" if name in ("MANIFEST.tsv", "DATA_DICTIONARY.md") else "?")


def finalize(a):
    root = os.path.join(a.out, a.cohort)
    if os.path.exists(a.out):
        shutil.rmtree(a.out)
    os.makedirs(root)
    for d in a.strata_dirs:
        shutil.copytree(d, os.path.join(root, os.path.basename(d.rstrip("/"))), symlinks=False)
    st = pd.read_csv(a.strata, sep="\t")
    st.insert(0, "cohort", a.cohort)
    st.to_csv(os.path.join(root, f"{a.cohort}_cohort_I_strata.tsv"), sep="\t", index=False)

    lines = [f"# Data dictionary: {a.cohort} xQTL summary export", "",
             f"Pipeline version {a.pipeline_version}; image sha256 {a.image_sha256}; restricted file version {a.restricted_version}.",
             "All coordinates GRCh38; variant IDs chr:pos:ref:alt; effects are per ALT allele.", "",
             "| Class | Content |", "|---|---|"]
    lines += [f"| {k} | {v} |" for k, v in DESCRIPTIONS.items()]
    for k, cols in COLUMNS.items():
        lines += ["", f"## {k} columns", "", "| column | meaning |", "|---|---|"] + [f"| {c} | {d} |" for c, d in cols]
    lines += ["", "B tables are tensorQTL map_cis / map_independent output (plus arm, qval, qvalue_method).",
              "C: one .npz (key R) per window with a matching .variants.tsv.gz giving the row order."]
    with open(os.path.join(root, "DATA_DICTIONARY.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    findings = scan_dir(a.out, load_ids(a.ids))
    if findings:
        pd.DataFrame(findings, columns=["where", "finding"]).to_csv(a.leak_report, sep="\t", index=False)
        shutil.rmtree(a.out)
        sys.exit(f"xqtl_export.py: LEAK SCAN FAILED with {len(findings)} finding(s); the export was deleted. "
                 f"See {a.leak_report} (local only, do not share).")

    recs = []
    for dp, _, fs in os.walk(root):
        for f in sorted(fs):
            p = os.path.join(dp, f)
            rel = os.path.relpath(p, a.out)
            cls = file_class(f)
            recs.append({"path": rel, "sha256": sha256(p), "bytes": os.path.getsize(p), "rows": count_rows(p),
                         "data_class": cls, "description": DESCRIPTIONS.get(cls, ""),
                         "image_sha256": a.image_sha256, "restricted_version": a.restricted_version,
                         "pipeline_version": a.pipeline_version})
    pd.DataFrame(recs).sort_values("path").to_csv(os.path.join(root, "MANIFEST.tsv"), sep="\t", index=False)
    pd.DataFrame(columns=["where", "finding"]).to_csv(a.leak_report, sep="\t", index=False)
    print(f"{a.out}: {len(recs)} files, leak scan passed")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("stratum")
    s.add_argument("--anc", required=True)
    s.add_argument("--cohort", required=True)
    s.add_argument("--in", dest="input", required=True)
    s.add_argument("--qc-metrics", required=True)
    s.add_argument("--cell-props")
    s.add_argument("--format", choices=["parquet", "tsv.gz"], default="parquet")
    s.add_argument("--round-digits", type=int)
    s.add_argument("--min-mac", type=int, default=0)
    s.add_argument("--out", required=True)
    f = sub.add_parser("finalize")
    f.add_argument("--cohort", required=True)
    f.add_argument("--strata", required=True)
    f.add_argument("--ids", nargs="+", required=True)
    f.add_argument("--image-sha256", default="unknown")
    f.add_argument("--restricted-version", default="unknown")
    f.add_argument("--pipeline-version", default="unknown")
    f.add_argument("--leak-report", default="leak_report.tsv")
    f.add_argument("--out", required=True)
    f.add_argument("strata_dirs", nargs="*")
    c = sub.add_parser("scan")
    c.add_argument("--ids", nargs="+", required=True)
    c.add_argument("--dir", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "stratum":
        stratum(a)
    elif a.cmd == "finalize":
        finalize(a)
    else:
        fnd = scan_dir(a.dir, load_ids(a.ids))
        for w, h in fnd:
            print(f"{w}\t{h}")
        sys.exit(1 if fnd else 0)


if __name__ == "__main__":
    main()

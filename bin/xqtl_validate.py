#!/usr/bin/env python3
"""Validate a placenta-xqtl-federated run before any work starts.

Checks, in order, and reports every problem found (not just the first):
  1. cohort parameters against assets/schema/cohort.schema.json (defaults filled in);
  2. mode-specific required inputs;
  3. samplesheet (sample_id,fastq_1,fastq_2): columns, unique IDs, files exist, layout;
  4. metadata TSV: required columns, sex coding, GA range, optional exposures;
  5. links between samplesheet and metadata;
  6. enabled GxE exposures have a column and enough non-missing values; enabled sensitivity
     covariates have a numeric column.

Writes the resolved parameters (defaults filled) as JSON and a TSV report.
Exit status 0 = valid, 1 = invalid (messages on stderr), 2 = usage error.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import sys

import jsonschema

REQUIRED_META = ["rnaseq_id", "genotype_id", "sex", "GA"]
OPTIONAL_META = ["ppBMI", "gdm", "ogtt", "maternal_age", "parity", "smoking"]
SEX_CODES = {"m": "M", "male": "M", "f": "F", "female": "F"}
RESTRICTED_HINT = (
    "The restricted file (xqtl-restricted-<version>.tar.gz) is supplied by MDACC by email "
    "after the data-sharing agreement is approved. If your site has not received it, request it "
    "from the scientific contact (Sean T. Bresnahan, stbresnahan@mdanderson.org), stating your "
    "cohort, the pipeline version and the approved agreement. See docs/SETUP.md."
)

MODE_REQUIRES = {
    "prepare_references": ["site_refs"],
    "imputation_prep": ["genotypes.raw"],
    # site_refs is not needed when phenotypes_dir supplies ready phenotype inputs (see main())
    "site": ["samplesheet", "metadata", "restricted", "site_refs", "genotypes.imputed_dir"],
    "phenotyping_only": ["samplesheet", "metadata", "restricted", "site_refs"],
    "export": ["samplesheet", "metadata"],
    "central": ["central.studies"],
    "smoke": [],
}


def fill_defaults(schema: dict, inst):
    """Recursively insert schema defaults for missing object properties.

    Map-type objects with a default (e.g. gxe.exposures) are merged key by key: the
    user's entries override or extend the defaults instead of replacing the whole map,
    so enabling ``ppBMI`` does not silently drop the default ``GA`` exposure.
    """
    if not isinstance(schema, dict):
        return inst
    if schema.get("type") == "object" or "properties" in schema or "additionalProperties" in schema:
        if inst is None:
            inst = copy.deepcopy(schema.get("default", {}))
        if not isinstance(inst, dict):
            return inst
        props = schema.get("properties", {})
        for key, sub in props.items():
            if not isinstance(sub, dict):
                continue
            is_obj = sub.get("type") == "object" or "properties" in sub or "additionalProperties" in sub
            if key not in inst:
                if "default" in sub:
                    inst[key] = copy.deepcopy(sub["default"])
                elif is_obj:
                    inst[key] = {}
            elif is_obj and isinstance(sub.get("default"), dict) and isinstance(inst[key], dict) \
                    and isinstance(sub.get("additionalProperties"), dict):
                merged = copy.deepcopy(sub["default"])
                for k, v in inst[key].items():
                    merged[k] = {**merged.get(k, {}), **v} if isinstance(v, dict) else v
                inst[key] = merged
            if key in inst and is_obj:
                inst[key] = fill_defaults(sub, inst[key])
        addl = schema.get("additionalProperties")
        if isinstance(addl, dict):
            for key in list(inst.keys()):
                if key not in props:
                    inst[key] = fill_defaults(addl, inst[key])
    return inst


def get_path(d: dict, dotted: str):
    cur = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or cur.get(part) in (None, ""):
            return None
        cur = cur[part]
    return cur


class Report:
    def __init__(self):
        self.rows: list[tuple[str, str, str]] = []

    def ok(self, check, detail=""):
        self.rows.append((check, "PASS", detail))

    def err(self, check, detail):
        self.rows.append((check, "FAIL", detail))

    def warn(self, check, detail):
        self.rows.append((check, "WARN", detail))

    @property
    def failed(self):
        return any(s == "FAIL" for _, s, _ in self.rows)

    def write(self, path):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["check", "status", "detail"])
            w.writerows(self.rows)


def read_table(path: str, sep: str) -> tuple[list[str], list[dict]]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh, delimiter=sep)
        header = [h.strip() for h in (rdr.fieldnames or [])]
        rows = [{(k or "").strip(): (v or "").strip() for k, v in r.items()} for r in rdr]
    return header, rows


def resolve(path: str, base: str) -> str:
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(base, path))


def check_samplesheet(path, params, rep, check_files=True):
    base = os.path.dirname(os.path.abspath(path))
    try:
        header, rows = read_table(path, ",")
    except OSError as e:
        rep.err("samplesheet.read", f"cannot read {path}: {e}")
        return []
    missing = [c for c in ("sample_id", "fastq_1", "fastq_2") if c not in header]
    if missing:
        rep.err("samplesheet.columns", f"missing column(s) {missing}; expected sample_id,fastq_1,fastq_2")
        return []
    ids = [r["sample_id"] for r in rows]
    dups = sorted({i for i in ids if ids.count(i) > 1})
    if dups:
        rep.err("samplesheet.unique", f"duplicate sample_id: {dups[:10]}")
    if any(not i for i in ids):
        rep.err("samplesheet.ids", "empty sample_id value(s)")
    layout = (params.get("rnaseq") or {}).get("read_layout", "paired")
    bad_layout, missing_files = [], []
    for r in rows:
        if layout == "paired" and not r["fastq_2"]:
            bad_layout.append(r["sample_id"])
        if layout == "single" and r["fastq_2"]:
            bad_layout.append(r["sample_id"])
        if check_files:
            for col in ("fastq_1", "fastq_2"):
                if r[col] and not os.path.exists(resolve(r[col], base)):
                    missing_files.append(r[col])
    if bad_layout:
        rep.err("samplesheet.layout", f"rnaseq.read_layout={layout} but {len(bad_layout)} row(s) disagree, e.g. {bad_layout[:5]}")
    if missing_files:
        rep.err("samplesheet.files", f"{len(missing_files)} FASTQ path(s) not found, e.g. {missing_files[:3]}")
    if not rep.failed:
        rep.ok("samplesheet", f"{len(rows)} samples, layout={layout}")
    return ids


def as_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def check_metadata(path, params, rep, sample_ids):
    try:
        header, rows = read_table(path, "\t")
    except OSError as e:
        rep.err("metadata.read", f"cannot read {path}: {e}")
        return
    missing = [c for c in REQUIRED_META if c not in header]
    if missing:
        rep.err("metadata.columns", f"missing required column(s) {missing}; required: {REQUIRED_META}")
        return
    rn = [r["rnaseq_id"] for r in rows]
    dups = sorted({i for i in rn if rn.count(i) > 1})
    if dups:
        rep.err("metadata.unique", f"duplicate rnaseq_id: {dups[:10]}")
    bad_sex = [r["rnaseq_id"] for r in rows if r["sex"] and r["sex"].lower() not in SEX_CODES]
    if bad_sex:
        rep.err("metadata.sex", f"sex must be M/F (or male/female); bad rows e.g. {bad_sex[:5]}")
    n_missing_sex = sum(1 for r in rows if not r["sex"])
    if n_missing_sex:
        rep.warn("metadata.sex_missing", f"{n_missing_sex} sample(s) with missing sex are imputed with the mode")
    bad_ga = [r["rnaseq_id"] for r in rows if as_float(r["GA"]) is None or not (20 <= as_float(r["GA"]) <= 45)]
    if bad_ga:
        rep.err("metadata.GA", f"GA (weeks at delivery) must be numeric in [20, 45]; bad rows e.g. {bad_ga[:5]}")
    extra = [c for c in header if c not in REQUIRED_META + OPTIONAL_META]
    if extra:
        rep.warn("metadata.extra_columns", f"ignored column(s): {extra[:10]}")
    if sample_ids:
        not_in_meta = sorted(set(sample_ids) - set(rn))
        if not_in_meta:
            rep.err("link.samplesheet_metadata", f"{len(not_in_meta)} samplesheet sample_id(s) missing from metadata rnaseq_id, e.g. {not_in_meta[:5]}")
        else:
            rep.ok("link.samplesheet_metadata", "every sample has metadata")
    # exposures
    gxe = params.get("gxe") or {}
    if gxe.get("enabled", True):
        for name, spec in (gxe.get("exposures") or {}).items():
            if not spec.get("enabled"):
                continue
            col = spec.get("column") or name
            if col not in header:
                rep.err(f"gxe.{name}", f"exposure '{name}' is enabled but column '{col}' is not in metadata")
                continue
            vals = [as_float(r[col]) for r in rows if r[col] not in ("", "NA", "na", "NaN")]
            if any(v is None for v in vals):
                rep.err(f"gxe.{name}", f"column '{col}' has non-numeric values")
                continue
            if spec.get("type") == "binary" and any(v not in (0.0, 1.0) for v in vals):
                rep.err(f"gxe.{name}", f"binary exposure '{col}' must be coded 0/1")
                continue
            n = len(vals)
            if n < gxe.get("min_n", 50):
                rep.warn(f"gxe.{name}", f"only {n} non-missing values; the exposure will be skipped in groups below gxe.min_n")
            else:
                rep.ok(f"gxe.{name}", f"{n} non-missing values")
        for name, spec in (gxe.get("sensitivity") or {}).items():
            if not spec.get("enabled") or name == "ancestry_pcs":
                continue
            col = spec.get("column") or name
            if col not in header:
                rep.err(f"gxe.sensitivity.{name}", f"sensitivity covariate '{name}' is enabled but column '{col}' is not in metadata")
            elif any(as_float(r[col]) is None for r in rows if r[col] not in ("", "NA", "na", "NaN")):
                rep.err(f"gxe.sensitivity.{name}", f"column '{col}' has non-numeric values")
            else:
                rep.ok(f"gxe.sensitivity.{name}", f"column '{col}' present")
    if not any(s == "FAIL" and c.startswith("metadata") for c, s, _ in rep.rows):
        rep.ok("metadata", f"{len(rows)} rows")


def _restricted_ok(path, params, rep):
    """Verify the restricted archive (version + checksums) with bin/xqtl_restricted.py."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import xqtl_restricted  # noqa: E402  (sibling script in bin/)

    ok, msgs, _ = xqtl_restricted.verify(path, params.get("expected_restricted_version"))
    if not ok:
        rep.err("required.restricted", f"{path}: " + "; ".join(msgs) + ". " + RESTRICTED_HINT)
    return ok


ISOFORM_MODALITIES = ("isoform_expression", "isoform_usage")


def check_phenotypes_dir(params, rep):
    d = params["phenotypes_dir"]
    if not os.path.isdir(d):
        rep.err("phenotypes_dir", f"directory not found: {d}")
        return
    mods = (params.get("mapping") or {}).get("modalities") or []
    need = [f"{m}.bed.gz" for m in mods] + ["qc_metrics.tsv"]
    missing = [f for f in need if not os.path.exists(os.path.join(d, f))]
    if missing:
        rep.err("phenotypes_dir", f"missing file(s) in {d}: {missing}")
    else:
        rep.ok("phenotypes_dir", f"{len(mods)} modality BED(s) in {d}")


SITE_REFS_FILES = ["MANIFEST.tsv", "genome/genome.fa", "genome/genome.fa.fai", "genome/genome.dict",
                   "annotation/annotation.gtf", "annotation/exonic.gtf", "annotation/intronic.gtf",
                   "annotation/refFlat.txt", "annotation/gene_gc_length.tsv", "star/SA",
                   "salmon/expression", "irfinder/IRFinder/ref-cover.bed"]


def check_site_refs(params, rep):
    """Indices written by --mode prepare_references (needed by phenotyping)."""
    d = params.get("site_refs")
    if not d:
        return
    if not os.path.isdir(d):
        rep.err("site_refs", f"directory not found: {d}. Build it with --mode prepare_references")
        return
    missing = [f for f in SITE_REFS_FILES if not os.path.exists(os.path.join(d, f))]
    if missing:
        rep.err("site_refs", f"incomplete site_refs {d} (missing {missing[:4]}); rerun --mode prepare_references")
        return
    with open(os.path.join(d, "MANIFEST.tsv")) as fh:
        man = dict(line.rstrip("\n").split("\t", 1) for line in fh if "\t" in line)
    rl = (params.get("rnaseq") or {}).get("read_length")
    if rl is not None and str(man.get("read_length")) != str(rl):
        rep.warn("site_refs.read_length", f"STAR index built for read length {man.get('read_length')}, "
                 f"rnaseq.read_length is {rl}; alignment still works but rebuilding is recommended")
    if not os.path.exists(os.path.join(d, "editing", "edit_sites.bed")):
        mods = (params.get("mapping") or {}).get("modalities") or []
        if "rna_editing" in mods:
            rep.warn("site_refs.editing", "no REDIportal editing sites in site_refs: rna_editing will be skipped")
    rep.ok("site_refs", d)


def _raw_files(raw: str):
    """Files making up genotypes.raw (VCF/BCF file or PLINK 1/2 prefix); None if not recognised."""
    if raw.endswith((".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")):
        return [raw]
    if os.path.exists(raw + ".bed") or os.path.exists(raw + ".bim"):
        return [raw + e for e in (".bed", ".bim", ".fam")]
    if os.path.exists(raw + ".pgen"):
        pvar = raw + ".pvar" if os.path.exists(raw + ".pvar") else raw + ".pvar.zst"
        return [raw + ".pgen", pvar, raw + ".psam"]
    return None


def check_imputation_prep(params, rep):
    """--mode imputation_prep: raw genotype files, GRCh38 FASTA and (GRCh37 input) liftover chain."""
    g = params.get("genotypes") or {}
    raw = g.get("raw")
    if raw:
        files = _raw_files(str(raw))
        if files is None:
            rep.err("imputation_prep.raw", f"{raw} is neither a VCF/BCF file nor a PLINK prefix "
                    "(.bed/.bim/.fam or .pgen/.pvar/.psam)")
        else:
            missing = [f for f in files if not os.path.exists(f)]
            if missing:
                rep.err("imputation_prep.raw", f"missing file(s): {missing}")
            else:
                rep.ok("imputation_prep.raw", ", ".join(os.path.basename(f) for f in files))
    refs = params.get("refs_dir") or "/opt/xqtl/refs"
    fa = params.get("genome_fasta") or os.path.join(refs, "genome", "GRCh38_no_alt_analysis_set.fa.gz")
    if os.path.exists(fa):
        rep.ok("imputation_prep.genome_fasta", fa)
    else:
        rep.err("imputation_prep.genome_fasta", f"GRCh38 FASTA not found: {fa} (set --genome_fasta or refs_dir)")
    if g.get("raw_build", "GRCh38") == "GRCh37":
        chain = params.get("liftover_chain") or os.path.join(refs, "liftover", "hg19ToHg38.over.chain.gz")
        if os.path.exists(chain):
            rep.ok("imputation_prep.liftover_chain", chain)
        else:
            rep.err("imputation_prep.liftover_chain", f"hg19ToHg38 chain not found: {chain} "
                    "(GRCh37 input is lifted to GRCh38; set --liftover_chain or refs_dir)")
    if g.get("source") == "wes":
        rep.warn("imputation_prep.source", "WES genotypes cover few off-target sites; TOPMed imputation "
                 "accuracy will be low outside exons. Array or WGS genotypes are preferred")


def check_central(params, rep):
    """--mode central: studies table (export folders with MANIFEST.tsv and I_strata), optional GWAS catalog
    and annotation tables. Staged files are matched by basename, so basenames must be unique."""
    c = params.get("central") or {}
    sp = c.get("studies")
    if sp and os.path.exists(sp):
        cols, rows = read_table(sp, "\t")
        if not {"study_id", "export_dir"} <= set(cols):
            rep.err("central.studies", f"{sp} needs columns study_id, export_dir (found {cols})")
        else:
            ids = [r["study_id"] for r in rows]
            dup = sorted({i for i in ids if ids.count(i) > 1})
            if dup:
                rep.err("central.studies", f"duplicate study_id: {dup}")
            if len(rows) < 1:
                rep.err("central.studies", "no studies listed")
            for r in rows:
                d = resolve(r["export_dir"], os.path.dirname(os.path.abspath(sp)))
                if not os.path.isdir(d):
                    rep.err(f"central.studies.{r['study_id']}", f"export_dir not found: {d}")
                elif not os.path.exists(os.path.join(d, "MANIFEST.tsv")):
                    rep.err(f"central.studies.{r['study_id']}", f"no MANIFEST.tsv in {d} (give the export folder "
                            "that holds MANIFEST.tsv)")
                elif not any(f.endswith("_cohort_I_strata.tsv") for f in os.listdir(d)):
                    rep.err(f"central.studies.{r['study_id']}", f"no *_cohort_I_strata.tsv in {d}")
                else:
                    rep.ok(f"central.studies.{r['study_id']}", d)
            if len(rows) == 1:
                rep.warn("central.studies", "only one study: the meta-analysis reduces to within-study "
                         "pooling across ancestries")
    elif sp:
        rep.err("central.studies", f"file not found: {sp}")
    for key, need in (("gwas_catalog", ["trait_id", "path", "ancestry", "n", "trait_type"]),
                      ("annotations", ["name", "bed"])):
        tp = c.get(key)
        if not tp:
            continue
        if not os.path.exists(tp):
            rep.err(f"central.{key}", f"file not found: {tp}")
            continue
        cols, rows = read_table(tp, "\t")
        miss = [x for x in need if x not in cols]
        if miss:
            rep.err(f"central.{key}", f"{tp} is missing column(s) {miss}")
            continue
        pcol = "path" if key == "gwas_catalog" else "bed"
        paths = [resolve(r[pcol], os.path.dirname(os.path.abspath(tp))) for r in rows]
        bases = [os.path.basename(x) for x in paths]
        dup = sorted({b for b in bases if bases.count(b) > 1})
        if dup:
            rep.err(f"central.{key}", f"file basenames must be unique (files are staged by name): {dup}")
        missing = [x for x in paths if not os.path.exists(x)]
        if missing:
            rep.err(f"central.{key}", f"missing file(s): {missing[:5]}")
        if key == "gwas_catalog":
            bad = [r["trait_type"] for r in rows if r["trait_type"] not in ("quantitative", "cc")]
            if bad:
                rep.err("central.gwas_catalog", f"trait_type must be quantitative or cc (got {sorted(set(bad))})")
            for r in rows:
                if r["trait_type"] == "cc":
                    try:
                        ok = 0 < float(r.get("prop_cases") or "nan") < 1
                    except ValueError:
                        ok = False
                    if not ok:
                        rep.err("central.gwas_catalog", f"{r['trait_id']}: trait_type cc needs prop_cases in (0, 1)")
            ids = [r["trait_id"] for r in rows]
            if len(set(ids)) != len(ids):
                rep.err("central.gwas_catalog", "duplicate trait_id")
        if not rep.failed:
            rep.ok(f"central.{key}", f"{len(rows)} entries")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", required=True, help="JSON dump of Nextflow params")
    ap.add_argument("--schema", required=True)
    ap.add_argument("--mode", required=True, choices=sorted(MODE_REQUIRES))
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--no-file-check", action="store_true", help="skip FASTQ existence checks")
    a = ap.parse_args(argv)

    with open(a.params) as fh:
        params = json.load(fh)
    with open(a.schema) as fh:
        schema = json.load(fh)
    params = fill_defaults(schema, params)
    rep = Report()

    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(params), key=lambda e: list(e.path))
    for e in errors:
        loc = ".".join(str(p) for p in e.path) or "(top level)"
        rep.err(f"schema.{loc}", e.message)
    if not errors:
        rep.ok("schema", "cohort parameters match the schema")

    required = list(MODE_REQUIRES[a.mode])
    if a.mode == "site" and params.get("phenotypes_dir"):
        required.remove("site_refs")
        check_phenotypes_dir(params, rep)
    elif a.mode in ("site", "phenotyping_only") and not a.no_file_check:
        check_site_refs(params, rep)
    if a.mode == "imputation_prep":
        check_imputation_prep(params, rep)
    if a.mode == "central":
        check_central(params, rep)
    mods = (params.get("mapping") or {}).get("modalities") or []
    if any(m in mods for m in ISOFORM_MODALITIES) and "expression" not in mods:
        rep.err("mapping.modalities", "isoform_expression / isoform_usage need 'expression' in mapping.modalities "
                "(expression connectivity outliers are removed from the isoform modalities)")
    for key in required:
        val = get_path(params, key)
        if val is None:
            msg = f"--mode {a.mode} needs '{key}' in cohort.yml"
            if key == "restricted":
                msg += ". " + RESTRICTED_HINT
            rep.err(f"required.{key}", msg)
        elif key in ("samplesheet", "metadata", "restricted") and not os.path.exists(val):
            msg = f"'{key}' path does not exist: {val}"
            if key == "restricted":
                msg += ". " + RESTRICTED_HINT
            rep.err(f"required.{key}", msg)
        elif key == "restricted" and not _restricted_ok(val, params, rep):
            pass
        elif key == "genotypes.imputed_dir" and not os.path.isdir(val):
            rep.err(f"required.{key}", f"directory not found: {val}")
        else:
            rep.ok(f"required.{key}", str(val))

    sample_ids = []
    if params.get("samplesheet") and os.path.exists(params["samplesheet"]):
        sample_ids = check_samplesheet(params["samplesheet"], params, rep, check_files=not a.no_file_check)
    if params.get("metadata") and os.path.exists(params["metadata"]):
        check_metadata(params["metadata"], params, rep, sample_ids)

    rep.write(a.out_report)
    with open(a.out_json, "w") as fh:
        json.dump(params, fh, indent=2, sort_keys=True)
    for check, status, detail in rep.rows:
        if status != "PASS":
            print(f"[{status}] {check}: {detail}", file=sys.stderr)
    if rep.failed:
        print("\nValidation FAILED. Fix the problems above and rerun (see docs/CONFIGURATION.md).", file=sys.stderr)
        return 1
    print(f"Validation passed ({sum(1 for r in rep.rows if r[1] == 'PASS')} checks).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

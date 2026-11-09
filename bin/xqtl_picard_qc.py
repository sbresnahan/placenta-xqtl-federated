#!/usr/bin/env python3
"""xqtl_picard_qc.py -- PicardTools sequencing-QC metrics for the HCP prior matrix.

Picard QC metrics and subject-specific bias covariates per sample, in four subcommands that
Nextflow calls:

  generate-refflat  --gtf GTF[.gz] --output refFlat.txt
  gene-gc-length    --gtf GTF[.gz] --fasta genome.fa --output gene_gc_length.tsv
  collect           --bam BAM --sample ID --refflat refFlat.txt --fasta genome.fa
                    --outdir RAW_DIR --output ID.qcmetrics.tsv [--java-mem 4g]
  merge             --rows *.qcmetrics.tsv --salmon-dir DIR --samples samples.txt
                    --gtf GTF[.gz] --gene-annot gene_gc_length.tsv --output qc_metrics.tsv

Panel (devBrain-style): every numeric field of CollectAlignmentSummaryMetrics,
CollectInsertSizeMetrics, CollectRnaSeqMetrics (REF_FLAT, STRAND=NONE), CollectGcBiasMetrics
(summary table) and MarkDuplicates, emitted as `AlignMetrics.` / `InsertSize.` / `RnaMetrics.` /
`GcBias.` / `DupMetrics.` + FIELD. Multi-row tables emit the PAIR (else UNPAIRED) row
unsuffixed and other categories suffixed. merge adds SubjectBias.GC and
SubjectBias.GENE_LENGTH (Mostafavi et al. 2013: per-sample Pearson correlation of
log2(TPM+1) with gene exonic GC content / exonic length) and fills missing cells with
cohort column medians. Zero-variance and redundant
columns are pruned downstream by xqtl_hcp.R.
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xqtl_gtf import _open  # noqa: E402


# =============================================================================
# PicardTools metric parsing
# =============================================================================

def parse_picard_metrics(metrics_file):
    """Parse a PicardTools metrics file into a dict of {column: value}.

    PicardTools metrics files have a header section, then a METRICS CLASS line,
    then a column header line, then one or more data rows. We extract the first
    data row (most tools produce one row per BAM).
    """
    if not os.path.exists(metrics_file):
        return {}

    with open(metrics_file) as f:
        lines = f.readlines()

    metrics_start = None
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            metrics_start = i
            break

    if metrics_start is None:
        return {}

    if metrics_start + 2 >= len(lines):
        return {}

    header = lines[metrics_start + 1].strip().split("\t")
    data = lines[metrics_start + 2].strip().split("\t")

    if len(data) < len(header):
        data.extend([""] * (len(header) - len(data)))

    return dict(zip(header, data))


def coerce_numeric(val):
    """Try to convert a string to float; return None if not numeric."""
    if val is None or val == "" or val == "?":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# =============================================================================
# Full-panel extraction (devBrain-style)
# =============================================================================
# The original pipeline extracted a 24-metric whitelist. The devBrain schema
# (picard_merge_metrics.R) instead keeps the full MultiQC-style field set,
# dropping only identifier/categorical columns and letting downstream
# zero-variance + |r|>0.9 pruning handle redundancy. We mirror that here:
# every numeric field of each tool's summary table is emitted as
# `Tool.FIELD`, with the exclusions below.

QC_EXCLUDE_FIELDS = {
    "SAMPLE", "LIBRARY", "READ_GROUP", "CATEGORY",
    "ACCUMULATION_LEVEL", "READS_USED", "PAIR_ORIENTATION",
}

# Picard >= 2.27 renamed these fields; map the new raw names back to the
# column names the pipeline has always emitted (applied to the primary row
# only, so the original 24 column names are preserved).
FIELD_RENAME = {
    "PF_READS_ALIGNED": "READS_ALIGNED",
    "PCT_PF_READS_ALIGNED": "PCT_READS_ALIGNED",
}


def parse_picard_metrics_all(metrics_file):
    """Parse a PicardTools metrics file into a list of {column: value} dicts,
    one per data row. Multi-row tables (e.g. AlignmentSummaryMetrics) carry
    one row per CATEGORY. Parsing stops at the end of the metrics table
    (blank line or the next ## section, e.g. a HISTOGRAM_CLASS block)."""
    if not os.path.exists(metrics_file):
        return []

    with open(metrics_file) as f:
        lines = f.readlines()

    metrics_start = None
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            metrics_start = i
            break

    if metrics_start is None or metrics_start + 2 >= len(lines):
        return []

    header = lines[metrics_start + 1].strip().split("\t")
    rows = []
    for line in lines[metrics_start + 2:]:
        if not line.strip() or line.startswith("#"):
            break
        data = line.strip().split("\t")
        if len(data) < len(header):
            data.extend([""] * (len(header) - len(data)))
        rows.append(dict(zip(header, data)))
    return rows


def extract_full_panel(rows, prefix):
    """Flatten ALL numeric fields of a Picard metrics table into
    `prefix.FIELD` columns (devBrain-style full panel).

    - ID/categorical columns in QC_EXCLUDE_FIELDS are dropped.
    - The primary row (CATEGORY=PAIR when present, else UNPAIRED, else the
      first row) is emitted unsuffixed; additional rows (e.g.
      FIRST_OF_PAIR / SECOND_OF_PAIR, or extra LIBRARY rows) are suffixed
      with their category value.
    - FIELD_RENAME maps Picard >= 2.27 names back to the pipeline's
      metric column names on the primary row.

    For multi-category tables the primary row is the library-level PAIR row,
    so AlignMetrics.* values are library-level.
    """
    if not rows:
        return {}
    cats = [r.get("CATEGORY", "") for r in rows]
    primary_idx = 0
    for want in ("PAIR", "UNPAIRED"):
        if want in cats:
            primary_idx = cats.index(want)
            break
    multi = len(rows) > 1
    out = {}
    for i, row in enumerate(rows):
        suffix = ""
        if multi and i != primary_idx:
            cat = row.get("CATEGORY") or row.get("LIBRARY") or f"ROW{i + 1}"
            suffix = f".{cat}"
        for field, val in row.items():
            if field in QC_EXCLUDE_FIELDS:
                continue
            num = coerce_numeric(val)
            if num is None:
                continue
            if not suffix:
                field = FIELD_RENAME.get(field, field)
            out[f"{prefix}.{field}{suffix}"] = num
    return out


# =============================================================================
# PicardTools Collect* tool runners
# =============================================================================

def run_picard(cmd, label, sample):
    """Run a PicardTools command, capturing stderr. Returns True on success."""
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600, check=False
        )
        if result.returncode != 0:
            sys.stderr.write(
                f"  WARN [{label}] {sample}: picard exited {result.returncode}\n"
            )
            if result.stderr:
                err_lines = result.stderr.strip().split("\n")
                for line in err_lines:
                    sys.stderr.write(f"        {line}\n")
            return False
        return True
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"  WARN [{label}] {sample}: timed out (3600s)\n")
        return False
    except Exception as e:
        sys.stderr.write(f"  WARN [{label}] {sample}: {e}\n")
        return False


def collect_alignment_summary(picard_cmd, bam, out_dir, sample, parse_only=False):
    """CollectAlignmentSummaryMetrics -> alignment_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "alignment_metrics.txt")
    cmd = [*picard_cmd, "CollectAlignmentSummaryMetrics", f"I={bam}", f"O={out_file}"]
    if parse_only or run_picard(cmd, "AlignSummary", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "AlignMetrics")
    return {}


def collect_insert_size(picard_cmd, bam, out_dir, sample, parse_only=False):
    """CollectInsertSizeMetrics -> insert_size_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "insert_size_metrics.txt")
    hist_file = os.path.join(out_dir, "insert_size_histogram.pdf")
    cmd = [*picard_cmd, "CollectInsertSizeMetrics",
           f"I={bam}", f"O={out_file}", f"H={hist_file}"]
    if parse_only or run_picard(cmd, "InsertSize", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "InsertSize")
    return {}


def collect_rna_seq_metrics(picard_cmd, bam, out_dir, sample, refflat,
                            parse_only=False):
    """CollectRnaSeqMetrics -> rna_metrics.txt (full panel)"""
    out_file = os.path.join(out_dir, "rna_metrics.txt")
    cmd = [*picard_cmd, "CollectRnaSeqMetrics",
           f"I={bam}", f"O={out_file}", f"REF_FLAT={refflat}", "STRAND=NONE"]
    if parse_only or run_picard(cmd, "RnaSeqMetrics", sample):
        return extract_full_panel(parse_picard_metrics_all(out_file), "RnaMetrics")
    return {}


def collect_gc_bias(picard_cmd, bam, out_dir, sample, fasta=None,
                    parse_only=False):
    """CollectGcBiasMetrics -> gc_bias_metrics.txt (+ summary)

    R= (reference FASTA) is a REQUIRED argument for CollectGcBiasMetrics;
    without it the tool errors out and no GcBias columns are produced.
    The dropout/coverage summaries live in the S= summary metrics file,
    not the per-GC-bin detail table (O=).
    """
    out_file = os.path.join(out_dir, "gc_bias_metrics.txt")
    summary_file = os.path.join(out_dir, "gc_bias_summary.txt")
    chart_file = os.path.join(out_dir, "gc_bias_chart.pdf")
    cmd = [*picard_cmd, "CollectGcBiasMetrics",
           f"I={bam}", f"O={out_file}", f"S={summary_file}", f"CHART={chart_file}"]
    if fasta:
        cmd.append(f"R={fasta}")
    if parse_only or run_picard(cmd, "GcBias", sample):
        return extract_full_panel(parse_picard_metrics_all(summary_file), "GcBias")
    return {}


def mark_duplicates(picard_cmd, bam, out_dir, sample, parse_only=False):
    """MarkDuplicates -> dup_metrics.txt (metrics only; temp BAM discarded)"""
    out_file = os.path.join(out_dir, "dup_metrics.txt")
    tmp_bam = os.path.join(out_dir, "tmp_dedup.bam")
    cmd = [*picard_cmd, "MarkDuplicates",
           f"I={bam}", f"O={tmp_bam}", f"M={out_file}",
           "ASSUME_SORTED=true", "VALIDATION_STRINGENCY=LENIENT"]
    if parse_only or run_picard(cmd, "MarkDup", sample):
        result = extract_full_panel(parse_picard_metrics_all(out_file), "DupMetrics")
        if os.path.exists(tmp_bam):
            os.remove(tmp_bam)
        return result
    if os.path.exists(tmp_bam):
        os.remove(tmp_bam)
    return {}


# =============================================================================
# Subject-specific GC and gene-length bias (from HCP paper)
# =============================================================================

def compute_subject_specific_bias(tpm_matrix, gene_annot):
    """Compute subject-specific GC bias and gene-length bias.

    Following Mostafavi et al. 2013:
      - GC bias: Pearson corr(log2(TPM+1), gene GC content) per sample
      - Gene-length bias: Pearson corr(log2(TPM+1), gene length) per sample
    """
    common_genes = tpm_matrix.index.intersection(gene_annot.index)
    if len(common_genes) == 0:
        sys.stderr.write(
            "  WARN: no common genes between TPM matrix and gene annotation\n")
        # Keep 'sample' as the index (not a column) so the downstream join
        # aligns and no duplicate 'sample' column appears in the output
        return pd.DataFrame(
            columns=["SubjectBias.GC", "SubjectBias.GENE_LENGTH"]
        ).rename_axis("sample")

    tpm_sub = tpm_matrix.loc[common_genes]
    gc = gene_annot.loc[common_genes, "gc_content"].values
    length = gene_annot.loc[common_genes, "gene_length"].values
    log_tpm = np.log2(tpm_sub.values + 1)

    results = []
    for i, sample in enumerate(tpm_sub.columns):
        expr = log_tpm[:, i]
        gc_bias = np.corrcoef(expr, gc)[0, 1] if np.std(expr) > 0 and np.std(gc) > 0 else np.nan
        len_bias = np.corrcoef(expr, length)[0, 1] if np.std(expr) > 0 and np.std(length) > 0 else np.nan
        results.append({"sample": sample, "SubjectBias.GC": gc_bias,
                        "SubjectBias.GENE_LENGTH": len_bias})
    return pd.DataFrame(results).set_index("sample")


def load_tx_to_gene(gtf_path):
    """Parse transcript_id -> gene_id mapping from GTF transcript features.

    Needed because HPLRv2/SQANTI3 transcript IDs do not share a parseable
    prefix with gene IDs, so heuristic string-splitting of Salmon Names
    produces gene IDs that never match the GTF gene_id space.
    """
    tx_to_gene = {}
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', fields[8]):
                attr_dict[m.group(1)] = m.group(2)
            tid = attr_dict.get("transcript_id")
            gid = attr_dict.get("gene_id")
            if tid and gid:
                tx_to_gene[tid] = gid
    return tx_to_gene


def load_tpm_matrix(salmon_dir, samples, tx_to_gene=None):
    """Load gene-level TPM matrix from per-sample Salmon quant.sf files.

    Aggregates transcript-level TPM to gene-level by summing transcripts per gene.
    If tx_to_gene (from load_tx_to_gene) is given, transcript Names are mapped
    exactly to GTF gene_ids; otherwise falls back to a heuristic
    (Name prefix before the first '.' or '__').
    """
    tpm_dfs = []
    for sample in samples:
        quant_file = os.path.join(salmon_dir, sample, "quant.sf")
        if not os.path.exists(quant_file):
            sys.stderr.write(f"  WARN: quant.sf not found for {sample}: {quant_file}\n")
            continue
        df = pd.read_csv(quant_file, sep="\t")
        if tx_to_gene:
            df["gene_id"] = df["Name"].map(tx_to_gene)
            df = df.dropna(subset=["gene_id"])
        else:
            df["gene_id"] = df["Name"].str.split(".").str[0]
            df["gene_id"] = df["gene_id"].str.split("__").str[0]
        gene_tpm = df.groupby("gene_id")["TPM"].sum()
        gene_tpm.name = sample
        tpm_dfs.append(gene_tpm)

    if not tpm_dfs:
        return pd.DataFrame()

    tpm_matrix = pd.concat(tpm_dfs, axis=1)
    tpm_matrix.index.name = "gene_id"
    return tpm_matrix


# =============================================================================
# refFlat generation from GTF
# =============================================================================

def generate_refflat(gtf_path, output_path):
    """Convert a GTF file to refFlat format for PicardTools CollectRnaSeqMetrics.

    refFlat format (tab-delimited, one line per transcript, 11 fields):
      geneName name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds

    Coordinates are 0-based. exonStarts and exonEnds are comma-separated lists
    ending with a comma.
    """
    genes = {}  # gene_name -> {transcript_id -> tx dict}
    tx_to_gene = {}  # transcript_id -> gene_name

    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            start = int(start) - 1  # GTF 1-based -> refFlat 0-based
            end = int(end)

            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', attrs):
                attr_dict[m.group(1)] = m.group(2)

            gene_name = attr_dict.get("gene_name") or attr_dict.get("gene_id", "")
            transcript_id = attr_dict.get("transcript_id", "")

            if feature == "transcript" and transcript_id:
                if gene_name not in genes:
                    genes[gene_name] = {}
                genes[gene_name][transcript_id] = {
                    "chrom": chrom, "strand": strand,
                    "txStart": start, "txEnd": end,
                    "cdsStart": start, "cdsEnd": end, "exons": [],
                }
                tx_to_gene[transcript_id] = gene_name

            elif feature == "exon" and transcript_id:
                gname = tx_to_gene.get(transcript_id)
                if gname and transcript_id in genes.get(gname, {}):
                    genes[gname][transcript_id]["exons"].append((start, end))

            elif feature == "CDS" and transcript_id:
                gname = tx_to_gene.get(transcript_id)
                if gname and transcript_id in genes.get(gname, {}):
                    tx = genes[gname][transcript_id]
                    if tx["cdsStart"] == tx["txStart"]:
                        tx["cdsStart"] = start
                    tx["cdsEnd"] = end

    with open(output_path, "w") as out:
        for gene_name in sorted(genes):
            txs = genes[gene_name]
            for tx_id in sorted(txs):
                tx = txs[tx_id]
                if not tx["exons"]:
                    continue
                exons = sorted(tx["exons"], key=lambda x: x[0])
                exon_count = len(exons)
                exon_starts = ",".join(str(e[0]) for e in exons) + ","
                exon_ends = ",".join(str(e[1]) for e in exons) + ","
                out.write("\t".join([
                    gene_name, tx_id, tx["chrom"], tx["strand"],
                    str(tx["txStart"]), str(tx["txEnd"]),
                    str(tx["cdsStart"]), str(tx["cdsEnd"]),
                    str(exon_count), exon_starts, exon_ends,
                ]) + "\n")

    n_genes = len(genes)
    n_tx = sum(len(txs) for txs in genes.values())
    print(f"refFlat written: {output_path} ({n_genes} genes, {n_tx} transcripts)")


# =============================================================================
# Gene GC/length annotation table generation
# =============================================================================

def generate_gene_gc_length(gtf_path, fasta_path, output_path):
    """Compute per-gene GC content and gene length from GTF + genome FASTA.

    Gene length = sum of exon lengths (union of merged exons per gene).
    GC content = fraction of G/C bases in the gene's exonic sequence.
    """
    try:
        import pysam
    except ImportError:
        sys.stderr.write(
            "ERROR: pysam required for --generate-gene-annot. "
            "Install with: uv pip install pysam\n")
        sys.exit(1)

    gene_exons = defaultdict(list)  # gene_id -> [(chrom, start, end)]

    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature != "exon":
                continue
            start = int(start) - 1
            end = int(end)

            attr_dict = {}
            for m in re.finditer(r'(\w+)\s+"([^"]*)"', attrs):
                attr_dict[m.group(1)] = m.group(2)
            gene_id = attr_dict.get("gene_id", "")
            if not gene_id:
                continue
            gene_exons[gene_id].append((chrom, start, end))

    print(f"Loaded {len(gene_exons)} genes from GTF")

    fasta = pysam.FastaFile(fasta_path)
    fasta_chroms = set(fasta.references)

    results = []
    n_done = 0
    for gene_id in sorted(gene_exons):
        exons = gene_exons[gene_id]
        by_chrom = defaultdict(list)
        for chrom, start, end in exons:
            by_chrom[chrom].append((start, end))

        total_length = 0
        total_gc = 0

        for chrom, intervals in by_chrom.items():
            if chrom not in fasta_chroms:
                continue
            intervals.sort()
            merged = []
            for start, end in intervals:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            for start, end in merged:
                seq = fasta.fetch(chrom, start, end).upper()
                total_length += len(seq)
                total_gc += seq.count("G") + seq.count("C")

        gc_content = total_gc / total_length if total_length > 0 else np.nan
        results.append({"gene_id": gene_id, "gc_content": gc_content,
                        "gene_length": total_length})

        n_done += 1
        if n_done % 2000 == 0:
            print(f"  Processed {n_done}/{len(gene_exons)} genes")

    df = pd.DataFrame(results).set_index("gene_id")
    df.to_csv(output_path, sep="\t")
    print(f"Gene GC/length annotation written: {output_path} ({len(df)} genes)")



# =============================================================================
# CLI
# =============================================================================

def cmd_collect(args):
    os.makedirs(args.outdir, exist_ok=True)
    picard = ["picard", f"-Xmx{args.java_mem}"] if args.java_mem else ["picard"]
    row = {}
    row.update(collect_alignment_summary(picard, args.bam, args.outdir, args.sample))
    if args.paired:
        row.update(collect_insert_size(picard, args.bam, args.outdir, args.sample))
    row.update(collect_rna_seq_metrics(picard, args.bam, args.outdir, args.sample,
                                       args.refflat))
    row.update(collect_gc_bias(picard, args.bam, args.outdir, args.sample,
                               fasta=args.fasta))
    row.update(mark_duplicates(picard, args.bam, args.outdir, args.sample))
    if not row:
        sys.exit(f"ERROR: no Picard metrics collected for {args.sample}")
    df = pd.DataFrame([row], index=pd.Index([args.sample], name="sample"))
    df.to_csv(args.output, sep="\t")
    print(f"{args.sample}: {df.shape[1]} metrics -> {args.output}")


def cmd_merge(args):
    rows = []
    for path in args.rows:
        df = pd.read_csv(path, sep="\t", index_col=0)
        if len(df) != 1:
            sys.exit(f"ERROR: {path} has {len(df)} rows (expected 1)")
        rows.append(df)
    merged = pd.concat(rows)
    merged.index = merged.index.astype(str)
    merged.index.name = "sample"
    if merged.index.duplicated().any():
        sys.exit("ERROR: duplicate samples in Picard rows: "
                 + ",".join(merged.index[merged.index.duplicated()]))

    with open(args.samples) as fh:
        samples = [s.strip() for s in fh if s.strip()]
    missing = sorted(set(samples) - set(merged.index))
    if missing:
        sys.exit(f"ERROR: Picard rows missing for {len(missing)} samples: {missing[:5]}")
    merged = merged.loc[samples]

    gene_annot = pd.read_csv(args.gene_annot, sep="\t", index_col=0)
    tpm = load_tpm_matrix(args.salmon_dir, samples, load_tx_to_gene(args.gtf))
    bias = compute_subject_specific_bias(tpm, gene_annot)
    merged = merged.join(bias, how="left")

    num_cols = merged.select_dtypes(include=[np.number]).columns
    n_missing = int(merged[num_cols].isna().sum().sum())
    if n_missing:
        medians = merged[num_cols].median().fillna(0.0)
        merged[num_cols] = merged[num_cols].fillna(medians)
        sys.stderr.write(f"  WARN: {n_missing} missing cells filled with column medians\n")
    merged.to_csv(args.output, sep="\t")
    print(f"qc_metrics: {merged.shape[0]} samples x {merged.shape[1]} metrics -> {args.output}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("generate-refflat")
    p.add_argument("--gtf", required=True)
    p.add_argument("--output", required=True)

    p = sub.add_parser("gene-gc-length")
    p.add_argument("--gtf", required=True)
    p.add_argument("--fasta", required=True, help="indexed (samtools faidx) genome FASTA")
    p.add_argument("--output", required=True)

    p = sub.add_parser("collect")
    p.add_argument("--bam", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--refflat", required=True)
    p.add_argument("--fasta", required=True,
                   help="genome FASTA with .fai and .dict (CollectGcBiasMetrics R=)")
    p.add_argument("--outdir", required=True, help="directory for raw Picard outputs")
    p.add_argument("--output", required=True)
    p.add_argument("--paired", action="store_true",
                   help="run CollectInsertSizeMetrics (paired-end libraries only)")
    p.add_argument("--java-mem", default=None, help="Java heap, e.g. 4g")

    p = sub.add_parser("merge")
    p.add_argument("--rows", nargs="+", required=True)
    p.add_argument("--salmon-dir", required=True, help="contains <sample>/quant.sf")
    p.add_argument("--samples", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--gene-annot", required=True)
    p.add_argument("--output", required=True)

    args = ap.parse_args()
    if args.cmd == "generate-refflat":
        generate_refflat(args.gtf, args.output)
    elif args.cmd == "gene-gc-length":
        generate_gene_gc_length(args.gtf, args.fasta, args.output)
    elif args.cmd == "collect":
        cmd_collect(args)
    else:
        cmd_merge(args)


if __name__ == "__main__":
    main()

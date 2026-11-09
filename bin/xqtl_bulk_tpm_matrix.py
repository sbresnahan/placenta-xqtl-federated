#!/usr/bin/env python3
"""Script 2: Build a gene-symbol x sample TPM matrix from Salmon quant.sf files.

Aggregates per-sample Salmon transcript TPMs to gene-level TPMs using the
gene_name attribute from the normalized HPLRv2 GTF (produced by prep_PANTRY_gtf.R).
This matrix is the bulk mixture input for MuSiC deconvolution.

Using gene symbols (not Ensembl gene_id) maximizes overlap with the Campbell
et al. 2023 scRNA-seq reference (GSE182381), which uses GENCODE gene symbols.

Inputs:
  - Salmon quant.sf files: <salmon_dir>/<sample>/quant.sf
    Columns: Name, Length, EffectiveLength, TPM, NumReads
  - Normalized HPLRv2 GTF (has gene_name + transcript_id attributes on
    transcript feature lines, added by prep_PANTRY_gtf.R)
  - Samples file: one sample ID per line

Output:
  - TSV: genes (symbols) x samples, TPM values, non-log linear space

Usage:
  python3 build_bulk_tpm_matrix.py \
    --gtf        /path/to/HPLRv2.0.annotated.PANTRY.gtf \
    --salmon-dir /path/to/cohort/salmon_expression \
    --samples    /path/to/samples.txt \
    --output     /path/to/cohort_bulk_tpm_symbols.tsv

"""

import argparse
import sys
import os
from pathlib import Path

import pandas as pd
import numpy as np

from xqtl_gtf import _open

# ---------------------------------------------------------------------------
# GTF parsing — extract transcript_id -> gene_name map
# ---------------------------------------------------------------------------

def build_tx_to_gene_name_map(gtf_path: Path) -> pd.DataFrame:
    """Parse GTF and build a transcript_id -> gene_name mapping table.

    Uses the 'transcript' feature lines, which carry both transcript_id and
    (after prep_PANTRY_gtf.R normalization) gene_name attributes.

    For transcripts missing gene_name, falls back to gene_id with the version
    suffix stripped (e.g. ENSG00000123456.1 -> ENSG00000123456).

    Returns a DataFrame with columns: transcript_id, gene_name
    """
    records = []
    with _open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'transcript':
                continue
            attrs = fields[8]
            # Parse the GTF attribute column (key "value"; pairs separated by ';')
            attr_dict = {}
            for pair in attrs.split(';'):
                pair = pair.strip()
                if not pair:
                    continue
                # Handle both: key "value"  and  key "value";
                parts = pair.split(' ', 1)
                if len(parts) != 2:
                    continue
                key = parts[0]
                val = parts[1].strip().strip('"')
                attr_dict[key] = val
            tx_id = attr_dict.get('transcript_id')
            gene_name = attr_dict.get('gene_name')
            gene_id = attr_dict.get('gene_id')
            if tx_id is None:
                continue
            if gene_name is None or gene_name == '':
                # Fallback: strip version suffix from gene_id
                if gene_id is not None:
                    gene_name = gene_id.split('.')[0]
                else:
                    continue
            records.append({'transcript_id': tx_id, 'gene_name': gene_name})

    df = pd.DataFrame(records)
    # Drop duplicate transcript_ids (keep first)
    df = df.drop_duplicates(subset='transcript_id', keep='first')
    print(f"  GTF: {len(df)} transcript_id -> gene_name mappings")
    if df['gene_name'].duplicated().any():
        n_dup = df['gene_name'].duplicated().sum()
        print(f"  GTF: {n_dup} transcripts share gene_name with another transcript (expected)")
    return df


# ---------------------------------------------------------------------------
# Salmon quant.sf aggregation
# ---------------------------------------------------------------------------

def aggregate_sample_tpm(quant_sf: Path, tx_map: pd.DataFrame) -> pd.Series:
    """Read one sample's quant.sf and aggregate transcript TPMs to gene-level.

    Sums TPMs across all transcripts mapping to the same gene_name.
    Returns a Series indexed by gene_name with TPM values.
    """
    df = pd.read_csv(quant_sf, sep='\t')
    # Salmon quant.sf: Name, Length, EffectiveLength, TPM, NumReads
    # 'Name' is the transcript ID
    df = df[['Name', 'TPM']].rename(columns={'Name': 'transcript_id', 'TPM': 'tpm'})
    # Merge with transcript->gene map
    df = df.merge(tx_map, on='transcript_id', how='inner')
    # Aggregate: sum TPM across transcripts per gene
    gene_tpm = df.groupby('gene_name')['tpm'].sum()
    return gene_tpm


def aggregate_multirun_tpm(quant_sfs: list, tx_map: pd.DataFrame) -> pd.Series:
    """Collapse technical-replicate runs into one gene-level TPM profile.

    Sums per-transcript NumReads across the runs, recomputes TPM from the
    summed counts (rate = NumReads / EffectiveLength, TPM = rate / sum(rate)
    * 1e6), then aggregates to gene level. TPMs themselves are never
    averaged. EffectiveLength is taken from the first run (identical across
    runs of the same cohort/protocol).
    """
    summed = None
    eff_len = None
    for q in quant_sfs:
        df = pd.read_csv(q, sep='\t')
        s = df.set_index('Name')['NumReads'].astype(float)
        if eff_len is None:
            eff_len = df.set_index('Name')['EffectiveLength'].astype(float)
        summed = s if summed is None else summed.add(s, fill_value=0.0)
    rate = summed / eff_len
    tpm = rate / rate.sum() * 1e6
    df = tpm.rename('tpm').reset_index().rename(columns={'Name': 'transcript_id'})
    df = df.merge(tx_map, on='transcript_id', how='inner')
    return df.groupby('gene_name')['tpm'].sum()


def build_bulk_matrix(salmon_dir: Path, samples: list, tx_map: pd.DataFrame,
                      collapse_map: dict = None) -> pd.DataFrame:
    """Build genes x samples TPM matrix from per-sample Salmon quant.sf files.

    Args:
        salmon_dir: Directory containing <sample>/quant.sf subdirectories.
        samples: List of sample IDs to process.
        tx_map: DataFrame with transcript_id, gene_name columns.
        collapse_map: Optional dict {sample_id: individual_id}. Samples
            sharing an individual_id are technical replicates and are
            collapsed to one column (summed counts, recomputed TPM) named
            by the individual_id.

    Returns:
        DataFrame: rows = gene_name, columns = sample IDs (or individual
        IDs when collapse_map is given), values = TPM.
    """
    if collapse_map:
        # Group samples by individual, preserving first-seen order
        ind_to_samples = {}
        for s in samples:
            ind_to_samples.setdefault(collapse_map.get(s, s), []).append(s)
        n_pairs = sum(1 for v in ind_to_samples.values() if len(v) > 1)
        print(f"  Collapse map: {len(samples)} samples -> "
              f"{len(ind_to_samples)} individuals ({n_pairs} replicated)")
        work_items = list(ind_to_samples.items())
    else:
        work_items = [(s, [s]) for s in samples]

    matrices = []
    missing = []
    for i, (name, group) in enumerate(work_items):
        quant_paths = [salmon_dir / s / 'quant.sf' for s in group]
        present = [q for q in quant_paths if q.exists()]
        if not present:
            missing.append(name)
            continue
        if len(present) < len(quant_paths):
            print(f"  WARNING: {name}: only {len(present)}/{len(quant_paths)} "
                  f"run quant.sf files present")
        if len(present) == 1:
            gene_tpm = aggregate_sample_tpm(present[0], tx_map)
        else:
            gene_tpm = aggregate_multirun_tpm(present, tx_map)
        gene_tpm.name = name
        matrices.append(gene_tpm)
        if (i + 1) % 50 == 0 or (i + 1) == len(work_items):
            print(f"  Processed {i+1}/{len(work_items)} samples")

    if missing:
        print(f"  WARNING: {len(missing)} samples missing quant.sf: {missing[:5]}{'...' if len(missing)>5 else ''}")

    if not matrices:
        raise ValueError(f"No quant.sf files found in {salmon_dir}")

    # Outer join to keep all genes (fill missing with 0)
    bulk = pd.concat(matrices, axis=1)
    bulk = bulk.fillna(0.0)

    # Handle duplicate gene names (sum them — standard for gene-level aggregation)
    if bulk.index.duplicated().any():
        n_before = len(bulk)
        bulk = bulk.groupby(level=0).sum()
        print(f"  Summed duplicate gene names: {n_before} -> {len(bulk)} unique genes")

    return bulk


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build gene-symbol x sample TPM matrix from Salmon quant.sf for MuSiC deconvolution"
    )
    parser.add_argument('--gtf', type=str, default=None,
                        help='Path to normalized HPLRv2 GTF (with gene_name attribute)')
    parser.add_argument('--salmon-dir', type=str, default=None,
                        help='Directory containing <sample>/quant.sf subdirectories')
    parser.add_argument('--samples', type=str, default=None,
                        help='File with one sample ID per line')
    parser.add_argument('--output', type=str, required=True,
                        help='Output TSV path (genes x samples, TPM values)')
    parser.add_argument('--collapse-map', type=str, default=None,
                        help='Optional TSV with rnaseq_id, array_id columns. '
                             'Runs sharing an array_id are technical '
                             'replicates and are collapsed to one column '
                             '(summed counts, recomputed TPM) per individual.')
    args = parser.parse_args()

    if not (args.gtf and args.salmon_dir and args.samples):
        parser.error("--gtf, --salmon-dir and --samples are required")
    gtf_path = Path(args.gtf)
    salmon_dir = Path(args.salmon_dir)
    samples_file = Path(args.samples)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"=== Building bulk TPM matrix ===")
    print(f"GTF:        {gtf_path}")
    print(f"Salmon dir: {salmon_dir}")
    print(f"Samples:    {samples_file}")
    print(f"Output:     {output_path}")
    print()

    # --- Load samples ---
    with open(samples_file) as fh:
        samples = [line.strip() for line in fh if line.strip()]
    print(f"Samples to process: {len(samples)}")

    # --- Build transcript -> gene_name map ---
    print(f"\nParsing GTF for transcript_id -> gene_name map...")
    tx_map = build_tx_to_gene_name_map(gtf_path)

    # --- Optional replicate collapse map ---
    collapse_map = None
    if args.collapse_map:
        cm = pd.read_csv(args.collapse_map, sep='\t')
        if not {'rnaseq_id', 'array_id'} <= set(cm.columns):
            parser.error("--collapse-map TSV must have rnaseq_id and array_id columns")
        collapse_map = dict(zip(cm['rnaseq_id'], cm['array_id']))
        print(f"Collapse map: {args.collapse_map} ({len(collapse_map)} runs)")

    # --- Build bulk matrix ---
    print(f"\nAggregating Salmon TPMs to gene level...")
    bulk = build_bulk_matrix(salmon_dir, samples, tx_map, collapse_map=collapse_map)
    print(f"\nBulk matrix: {bulk.shape[0]} genes x {bulk.shape[1]} samples")

    # --- Report gene overlap with reference ---
    # The reference (GSE182381) uses GENCODE gene symbols. We report overlap
    # if a reference gene list is available alongside.
    ref_genes_path = output_path.parent / 'reference_gene_symbols.txt'
    if ref_genes_path.exists():
        ref_genes = set(open(ref_genes_path).read().splitlines())
        bulk_genes = set(bulk.index)
        shared = bulk_genes & ref_genes
        print(f"\nGene overlap with reference:")
        print(f"  Bulk genes:    {len(bulk_genes)}")
        print(f"  Reference genes: {len(ref_genes)}")
        print(f"  Shared:        {len(shared)} ({100*len(shared)/len(ref_genes):.1f}% of reference)")
    else:
        print(f"\n(No reference_gene_symbols.txt found — skipping overlap report)")
        print(f"  To enable: extract gene symbols from GSE182381_reference_sample.txt.gz")
        print(f"  and save as {ref_genes_path}")

    # --- Write output ---
    print(f"\nWriting {output_path}...")
    # Write TSV: first column = gene symbol, then sample columns
    bulk_out = bulk.copy()
    bulk_out.index.name = 'gene_symbol'
    bulk_out.to_csv(output_path, sep='\t', float_format='%.6g')
    print(f"Done. Size: {output_path.stat().st_size / 1e6:.1f} MB")
    print(f"\n=== Bulk TPM matrix complete ===")


if __name__ == '__main__':
    main()

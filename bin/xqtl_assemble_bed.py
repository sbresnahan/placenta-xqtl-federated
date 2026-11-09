#!/usr/bin/env python3
"""xqtl_assemble_bed.py -- assemble per-sample quantifications into phenotype BED files.

Modified from Pantry's assemble_bed.py (github.com/PejLab/Pantry, MIT licence):
  * Salmon quant.sf input in place of kallisto;
  * GTF parsing uses bin/xqtl_gtf.py (no gtfparse dependency);
  * splicing clusters are keyed by the full cluster token, so frozen clusters (clu_N_strand, from
    the restricted file) and site-local novel clusters (nclu_N_strand) never collide;
  * intron retention reads IRFinder-S per-sample tables (see docs/PHENOTYPING.md);
  * no LaDDR 'latent' subcommand.
"""

import argparse
from pathlib import Path
from xqtl_gtf import read_gtf
import numpy as np
import pandas as pd

def load_tss(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations from GTF file

    Returns TSS as the first four columns of the BED format, meaning the
    coordinates are 0-based and chromEnd is just chromStart + 1.
    """
    anno = read_gtf(ref_anno, features=['gene'])
    anno['chromEnd'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['chromStart'] = anno['chromEnd'] - 1 # BED coordinates are 0-based
    anno['#chrom'] = anno['seqname']
    anno = anno.sort_values(['#chrom', 'chromStart'])
    anno = anno[['#chrom', 'chromStart', 'chromEnd', 'gene_id']]
    # Rename columns for tensorQTL:
    anno.columns = ['#chr', 'start', 'end', 'gene_id']
    return anno

def load_tss_with_gene_id_base(ref_anno: Path) -> pd.DataFrame:
    """Load TSS annotations with versionless gene IDs for external tools."""
    anno = load_tss(ref_anno)
    anno['gene_id_base'] = anno['gene_id'].astype(str).str.split('.').str[0]
    return anno

def load_exons(ref_anno: Path) -> pd.DataFrame:
    """Load exon annotations

    Returns exons with start and end oriented on the gene's strand.
    """
    anno = read_gtf(ref_anno, features=['exon'])
    anno['chrom'] = anno['seqname']
    anno['exonStart'] = np.where(anno['strand'] == '+', anno['start'], anno['end'])
    anno['exonEnd'] = np.where(anno['strand'] == '+', anno['end'], anno['start'])
    return anno[['gene_id', 'chrom', 'exonStart', 'exonEnd']]

def transcript_to_gene_map(ref_anno: Path) -> pd.DataFrame:
    """Load transcript IDs and corresponding gene IDs from GTF file"""
    anno = read_gtf(ref_anno, features=['transcript'])
    return anno[['gene_id', 'transcript_id']]

def map_introns_to_genes(introns: list, exons: pd.DataFrame) -> pd.DataFrame:
    """Map de novo splice junctions to genes based on known exon boundaries"""
    exons['exonStart'] = exons['exonStart'].astype(str)
    exons['exonEnd'] = exons['exonEnd'].astype(str)
    df = pd.DataFrame({'intron': introns})
    # Expected intron format: chrom:start:end:clu_ _
    # Some chromosome IDs contain underscores (e.g. "NC_000001.11"), so split
    # first by colon and then split the last field by underscores.
    df[['chrom', 'chr_start', 'chr_end', 'cluster_info']] = df['intron'].str.split(':', expand=True)
    # cluster_info is the full cluster token: clu_N_strand (frozen) or nclu_N_strand (novel)
    df['cluster'] = df['cluster_info']
    df['strand'] = df['cluster_info'].str[-1]
    df['start'] = np.where(df['strand'] == '+', df['chr_start'], df['chr_end'])
    df['end'] = np.where(df['strand'] == '+', df['chr_end'], df['chr_start'])
    start_matches = df.merge(
        exons[['chrom', 'exonEnd', 'gene_id']].rename(columns={'exonEnd': 'start'}),
        on=['chrom', 'start'],
        how='inner'
    )
    end_matches = df.merge(
        exons[['chrom', 'exonStart', 'gene_id']].rename(columns={'exonStart': 'end'}),
        on=['chrom', 'end'],
        how='inner'
    )
    df = pd.concat([start_matches, end_matches])
    clust_genes = df.groupby('cluster', group_keys=False).agg({'gene_id': pd.Series.unique})
    clust_genes = clust_genes.reset_index().explode('gene_id')
    return clust_genes

# =============================================================================
# Salmon reader (Pantry's original reads kallisto abundance.tsv)
# =============================================================================
# Original kallisto version read abundance.tsv with:
#   pd.read_csv(fname, sep='\t', index_col='target_id')
#   d = d[units]  # 'tpm' or 'est_counts'
#
# Salmon quant.sf has columns: Name, Length, EffectiveLength, TPM, NumReads
# We map: Name -> target_id (index), TPM -> tpm, NumReads -> est_counts,
# EffectiveLength -> eff_length
# =============================================================================

# Map Salmon quant.sf column names to the kallisto-compatible names used
# by the rest of the script:
SALMON_COLUMN_MAP = {
    'Name': 'target_id',
    'TPM': 'tpm',
    'NumReads': 'est_counts',
    'EffectiveLength': 'eff_length',
}

def load_salmon(sample_ids: list, salmon_dir: Path, units: str) -> pd.DataFrame:
    """Assemble Salmon TPM or NumReads outputs into a table.

    Reads quant.sf from each sample's output directory.

    Args:
        sample_ids: List of sample IDs.
        salmon_dir: Directory containing per-sample subdirectories with quant.sf.
        units: 'tpm' or 'est_counts' (mapped from Salmon 'TPM' or 'NumReads'),
            or 'tpm_from_counts': length-normalized TPM recomputed per sample as
            (NumReads/EffectiveLength) / sum(NumReads/EffectiveLength) * 1e6.
            Use 'tpm_from_counts' for the QU-corrected quants (expression_qu/):
            qu_correct_salmon.R writes count-proportional CPM values into the
            'TPM' column there (drop-in quant.sf compatibility), so the 'tpm'
            mode would read CPM-scale values, not length-normalized TPM.
    """
    counts = []
    for i, sample in enumerate(sample_ids):
        fname = salmon_dir / sample / 'quant.sf'
        d = pd.read_csv(fname, sep='\t')
        # Rename Salmon columns to kallisto-compatible names
        d = d.rename(columns=SALMON_COLUMN_MAP)
        d = d.set_index('target_id')
        if units == 'tpm_from_counts':
            eff = pd.to_numeric(d['eff_length'], errors='coerce')
            cnt = pd.to_numeric(d['est_counts'], errors='coerce').fillna(0.0)
            x = cnt / eff
            x[(eff.isna()) | (eff <= 0)] = 0.0  # undefined length -> no abundance
            x = x.fillna(0.0)
            total = x.sum()
            d = x / total * 1e6 if total > 0 else x
        else:
            d = d[units]  # 'tpm' or 'est_counts'
        d.name = sample
        counts.append(d)
    return pd.concat(counts, axis=1)


def assemble_alt_TSS_polyA(sample_ids: list, group1_dir: Path, group2_dir: Path, units: str, ref_anno: Path, bed: Path, min_frac: float = 0.05, max_frac: float = 0.95, skip_feature_filter: bool = False):
    """Assemble txrevise-based Salmon outputs into BED file

    By default, txrevise produces two sets of annotations per annotation type.
    Relative TSS/polyA site usage should be calculated within each set, and
    can then be combined to produce a single BED file.
    """
    df1 = load_salmon(sample_ids, group1_dir, units)
    df2 = load_salmon(sample_ids, group2_dir, units)

    # This assumes target_id is e.g. {gene_id}.grp_1.downstream.{transcript_id}
    df1['gene_id'] = df1.index.str.split('.grp_').str[0]
    df2['gene_id'] = df2.index.str.split('.grp_').str[0]

    # Calculate proportion of each transcript in each gene_id:
    gene_ids = df1['gene_id'] # groupby/apply removes gene_id, so save it to add back
    df1 = df1.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    if skip_feature_filter:
        # Union pooling: keep every quantified site; drop only all-zero rows.
        # The pooled devBrain detection filter decides at normalization.
        df1 = df1[df1.fillna(0).sum(axis=1) > 0]
    else:
        # Remove sites with mean relative usage < `min_frac` or > `max_frac`:
        df1 = df1[(df1.mean(axis=1) >= min_frac) & (df1.mean(axis=1) <= max_frac)]
    df1 = df1.join(gene_ids, how='left')

    gene_ids = df2['gene_id']
    df2 = df2.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    if skip_feature_filter:
        df2 = df2[df2.fillna(0).sum(axis=1) > 0]
    else:
        # Remove sites with mean relative usage < `min_frac` or > `max_frac`:
        df2 = df2[(df2.mean(axis=1) >= min_frac) & (df2.mean(axis=1) <= max_frac)]
    df2 = df2.join(gene_ids, how='left')

    df = pd.concat([df1, df2], axis=0)
    # Use __ to separate gene ID from other info in phenotype IDs
    df.index = df.index.str.replace('.grp_', '__grp_', 1)
    # Replace period separators for easier parsing in case of periods in transcript IDs
    df.index = df.index.str.replace('.upstream.', '_upstream_', 1)
    df.index = df.index.str.replace('.downstream.', '_downstream_', 1)
    df.index = df.index.rename('phenotype_id')
    df = df.reset_index()

    # Use gene's TSS for all of its isoforms:
    anno = load_tss(ref_anno)
    df = anno.merge(df.reset_index(), on='gene_id', how='inner')
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def assemble_expression(sample_ids: list, salmon_dir: Path, units: str, ref_anno: Path, bed_iso: Path, bed_gene: Path, min_count: int = 10, min_frac: float = 0.05, max_frac: float = 0.95, log2_expr: bool = False, bed_iso_expr: Path = None, skip_feature_filter: bool = False):
    """Assemble Salmon TPM or NumReads outputs into isoform- and gene-level BED files

    Isoform values are normalized to relative abundance in each gene. Isoforms
    with fewer than `min_count` reads on average are excluded
    (`est_counts` read counts are always used for this filtering).

    `bed_iso`, `bed_gene`, and `bed_iso_expr` are each optional (may be None):
    only the non-None BEDs are written. This lets the caller assemble the
    gene-level BED from one Salmon directory (e.g. original counts) and the
    isoform BEDs from another (e.g. QU-corrected counts) in two separate calls.
    At least one must be set.

    `bed_iso_expr` is the isoform *expression* (abundance) BED: the same
    per-transcript values as `bed_iso` but BEFORE the within-gene ratio
    division, filtered only by the `min_count` floor (the devBrain-style
    TPM > 0.1 in > 25% of samples filter is applied downstream, at the
    normalization stage). `bed_iso` remains the isoform *usage* (ratio) BED.

    Units note: for QU-corrected quants (expression_qu/), pass
    units='tpm_from_counts'. The adjusted quant.sf 'TPM' column is CPM-scale
    (count-proportional; see qu_correct_salmon.R), so 'tpm' would yield
    read-fraction abundances/ratios; 'tpm_from_counts' recomputes
    length-normalized TPM from NumReads x EffectiveLength, giving
    molecule-fraction abundances and usage ratios (PANTRY/devBrain
    convention).
    """
    if bed_iso is None and bed_gene is None and bed_iso_expr is None:
        raise ValueError("assemble_expression: at least one of bed_iso / bed_gene / bed_iso_expr must be provided")

    df_iso = load_salmon(sample_ids, salmon_dir, units)

    # Record isoforms with any value >= `min_count` to keep:
    if units == 'est_counts':
        df_counts = df_iso
    else:
        df_counts = load_salmon(sample_ids, salmon_dir, 'est_counts')
    if skip_feature_filter:
        # Union pooling: keep every isoform with any reads; drop only
        # all-zero rows. The pooled devBrain detection filter decides at
        # normalization.
        iso_enough_counts = df_counts[df_counts.fillna(0).sum(axis=1) > 0].index
    else:
        iso_enough_counts = df_counts[df_counts.mean(axis=1) >= min_count].index

    df_iso.index = df_iso.index.rename('transcript_id')
    gene_map = transcript_to_gene_map(ref_anno).set_index('transcript_id')
    df_iso = df_iso.join(gene_map, how='inner')
    assert df_iso.shape[0] > 0, 'No matching isoforms in Salmon output and reference annotation'

    # Also get gene-level expression:
    df_gene = df_iso.reset_index().drop('transcript_id', axis=1).groupby('gene_id', group_keys=True).sum()
    if log2_expr:
        df_gene = np.log2(df_gene + 1)

    # Capture isoform expression (abundance) BEFORE the within-gene ratio
    # division; same min_count floor as the ratio BED, no min/max-frac filter.
    if bed_iso_expr is not None:
        df_iso_expr = df_iso[df_iso.index.isin(iso_enough_counts)]

    # Calculate proportion of each transcript in each gene_id:
    df_iso = df_iso.groupby('gene_id', group_keys=False).apply(lambda x: x / x.sum(axis=0))
    # Remove isoforms with mean read count < `min_count`:
    df_iso = df_iso[df_iso.index.isin(iso_enough_counts)]
    if not skip_feature_filter:
        # Remove isoforms with mean relative abundance < `min_frac` or > `max_frac`:
        df_iso = df_iso[(df_iso.mean(axis=1) >= min_frac) & (df_iso.mean(axis=1) <= max_frac)]
    # groupby/apply removes gene_id, so add them back
    df_iso = df_iso.join(gene_map, how='left')

    anno = load_tss(ref_anno)
    if bed_iso is not None:
        df_iso_out = anno.merge(df_iso.reset_index(), on='gene_id', how='inner')
        df_iso_out['phenotype_id'] = df_iso_out['gene_id'] + '__' + df_iso_out['transcript_id']
        df_iso_out = df_iso_out[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
        df_iso_out.to_csv(bed_iso, sep='\t', index=False, float_format='%g')

    if bed_iso_expr is not None:
        df_expr_out = anno.merge(df_iso_expr.reset_index(), on='gene_id', how='inner')
        df_expr_out['phenotype_id'] = df_expr_out['gene_id'] + '__' + df_expr_out['transcript_id']
        df_expr_out = df_expr_out[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
        df_expr_out.to_csv(bed_iso_expr, sep='\t', index=False, float_format='%g')

    if bed_gene is not None:
        df_gene_out = anno.merge(df_gene.reset_index(), on='gene_id', how='inner')
        df_gene_out = df_gene_out.rename(columns={'gene_id': 'phenotype_id'})
        df_gene_out = df_gene_out[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
        df_gene_out.to_csv(bed_gene, sep='\t', index=False, float_format='%g')

def assemble_RNA_editing(data: Path, ref_anno: Path, bed: Path):
    """Convert prepared RNA editing phenotype matrix into BED file."""
    df = pd.read_csv(data, sep='\t')
    sample_ids = list(df.columns[1:])
    df = df.rename(columns={df.columns[0]: 'phenotype_id'})
    df['gene_id'] = df['phenotype_id'].str.replace(r'__.*$', '', regex=True)

    anno = load_tss(ref_anno)
    df = anno.merge(df, on='gene_id', how='inner')
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

# full cluster token at the end of a leafCutter intron ID chrom:start:end:<token>
CLUSTER_TOKEN = r':(n?clu_\d+_[+\-?])$'


def assemble_splicing(counts: Path, ref_anno: Path, bed: Path, min_frac: float = 0.05, max_frac: float = 0.95, skip_feature_filter: bool = False):
    """Convert leafcutter output into splicing BED file"""
    df = pd.read_csv(counts, sep=' ')
    sample_ids = list(df.columns)
    cluster = df.index.str.extract(CLUSTER_TOKEN, expand=False)
    df = df.groupby(cluster, group_keys=False).apply(lambda g: g / g.sum(axis=0))
    if skip_feature_filter:
        # Union pooling: keep every quantified junction; drop only all-zero
        # rows. The pooled devBrain detection filter decides at normalization.
        df = df[df.fillna(0).sum(axis=1) > 0]
    else:
        # Remove junctions with mean relative usage < `min_frac` or > `max_frac`:
        df = df[(df.mean(axis=1) >= min_frac) & (df.mean(axis=1) <= max_frac)]
    df['cluster'] = df.index.str.extract(CLUSTER_TOKEN, expand=False)
    df.index = df.index.rename('intron')
    df = df.reset_index()

    exons = load_exons(ref_anno)
    genes = map_introns_to_genes(df['intron'], exons)
    df = df.merge(genes, on='cluster', how='left')

    anno = load_tss(ref_anno)
    df = anno.merge(df, on='gene_id', how='inner')
    df['intron'] = df['intron'].str.replace(':', '_')
    df['phenotype_id'] = df['gene_id'] + '__' + df['intron']
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

IRF_MIN_DEPTH = 10


def load_irfinder(sample_ids: list, ir_dir: Path, min_depth: int = IRF_MIN_DEPTH) -> pd.DataFrame:
    """IRratio of 'clean' introns from IRFinder-S <ir_dir>/<sample>.IRFinder-IR-*.txt.

    The ratio is set to NaN when IntronDepth + max(SpliceLeft, SpliceRight) < min_depth.
    Rows are indexed by (gene_id, chrom, start, end, strand); IRFinder intron starts are 0-based.
    """
    cols = []
    for sample in sample_ids:
        # stranded libraries: IRFinder-S writes IRFinder-IR-dir.txt (used) next to -nondir.txt
        hits = [ir_dir / f'{sample}.IRFinder-IR-{k}.txt' for k in ('dir', 'nondir')]
        hits = [h for h in hits if h.exists()]
        if not hits:
            raise FileNotFoundError(f'no IRFinder-IR-dir.txt / -nondir.txt for {sample} in {ir_dir}')
        d = pd.read_csv(hits[0], sep='\t', dtype={'Chr': str})
        d = d.rename(columns={'Chr': 'chrom', 'Start': 'start', 'End': 'end', 'Strand': 'strand'})
        parts = d['Name'].str.split('/', expand=True)
        d['gene_id'] = parts[1]
        d = d[parts[2] == 'clean']
        depth = d['IntronDepth'] + d[['SpliceLeft', 'SpliceRight']].max(axis=1)
        val = d['IRratio'].astype(float).where(depth >= min_depth)
        val.index = pd.MultiIndex.from_frame(d[['gene_id', 'chrom', 'start', 'end', 'strand']])
        val = val[~val.index.duplicated()]
        val.name = sample
        cols.append(val)
    return pd.concat(cols, axis=1)


def assemble_intron_retention(sample_ids: list, ir_dir: Path, ref_anno: Path, bed: Path):
    """IRFinder-S intron retention ratios -> BED (phenotype GENE__IR_chrom_start_end_strand)."""
    df = load_irfinder(sample_ids, ir_dir)
    df = df[df.notnull().any(axis=1)].reset_index()
    df['phenotype_id'] = (df['gene_id'] + '__IR_' + df['chrom'].astype(str) + '_' + df['start'].astype(str)
                          + '_' + df['end'].astype(str) + '_' + df['strand'].astype(str))
    anno = load_tss(ref_anno)
    df = anno.merge(df.drop(columns=['chrom', 'start', 'end', 'strand']), on='gene_id', how='inner')
    df = df[['#chr', 'start', 'end', 'phenotype_id'] + sample_ids]
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def load_featureCounts(sample_ids: list, counts_dir: Path, feature: str, min_count: int = 10) -> pd.DataFrame:
    """Assemble featureCounts outputs into a table"""
    counts = []
    for i, sample in enumerate(sample_ids):
        fname_ex = counts_dir / f'{sample}.{feature}.counts.txt'
        d = pd.read_csv(fname_ex, sep='\t', index_col='Geneid', skiprows=1)
        d = d.iloc[:, 5]
        d.name = sample
        d[d < min_count] = np.nan
        counts.append(d)
    return pd.concat(counts, axis=1)

def assemble_stability(sample_ids: list, stab_dir: Path, ref_anno: Path, bed: Path, skip_feature_filter: bool = False):
    """Assemble exon to intron read ratios into mRNA stability BED file"""
    exon = load_featureCounts(sample_ids, stab_dir, 'exonic')
    intron = load_featureCounts(sample_ids, stab_dir, 'intronic')
    genes = exon.index[np.isin(exon.index, intron.index)]
    assert exon.loc[genes, :].index.equals(intron.loc[genes, :].index)
    assert exon.columns.equals(intron.columns)
    df = exon.loc[genes, :] / intron.loc[genes, :]
    if skip_feature_filter:
        # Union pooling: drop only features with no quantification at all;
        # the pooled >=40% detection filter decides at normalization.
        df = df[df.notnull().any(axis=1)]
    else:
        df = df[df.isnull().mean(axis=1) <= 0.5]

    anno = load_tss(ref_anno)
    anno = anno.rename(columns={'gene_id': 'phenotype_id'})
    df.index = df.index.rename('phenotype_id')
    df = anno.merge(df.reset_index(), on='phenotype_id', how='inner')
    df.to_csv(bed, sep='\t', index=False, float_format='%g')

def _load_samples(path: Path) -> list[str]:
    return pd.read_csv(path, sep='\t', header=None, dtype=str)[0].tolist()

def main():
    parser = argparse.ArgumentParser(description='Assemble data into an RNA phenotype BED file')
    sub = parser.add_subparsers(dest='cmd', required=True)

    # alt_TSS_polyA
    p_tss = sub.add_parser('alt-tss-polya', help='Assemble txrevise alt TSS or polyA usage')
    p_tss.add_argument('--samples', type=Path, required=True, help='Path to sample IDs file (one sample ID per line)')
    p_tss.add_argument('--group1-dir', type=Path, required=True, help='Group 1 directory (path of grp_1.upstream for TSS, grp_1.downstream for polyA) containing Salmon output directories named by sample')
    p_tss.add_argument('--group2-dir', type=Path, required=True, help='Group 2 directory (path of grp_2.upstream for TSS, grp_2.downstream for polyA) containing Salmon output directories named by sample')
    p_tss.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_tss.add_argument('--output', type=Path, required=True, help='Output BED file')
    p_tss.add_argument('--min-frac', type=float, default=0.05, help='Minimum mean fraction for TSS/polyA site to be included')
    p_tss.add_argument('--max-frac', type=float, default=0.95, help='Maximum mean fraction for TSS/polyA site to be included')
    p_tss.add_argument('--skip-feature-filter', action='store_true',
                        help=('Union-pooling mode: emit all quantified features (drop only '
                              'all-zero/all-NaN rows) and let the pooled devBrain detection '
                              'filters at the normalization stage decide the final feature set'))

    # expression
    p_expr = sub.add_parser('expression', help='Assemble isoform- and gene-level expression')
    p_expr.add_argument('--samples', type=Path, required=True, help='Path to sample IDs file (one sample ID per line)')
    p_expr.add_argument('--input-dir', type=Path, required=True, help='Directory containing Salmon output directories named by sample')
    p_expr.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_expr.add_argument('--output-isoforms', type=Path, default=None, help='Isoform ratio (usage) BED (optional; at least one of --output-isoforms/--output-expression/--output-isoform-expr is required)')
    p_expr.add_argument('--output-expression', type=Path, default=None, help='Gene-level expression BED (optional; at least one of --output-isoforms/--output-expression/--output-isoform-expr is required)')
    p_expr.add_argument('--output-isoform-expr', dest='output_isoform_expr', type=Path, default=None, help='Isoform expression (abundance, pre-ratio) BED (optional)')
    p_expr.add_argument('--min-count', type=int, default=10, help='Minimum mean count for isoform to be included in isoform-level BED')
    p_expr.add_argument('--min-frac', type=float, default=0.05, help='Minimum mean relative abundance for isoform to be included in isoform-level BED')
    p_expr.add_argument('--max-frac', type=float, default=0.95, help='Maximum mean relative abundance for isoform to be included in isoform-level BED')
    p_expr.add_argument('--units', choices=['tpm', 'est_counts', 'tpm_from_counts'], default='tpm',
                        help="Units to read from Salmon quant.sf (tpm, est_counts, or tpm_from_counts). "
                             "Defaults to tpm. Use tpm_from_counts for QU-corrected quants (expression_qu/): "
                             "their 'TPM' column is CPM-scale, so length-normalized TPM is recomputed from "
                             "NumReads and EffectiveLength instead")
    p_expr.add_argument('--log2-expr', action='store_true',
                        help='If set, gene-level expression values are log2 transformed (log2(x + 1)) before being written to BED file')
    p_expr.add_argument('--skip-feature-filter', action='store_true',
                        help=('Union-pooling mode: emit all quantified features (drop only '
                              'all-zero/all-NaN rows) and let the pooled devBrain detection '
                              'filters at the normalization stage decide the final feature set'))

    # RNA editing
    p_edit = sub.add_parser('rna-editing', help='Assemble RNA editing phenotypes')
    p_edit.add_argument('--input', type=Path, required=True, help='Prepared RNA editing phenotype matrix')
    p_edit.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_edit.add_argument('--output', type=Path, required=True, help='Output BED file')

    # splicing
    p_sp = sub.add_parser('splicing', help='Assemble leafcutter splicing phenotypes')
    p_sp.add_argument('--input', type=Path, required=True, help='leafcutter per-junction counts')
    p_sp.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_sp.add_argument('--output', type=Path, required=True, help='Output BED file')
    p_sp.add_argument('--min-frac', type=float, default=0.05, help='Minimum mean fraction for junction to be included')
    p_sp.add_argument('--max-frac', type=float, default=0.95, help='Maximum mean fraction for junction to be included')
    p_sp.add_argument('--skip-feature-filter', action='store_true',
                        help=('Union-pooling mode: emit all quantified features (drop only '
                              'all-zero/all-NaN rows) and let the pooled devBrain detection '
                              'filters at the normalization stage decide the final feature set'))

    # intron retention
    p_ir = sub.add_parser('intron-retention', help='Assemble IRFinder-S intron retention ratios')
    p_ir.add_argument('--samples', type=Path, required=True, help='Path to sample IDs file (one sample ID per line)')
    p_ir.add_argument('--input-dir', type=Path, required=True, help='Directory of <sample>.IRFinder-IR-dir.txt (or -nondir) files')
    p_ir.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_ir.add_argument('--output', type=Path, required=True, help='Output BED file')

    # stability
    p_stab = sub.add_parser('stability', help='Assemble mRNA stability (exon/intron ratios)')
    p_stab.add_argument('--samples', type=Path, required=True, help='Path to sample IDs file (one sample ID per line)')
    p_stab.add_argument('--input-dir', type=Path, required=True, help='Directory containing {sample}.exonic.counts.txt and {sample}.intronic.counts.txt files')
    p_stab.add_argument('--ref-anno', dest='ref_anno', type=Path, required=True, help='Reference annotation GTF file')
    p_stab.add_argument('--output', type=Path, required=True, help='Output BED file')
    p_stab.add_argument('--skip-feature-filter', action='store_true',
                        help=('Union-pooling mode: emit all quantified features (drop only '
                              'all-zero/all-NaN rows) and let the pooled devBrain detection '
                              'filters at the normalization stage decide the final feature set'))

    args = parser.parse_args()

    if args.cmd == 'alt-tss-polya':
        samples = _load_samples(args.samples)
        assemble_alt_TSS_polyA(samples, args.group1_dir, args.group2_dir, 'tpm',
                               args.ref_anno, args.output,
                               min_frac=args.min_frac, max_frac=args.max_frac,
                               skip_feature_filter=args.skip_feature_filter)
    elif args.cmd == 'expression':
        samples = _load_samples(args.samples)
        if args.output_isoforms is None and args.output_expression is None and args.output_isoform_expr is None:
            parser.error("expression: at least one of --output-isoforms / --output-expression / --output-isoform-expr is required")
        assemble_expression(samples, args.input_dir, args.units,
                            args.ref_anno, args.output_isoforms, args.output_expression,
                            min_count=args.min_count,
                            min_frac=args.min_frac, max_frac=args.max_frac,
                            log2_expr=args.log2_expr,
                            bed_iso_expr=args.output_isoform_expr,
                            skip_feature_filter=args.skip_feature_filter)
    elif args.cmd == 'intron-retention':
        assemble_intron_retention(_load_samples(args.samples), args.input_dir, args.ref_anno, args.output)
    elif args.cmd == 'rna-editing':
        assemble_RNA_editing(args.input, args.ref_anno, args.output)
    elif args.cmd == 'splicing':
        assemble_splicing(args.input, args.ref_anno, args.output,
                          min_frac=args.min_frac, max_frac=args.max_frac,
                          skip_feature_filter=args.skip_feature_filter)
    elif args.cmd == 'stability':
        samples = _load_samples(args.samples)
        assemble_stability(samples, args.input_dir, args.ref_anno, args.output,
                           skip_feature_filter=args.skip_feature_filter)

if __name__ == '__main__':
    main()

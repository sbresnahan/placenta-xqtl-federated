# Export specification

This page describes the folder a site shares with MDACC after `--mode site` (or `--mode export`). It is written to `<outdir>/export/` and holds only summary statistics. Everything else under `<outdir>/` stays at the site, including `site/`, `pipeline_info/` and the work directory.

Files A–F, I and J are produced from pipeline version 1.0.0 onwards. Files G (cell-type interaction, G×C) and H (exposure interaction, GxE) are added by later milestones. The export format is versioned with the pipeline (see the [changelog](../CHANGELOG.md)).

## Layout

```
export/<cohort_id>/
  MANIFEST.tsv                         J  every file: path, sha256, bytes, rows, class, versions
  DATA_DICTIONARY.md                   J  column definitions
  <cohort_id>_cohort_I_strata.tsv      I  ancestry strata: n individuals, n RNA-seq samples, status
  <ANC>/                                  one folder per mapped ancestry stratum (n >= ancestry.min_n)
    A_nominal/<cohort>_<ANC>_A_nominal_<modality>_<chrom>.parquet
    B_mapping/<cohort>_<ANC>_B_permutation_<arm>.parquet
    B_mapping/<cohort>_<ANC>_B_independent_<arm>.parquet
    C_ld/<cohort>_<ANC>_C_ld_<chrom>_<start>_<end>.npz  (+ .variants.tsv.gz, C_ld.windows.tsv)
    <cohort>_<ANC>_D_variant_qc.parquet
    <cohort>_<ANC>_E_features.parquet
    <cohort>_<ANC>_F_model_settings.json
    <cohort>_<ANC>_I_summary.parquet
```

With `export.format: tsv.gz`, tables are written as gzipped TSV instead of Parquet. Variant IDs are GRCh38 `chr:pos:ref:alt`. Effects are per copy of the ALT allele, in normalized phenotype units.

**Arms.** A modality with feature groups (for example `splicing`, where IDs are `GENE__feature`) is mapped twice:

- grouped, one test per gene, arm `<modality>`;
- ungrouped, one test per feature, arm `<modality>_ungrouped`.

Modalities without groups have one arm. The `combined` arm maps all modalities of a gene together as one group (PANTRY-style).

## File classes

| Class | Content | Key columns |
|---|---|---|
| A | all cis variant–phenotype pairs within `mapping.cis_window` (ungrouped, per modality and chromosome) | phenotype_id, variant_id, chrom, pos, ref, alt, effect_allele, other_allele, tss_distance, eaf, maf, ma_count, n, beta, se, pval |
| B | tensorQTL `map_cis` lead per phenotype or group, with beta-approximated permutation p-values and q-values; conditionally independent signals from `map_independent` | arm, phenotype_id, variant_id, slope, slope_se, pval_nominal, pval_beta, qval, qvalue_method, (group_id, group_size), rank |
| C | in-sample LD: signed r stored as float16, for variants with MAF ≥ `ld.maf_min` and imputation R² ≥ `ld.rsq_min`, in `ld.window_mb` windows with a `ld.step_mb` step. Each `.npz` holds exactly one square matrix `R`; the row order is in the matching `.variants.tsv.gz` | R; variant_id, chrom, pos, ref, alt, eaf |
| D | per-variant QC in the stratum | variant_id, eaf, maf, n, hwe_p, rsq, typed |
| E | phenotype features: coordinates (gene TSS), group, detection fraction, filter status | modality, phenotype_id, chr, start, end, detection_fraction, status, group_id |
| F | model settings: normalization, HCP grid with eGene counts, selected k, covariate **names** (no values), permutations, window, MAF threshold, seed | JSON |
| I | aggregate counts (genotype QC steps, normalization), and 5/25/50/75/95 % quantiles, mean and SD of sequencing metrics and cell-type proportions (computed only when n ≥ 10) | table, item, statistic, value |
| J | manifest with sha256 and the image sha256, restricted-file version and pipeline version; data dictionary | |

## What is never exported

- Sample, RNA-seq or genotype identifiers.
- Per-sample values of any kind: genotypes, dosages, phenotypes, covariates, genotype or HCP PCs, ancestry assignments, outlier lists.
- Minimum or maximum values.

The site keeps these under `<outdir>/site/` for its own records.

## Leak scan

The leak scan runs at the end of every export and cannot be disabled. It collects every identifier from:

- the metadata columns `rnaseq_id` and `genotype_id`;
- the samplesheet `sample_id` column;
- all sample IDs in the imputed VCFs.

It then checks every exported file. The export fails if any of the following is found:

- an identifier in a file name, a column name or any text value. Identifiers are matched as whole tokens, so `S001` matches `reads_S001` but not `S0010`. Purely numeric or very short IDs are matched only as a complete value. Variant IDs are exempt;
- an `.npz` file that holds anything other than one square matrix `R`;
- a file type other than Parquet, TSV, gzipped TSV, NPZ, JSON or Markdown.

When the scan fails, the export folder is deleted and the run stops. The findings are written to `leak_report.tsv` in the task work directory; it is local only and must not be shared. Fix the cause, then rebuild the export with `--mode export`.

## Rebuilding the export

`--mode export` rebuilds `export/` from a finished site run without recomputing anything, for example after changing `export.format` or `export.round_digits`. It reads:

- `<outdir>/site/ancestry/cohort_strata.tsv`;
- the per-stratum folders `<outdir>/site/<ANC>/`;
- `<outdir>/site/genotype_qc/cohort.vcf_samples.txt`;
- `<outdir>/site/restricted_version.txt`.

Any previous `export/` folder is removed first.

## Phenotype inputs (`phenotypes_dir`)

`--mode site` normally builds phenotypes from FASTQs. If `phenotypes_dir` is set, it maps ready inputs instead, and `site_refs` is not required. The directory must contain:

| File | Content |
|---|---|
| `<modality>.bed.gz` | one per modality in `mapping.modalities`. BED columns `#chr start end phenotype_id`, then one column per RNA-seq sample. Coordinates are the gene TSS (`start = TSS - 1`, `end = TSS`), so all features of a gene share one cis window. Grouped modalities use `GENE__feature` IDs |
| `expression.counts.bed.gz` | raw gene counts in the same layout (used for TMM-VST) |
| `qc_metrics.tsv` | `sample` column + sequencing QC metrics (used for HCP and file I) |
| `cell_proportions.tsv` | `sample` column + estimated cell-type proportions (optional) |

## Site-side processing behind the export

| Step | Method (atlas-compatible) |
|---|---|
| Genotype QC | TOPMed R² ≥ `qc.rsq_min`, biallelic. Dosages are kept and hard calls are the nearest genotype. Per-stratum MAF, MAC, HWE and missingness filters, then PC outliers (6 SD on PC1–5) are removed |
| Ancestry | PC-AiR on KING-robust kinship against the 1000 Genomes panel. Each sample gets the nearest superpopulation centroid on PCs 1–10. Samples are "unassigned" if the QDA posterior for that label is below `ancestry.assign_prob_min`. Strata below `ancestry.min_n` are counted but not mapped |
| Normalization | Expression: TMM-VST, then removal of connectivity outliers (bicor z < −3). Other modalities: quantile normalization + inverse normal transform. Technical replicates are averaged after normalization |
| Covariates | Sex, genotype PCs 1–5, cell-type proportions (dominant type dropped, `gxc.exclude_cell_types` removed, arcsinh) and HCP factors. Pruned at \|r\| > 0.9, keeping the higher-priority covariate |
| k selection | The HCP grid is `mapping.hcp_grid`, or by default 0–60 by 10 (n < 350) or 0–100 by 10. Values of k that leave fewer than 20 residual degrees of freedom are skipped. k* gives the most eGenes on `mapping.eval_chrom`, or genome-wide if that chromosome has fewer than `mapping.eval_min_phenotypes` phenotypes. Ties go to the smaller k. The grouped and ungrouped arms share k* |
| Mapping | tensorQTL `map_cis` (`mapping.nperm` permutations, Storey q-values with a BH fallback recorded in `qvalue_method`), `map_independent` in chunks of `mapping.independent_chunk_size`, and `map_nominal` |

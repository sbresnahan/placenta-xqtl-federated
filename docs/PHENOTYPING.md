# Phenotyping

This page describes how a site's RNA-seq FASTQs become the nine molecular phenotypes that are mapped by `--mode site`. The steps follow the placenta xQTL atlas (scripts `02_rnaseq_quantification` and `03_phenotyping`), which builds on [Pantry](https://github.com/PejLab/Pantry) (MIT licence, `LICENSES/PANTRY-MIT.txt`). Places where the federated pipeline differs from the atlas are listed at the end.

Phenotyping runs in two modes:

- `--mode site`, when `phenotypes_dir` is not set: FASTQs → phenotypes → genotype QC, mapping and export.
- `--mode phenotyping_only`: FASTQs → phenotype BEDs and normalized matrices, without genotypes. Nothing is exported. Use it to check RNA-seq processing before genotypes are ready, or to look at the cohort's phenotypes locally.

## Before the first run: site indices

```bash
run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode prepare_references
```

This builds the aligner and quantifier indices once, from the public references embedded in the image, and writes them to `site_refs:`. Expect 1–2 hours; STAR's genome generation needs about 32–40 GB of memory.

| Directory in `site_refs` | Content |
|---|---|
| `genome/` | GRCh38 analysis set FASTA, `.fai`, Picard `.dict` |
| `annotation/` | HPLRv2 annotation (header lines removed), exonic / intronic GTFs for stability, Picard `refFlat.txt`, gene GC content and length |
| `star/` | STAR index with the annotation's junctions, `--sjdbOverhang` = `rnaseq.read_length` − 1 |
| `salmon/expression/` | Salmon index of the annotation's transcripts (`gffread`) |
| `salmon/txrevise.grp_{1,2}.{upstream,downstream}/` | Salmon indices of the txrevise alternative transcript-start and transcript-end events |
| `irfinder/IRFinder/` | IRFinder-S reference (from the STAR index, with mapability exclusion) |
| `editing/` | REDIportal sites on the annotation's chromosomes and their genes (only if the REDIportal table is present) |
| `MANIFEST.tsv` | pipeline version, read length, sha256 of every input, tool versions |

Rebuild the indices if `rnaseq.read_length` changes a lot; the validator warns when the two differ. Alignment still works with a different read length.

## Inputs

- `samplesheet`: CSV with `sample_id,fastq_1,fastq_2`. Leave `fastq_2` empty for single-end data. Relative paths are resolved from the samplesheet's directory. `sample_id` must match `rnaseq_id` in the metadata.
- `rnaseq` in `cohort.yml`: `read_layout` (paired / single), `strandedness` (reverse / forward / unstranded), `read_length` and `library_type` (total / polyA / 3prime).
- `restricted`: the restricted file from MDACC. Phenotyping uses its frozen splicing clusters.

FASTQs are used as delivered: there is no adapter or quality trimming. STAR soft-clips adapter read-through.

## Per-sample steps

| Step | Tool and settings |
|---|---|
| Alignment | STAR 2-pass (`--twopassMode Basic --outSAMstrandField intronMotif`), coordinate-sorted BAM |
| Gene / transcript quantification | Salmon `quant -l A --validateMappings --numBootstraps 20` |
| txrevise events | Salmon against each of the 4 txrevise indices |
| Junctions | regtools `junctions extract -a 8 -m 50 -M 500000 -s XS` |
| Exonic / intronic reads | featureCounts on exonic (`-t exon --fracOverlap 1`) and intronic (`-t intron`) GTFs, strand from `rnaseq.strandedness`, `-p` for paired-end |
| Intron retention | IRFinder-S `BAM` mode |
| RNA editing | per-site editing levels at REDIportal sites (MAPQ 255, base quality ≥ 20) |
| Sequencing QC | Picard CollectAlignmentSummaryMetrics, MarkDuplicates metrics, CollectGcBiasMetrics, CollectInsertSizeMetrics (paired-end) and CollectRnaSeqMetrics |

## Cohort-level phenotypes

| Modality | How it is built |
|---|---|
| `expression` | Salmon transcript counts (bootstrap-based quantification-uncertainty correction) summarised to genes with tximport. Gene TPM goes into `expression.bed.gz`; raw counts into `expression.counts.bed.gz`, which normalization uses for TMM-VST |
| `isoform_expression` | transcript TPM, grouped by gene. Raw transcript counts are in `isoform_expression.counts.bed.gz` |
| `isoform_usage` | transcript TPM / gene TPM, grouped by gene |
| `alt_tss` | txrevise upstream event usage (groups 1 and 2) |
| `alt_polya` | txrevise downstream event usage (groups 1 and 2) |
| `splicing` | leafCutter intron excision ratios in frozen and novel clusters (below) |
| `intron_retention` | IRFinder-S IR ratios of clean introns, grouped by gene |
| `rna_editing` | editing levels at sites covered ≥ 20 reads in at least max(10, 10 % of samples), collapsing sites that correlate above 0.9 |
| `stability` | exonic / intronic read ratio per gene |

All BEDs use the gene TSS as coordinates (`start = TSS − 1`, `end = TSS`), and features of grouped modalities are named `GENE__feature`. `qc_metrics.tsv` has the Picard metrics plus `SubjectBias.GC` and `SubjectBias.GENE_LENGTH` (per-sample Pearson correlation of log2(TPM + 1) with gene exonic GC content and length; Mostafavi et al. 2013). As in the atlas, a metric missing for a sample is replaced by the cohort median.

### Splicing: frozen and novel clusters

Splicing phenotypes are only comparable between studies if every site uses the same leafCutter clusters. The restricted file therefore holds the atlas clusters (`splicing/frozen_clusters.tsv.gz`, IDs `clu_<n>_<strand>`). Each site:

1. counts its junction reads in the frozen clusters (leafCutter with the frozen clusters as the refined cluster file);
2. clusters its junctions locally with the atlas settings (`--maxintronlen 100000 --minclureads 30 --mincluratio 0.001`);
3. keeps local clusters that share no junction with the frozen set, renamed `nclu_<n>_<strand>`.

Phenotype IDs are `GENE__chrom_start_end_clu_<n>_<strand>`; novel clusters carry `nclu_`. Novel-cluster junctions are mapped and exported like the others, and the `novel` column of the splicing features file flags them. Only junctions in frozen clusters can be meta-analysed across studies; novel clusters are cohort-specific. `splicing.clusters_summary.tsv` counts the frozen, matched and novel clusters.

### Cell-type proportions

MuSiC estimates cell-type proportions for each RNA-seq sample from gene TPM, using the placental single-cell reference embedded in the image. `cell_proportions.tsv` has 8 collapsed cell groups (used as covariates and for G×C); `cell_proportions_full.tsv` has every reference cell type, and `maternal_flag.tsv` marks samples whose maternal-cell fraction is above 10 %.

### Library types

`3prime` libraries only produce `expression` and `stability`. polyA libraries produce all modalities, but intronic signal (stability, intron retention) is weaker than in total RNA. `rna_editing` is skipped when `site_refs` has no REDIportal sites.

## Outputs

| Path | Content |
|---|---|
| `<outdir>/phenotypes/` | `<modality>.bed.gz` (+ `.tbi`), the two count BEDs, `qc_metrics.tsv`, cell proportions, `splicing.clusters_summary.tsv`, `PHENOTYPES_SUMMARY.tsv` (features and samples per file) |
| `<outdir>/normalized/` | `--mode phenotyping_only` only: `ALL_<modality>.bed.gz` normalized across all samples (columns are genotype IDs), with `phenotype_groups.txt` for grouped modalities |
| `<outdir>/site/ALL/` | `--mode phenotyping_only` only: normalization diagnostics (outliers, feature lists) |

`<outdir>/phenotypes/` has the layout `phenotypes_dir` expects, so a later `--mode site` run can map it directly. None of these files leave the site.

## Differences from the atlas

| Area | Atlas | Federated pipeline | Why |
|---|---|---|---|
| Intron retention | MAJIQ | IRFinder-S | MAJIQ's licence does not allow redistribution in the image |
| Splicing clusters | clusters from the pooled atlas cohorts | frozen atlas clusters + site-specific novel clusters | sites cannot pool junctions |
| BAMs for regtools / featureCounts | shrunk BAMs | full BAMs | one BAM per sample, fewer steps |
| Salmon expression index | atlas-built index | annotation transcripts from `gffread`, no decoy sequences | rebuilt from public files at each site |
| Picard | in the main environment | separate Java environment in the image | dependency conflict |
| IRFinder-S | stock | memory check lowered (`IRFINDER_MIN_RAM_MB`) | the stock check refuses to run below 32 GB even when far less is needed |
| txrevise sets | upstream / downstream | the same; the `contained` sets are shipped but not indexed | matches the atlas |
| Cell proportions | technical replicates collapsed | one estimate per RNA-seq sample | consistent with the other phenotypes |
| Small feature sets | — | VST dispersion fit on all features when fewer than 1,000 have mean count > 5 | only affects tests and single-chromosome runs |

## Test fixture

`tests/scripts/test_phenotyping.sh` runs `prepare_references` and `phenotyping_only` on 8 public SNUH placenta samples (PRJNA820329), with reads cut to chr22 (`tests/fixtures/phenotyping/`, built by `tests/fixtures/make_fastq_fixture.sh`). It uses GENCODE v45 chr22, processed with Pantry's txrevise steps, as a public stand-in for the HPLRv2 annotation (`tests/fixtures/refs_chr22/`, built by `tests/fixtures/make_refs_fixture.sh`). The test's restricted file freezes every other leafCutter cluster of these samples, so both frozen and novel clusters are exercised (`tests/fixtures/make_phenotyping_restricted.sh`). The fixture FASTQs have constant base qualities (they were rebuilt from BAMs).

On the chr22 fixture, IRFinder-S writes an empty `intergenic.ROI.bed`, because it only keeps intergenic regions on chromosome names of up to two characters. This does not affect IR ratios.

# Configuration

Cohort settings are written in one YAML file, `cohort.yml`, and passed with `-params-file`. Cluster settings (queues, memory, bind paths) belong in a site profile, not here; see [HPC configuration](HPC.md).

```bash
./xqtl/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode <mode>
```

Any setting can also be given on the command line (`--outdir /scratch/run2`), which overrides the file. Nested settings are easier to keep in the file.

Every mode starts by checking `cohort.yml` and the input files against the schema (`assets/schema/cohort.schema.json`). The run stops before any task is submitted if a check fails. The report is written to `<outdir>/pipeline_info/validation_report.tsv`, and the complete settings with every default filled in are written to `resolved_params.json` in the same folder. To run only the checks:

```bash
./xqtl/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode site --validate_only
```

[Troubleshooting](TROUBLESHOOTING.md#input-validation-messages) explains each message.

## A minimal `cohort.yml`

```yaml
cohort_id: COHORT1                       # letters, digits, _ or -; used in file names
outdir: /scratch/xqtl/COHORT1            # visible to compute nodes
samplesheet: samplesheet.csv
metadata: metadata.tsv
restricted: /secure/xqtl-restricted-1.tar.gz
site_refs: /project/xqtl/site_refs       # written once by --mode prepare_references

genotypes:
  raw: /project/geno/array_merged        # PLINK prefix, VCF or BCF (only for --mode imputation_prep)
  raw_build: GRCh37
  source: array
  imputed_dir: /project/geno/topmed_r3   # decrypted TOPMed output (for --mode site)

rnaseq:
  read_layout: paired
  strandedness: reverse
  read_length: 150
  library_type: total
```

Everything else has a default. Relative paths in `cohort.yml` are resolved from the launch directory; relative FASTQ paths in the samplesheet are resolved from the samplesheet's folder.

## Input tables

**Samplesheet** (`samplesheet`, CSV): one row per RNA-seq library.

```text
sample_id,fastq_1,fastq_2
PL001,/data/fastq/PL001_R1.fastq.gz,/data/fastq/PL001_R2.fastq.gz
PL002,/data/fastq/PL002_R1.fastq.gz,/data/fastq/PL002_R2.fastq.gz
```

`fastq_2` is empty for single-end data. Every `sample_id` must appear as an `rnaseq_id` in the metadata.

**Metadata** (`metadata`, TSV): one row per RNA-seq sample.

| Column | Required | Content |
|---|---|---|
| `rnaseq_id` | yes | Matches `sample_id` in the samplesheet |
| `genotype_id` | yes | Sample ID in the imputed VCFs |
| `sex` | yes | Fetal sex: `M`/`F` (or `male`/`female`); missing values are imputed with the most common value |
| `GA` | yes | Gestational age at delivery, weeks, 20–45 |
| `ppBMI`, `gdm`, `ogtt` | no | Exposures for GxE (`gdm` coded 0/1) |
| `maternal_age`, `parity`, `smoking` | no | Sensitivity covariates for GxE |

Other columns are ignored with a warning. None of these values leave the site.

## Settings reference

### Run and inputs

| Key | Default | Meaning |
|---|---|---|
| `cohort_id` | required | Short cohort label used in file names and in the export |
| `mode` | — | `prepare_references`, `imputation_prep`, `site`, `phenotyping_only`, `export`, `central` or `smoke`; usually given as `--mode` |
| `outdir` | `results` | Output folder |
| `release` | `<cohort_id>_<yyyymmdd>` | Export release tag |
| `samplesheet` | — | Samplesheet CSV (above) |
| `metadata` | — | Metadata TSV (above) |
| `restricted` | — | `xqtl-restricted-<version>.tar.gz`, sent by MDACC by email ([Setup](SETUP.md)) |
| `site_refs` | — | Folder written by `--mode prepare_references` |
| `phenotypes_dir` | — | Ready phenotype BEDs; `--mode site` then skips phenotyping ([Export specification](EXPORT_SPEC.md#phenotype-inputs-phenotypes_dir)) |
| `image_sha256` | set by `run_xqtl.sh` | Checksum of the image, recorded in the export manifest |

### `genotypes`

| Key | Default | Meaning |
|---|---|---|
| `genotypes.raw` | — | Pre-imputation genotypes for `--mode imputation_prep`: bgzipped VCF/BCF, or a PLINK prefix (`.bed/.bim/.fam` or `.pgen/.pvar/.psam`) |
| `genotypes.raw_build` | `GRCh38` | Build of `raw`. `GRCh37` data are lifted over to GRCh38 |
| `genotypes.source` | `array` | `array`, `wes` or `wgs`. WES gives a warning: exome-only content imputes poorly outside exons |
| `genotypes.imputed_dir` | — | Folder with the decrypted TOPMed output (`chr*.dose.vcf.gz`) |
| `genotypes.imputed_pattern` | `chr{chr}.dose.vcf.gz` | File name pattern inside `imputed_dir`; `{chr}` is replaced by 1–22 |

### `imputation_prep`

Filters applied before upload to the TOPMed Imputation Server.

| Key | Default | Meaning |
|---|---|---|
| `imputation_prep.geno_missing_max` | 0.05 | Drop variants with a higher missing-call rate |
| `imputation_prep.maf_min` | 0.01 | Drop variants with a lower MAF |
| `imputation_prep.hwe_min` | 1e-6 | Drop variants with a lower Hardy–Weinberg exact-test p-value |
| `imputation_prep.sample_missing_max` | 0.05 | Report samples with higher autosomal missingness (not removed) |
| `imputation_prep.ref_match_min` | 0.95 | Stop if fewer SNPs than this match the GRCh38 reference (wrong build or strand) |

### `rnaseq`

| Key | Default | Meaning |
|---|---|---|
| `rnaseq.read_layout` | `paired` | `paired` or `single`; must agree with the samplesheet |
| `rnaseq.strandedness` | `reverse` | `reverse` (dUTP, Illumina TruSeq stranded), `forward` or `unstranded`; used by featureCounts, regtools and Salmon |
| `rnaseq.read_length` | 150 | Read length; the STAR index is built with `sjdbOverhang` = read length − 1. A mismatch with `site_refs` gives a warning |
| `rnaseq.library_type` | `total` | `total`, `polyA` or `3prime`. `3prime` disables the isoform, splicing, transcript-end and editing modalities |

### `qc` (post-imputation genotype QC)

| Key | Default | Meaning |
|---|---|---|
| `qc.rsq_min` | 0.4 | Minimum TOPMed imputation R² |
| `qc.maf_min` | 0.01 | Minimum MAF within an ancestry stratum (also the mapping MAF threshold) |
| `qc.mac_min` | 5 | Minimum minor-allele count within a stratum |
| `qc.hwe_p_min` | 1e-6 | Minimum Hardy–Weinberg p-value within a stratum |
| `qc.geno_missing_max` | 0.05 | Maximum variant missingness |
| `qc.sample_missing_max` | 0.05 | Maximum sample missingness; samples above it are removed |

### `ancestry`

| Key | Default | Meaning |
|---|---|---|
| `ancestry.groups` | `[AFR, AMR, EAS, EUR, SAS]` | Ancestry groups samples can be assigned to (1000 Genomes super-populations) |
| `ancestry.assign_prob_min` | 0.8 | Minimum classifier posterior probability for assignment; other samples are not assigned and not mapped |
| `ancestry.min_n` | 50 | Groups with fewer analysed samples are not mapped or exported (cannot be set below 50) |

### `mapping`

| Key | Default | Meaning |
|---|---|---|
| `mapping.modalities` | all nine | Modalities to map: `expression`, `isoform_expression`, `isoform_usage`, `splicing`, `intron_retention`, `alt_tss`, `alt_polya`, `rna_editing`, `stability` |
| `mapping.cis_window` | 1000000 | cis window in bp on each side of the feature TSS |
| `mapping.fdr` | 0.05 | q-value threshold for significant features (conditional analysis, signal sets) |
| `mapping.nperm` | 10000 | tensorQTL permutations |
| `mapping.hcp_grid` | automatic | Numbers of HCP factors to try; automatic is 0–60 by 10 for n < 350, otherwise 0–100 by 10 |
| `mapping.hcp_max_phenotypes` | 20000 | Phenotypes sampled (seeded) to fit HCP factors |
| `mapping.eval_chrom` | `chr1` | Chromosome on which eGenes are counted to choose the number of factors |
| `mapping.eval_min_phenotypes` | 300 | If `eval_chrom` has fewer phenotypes, the number of factors is chosen genome-wide |
| `mapping.independent_chunk_size` | 100 | Significant phenotypes per conditional-analysis task |
| `mapping.seed` | 20261029 | Seed for permutations and sampling |

### `ld` (file C)

| Key | Default | Meaning |
|---|---|---|
| `ld.maf_min` | 0.05 | Minimum MAF within the stratum for variants in the LD export |
| `ld.rsq_min` | 0.8 | Minimum imputation R² for variants in the LD export |
| `ld.window_mb` | 3 | LD window size, Mb |
| `ld.step_mb` | 1 | Step between windows, Mb (windows overlap by 2 Mb) |

### `gxc` (cell-type interaction)

| Key | Default | Meaning |
|---|---|---|
| `gxc.enabled` | true | Run G×C |
| `gxc.exclude_cell_types` | `[Maternal]` | Cell types never tested |
| `gxc.fdr` | 0.05 | BH q threshold for a pair's primary cell type |
| `gxc.min_n` | 20 | Minimum complete samples per fit |
| `gxc.maternal_threshold` | 0.1 | The sensitivity refit drops samples with a higher maternal proportion |

### `gxe` (exposure interaction)

| Key | Default | Meaning |
|---|---|---|
| `gxe.enabled` | true | Run GxE |
| `gxe.exposures` | GA on; `ppBMI`, `gdm`, `ogtt` off | One entry per exposure (below) |
| `gxe.sensitivity` | all off | SNP × covariate sensitivity refits of tier-1 hits (below) |
| `gxe.maf_interaction` | 0.05 | Minimum MAF for interaction tests |
| `gxe.min_n` | 50 | Minimum samples with the exposure in a stratum |
| `gxe.min_cases` | 20 | Binary exposures need at least this many cases and controls |
| `gxe.tier1_permutation_blocks` | `[100, 400, 500, 9000]` | Adaptive permutation blocks (cumulative 100, 500, 1000, 10000); a block runs only while needed |
| `gxe.tier2` | `all_cis_union` | `all_cis_union`: every cis variant of the tier-2 phenotypes; `pairs`: only the lead variants (smaller, a fallback) |
| `gxe.dosages` | false | Use imputed dosages instead of hard calls |
| `gxe.fdr` | 0.05 | q threshold for tier-1 hits |
| `gxe.seed` | 12345 | Base permutation seed |

Each entry of `gxe.exposures` (name: settings) accepts `gxe.exposures.<name>.enabled`, `.column` (metadata column, default the name), `.type` (`continuous` or `binary`), `.transform` (`none`, `log`, `inverse_normal`), `.min` and `.max` (values outside the range are set to missing; binary exposures must be 0/1). For example, to add pre-pregnancy BMI:

```yaml
gxe:
  exposures:
    ppBMI: {enabled: true}
```

Each entry of `gxe.sensitivity` (`maternal_age`, `parity`, `smoking`, `GA`, `ancestry_pcs`) accepts `gxe.sensitivity.<name>.enabled` and `.column`. `ancestry_pcs` uses the genotype PCs already in the mapping covariates. See [Interaction analyses](INTERACTIONS.md).

### `export`

| Key | Default | Meaning |
|---|---|---|
| `export.format` | `parquet` | `parquet` or `tsv.gz` |
| `export.min_mac_export` | 0 | Drop file A rows with a lower minor-allele count |
| `export.round_digits` | none | Round numeric columns to this many decimal places (smaller files) |
| `export.leak_scan` | true | Cannot be switched off ([Export specification](EXPORT_SPEC.md#leak-scan)) |

### `central` (MDACC only)

Used by `--mode central`; sites leave this block out. See [Central analysis](CENTRAL_ANALYSIS.md) for the table formats.

| Key | Default | Meaning |
|---|---|---|
| `central.studies` | required | TSV with `study_id`, `export_dir` |
| `central.gwas_catalog` | — | TSV with `trait_id`, `path`, `ancestry`, `n`, `trait_type`, `prop_cases`; enables colocalization |
| `central.annotations` | — | TSV with `name`, `bed`; annotation BEDs for PIP enrichment |
| `central.fdr` | 0.05 | Feature-level and interaction significance threshold |
| `central.palindromic_maf_max` | 0.4 | Palindromic SNPs with differing orientations are dropped above this MAF |
| `central.freq_diff_max` | 0.2 | Frequency-outlier flag threshold |
| `central.verify_checksums` | true | Check every file against `MANIFEST.tsv` |
| `central.gxc_min_contributions` | 2 | Minimum contributions for a G×C IVW |
| `central.pip_threshold` | 0.9 | High-confidence PIP threshold |
| `central.finemap.enabled` | true | Run SuShiE fine-mapping (also needed for colocalization and FUSION weights) |
| `central.finemap.L` | 10 | Maximum single effects per locus |
| `central.finemap.min_variants` | 50 | Minimum shared variants with LD per locus |
| `central.finemap.ld_ridge` | 0.001 | Minimum LD shrinkage |
| `central.finemap.purity` | 0.5 | Minimum credible-set purity |
| `central.finemap.max_loci` | all | Fine-map at most this many loci |
| `central.finemap.shard_size` | 20 | Loci per task |
| `central.coloc.enabled` | true | Run colocalization when `gwas_catalog` is set |
| `central.coloc.pp_h4` | 0.7 | PP.H4 threshold |
| `central.coloc.L` | 10 | SuSiE single effects per trait |
| `central.coloc.min_variants` | 50 | Minimum variants shared with the GWAS |
| `central.coloc.colocboost` | true | Also run ColocBoost |
| `central.coloc.M` | 500 | ColocBoost iterations |
| `central.coloc.shard_size` | 50 | Colocalization tasks per job |
| `central.twas.features` | `significant` | Features in the OTTERS files: `significant` or `all` |

## Command-line parameters

These are usually set by `run_xqtl.sh`, a site profile, or once on the command line.

| Parameter | Default | Meaning |
|---|---|---|
| `--sif` | set by `run_xqtl.sh` | Path to the image |
| `--image` | the GHCR image of this version | Container image when `--sif` is not set (Docker) |
| `--refs_dir` | `/opt/xqtl/refs` | Public references inside the image |
| `--genome_fasta`, `--liftover_chain`, `--annotation_gtf`, `--txrevise_dir`, `--editing_sites`, `--music_reference` | the file under `refs_dir` | Override one reference file (e.g. a single-chromosome genome) |
| `--star_sa_index_nbases` | 14 | STAR `--genomeSAindexNbases` (lower for small genomes) |
| `--normalize_min_samples` | 10 | Fewest samples a modality may be normalized on |
| `--expected_restricted_version` | 1 | Restricted-file version this pipeline version needs |
| `--validate_only` | false | Validate inputs and stop |
| `--publish_mode` | `copy` | How outputs are written to `outdir` (`copy`, `link`, `symlink`) |
| `--bind_test_path`, `--smoke_mem_hog_gb` | — | Smoke test only ([Setup](SETUP.md#5-check-that-your-cluster-works-smoke-test)) |

Scheduler and resource parameters (`--max_cpus`, `--max_memory`, `--max_time`, `--queue`, `--queue_rules`, `--gpu_queue`, `--gpu_cluster_options`, `--cluster_options`, `--use_gpu`, `--bind_paths`, `--lsf_per_job_mem`, `--queue_size`, `--submit_rate_limit`, `--job_array_size`) are described in [HPC configuration](HPC.md).

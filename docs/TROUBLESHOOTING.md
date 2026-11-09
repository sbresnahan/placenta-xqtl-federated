# Troubleshooting

Problems with the scheduler, memory, queues, bind paths or the container engine are covered in [HPC troubleshooting](HPC.md#hpc-troubleshooting). This page covers inputs, imputation, the site analysis, the export and the central analysis.

When something fails, first look at:

- `<outdir>/pipeline_info/validation_report.tsv` (input checks);
- the error block Nextflow prints, which names the failed task and its work directory;
- `.command.err` and `.command.log` in that work directory;
- `.nextflow.log` in the launch folder.

To ask MDACC for help, send the validation report, the Nextflow error block and the last lines of `.command.err`. These hold paths and messages, not data; check them for sample IDs before sending. Never send files from `<outdir>/site/`, `<outdir>/phenotypes/` or the work directory.

## Input validation messages

Every mode checks the inputs first and stops before submitting tasks if a check fails. Each line of the report has a check name, a status (`PASS`, `WARN`, `FAIL`) and a detail.

| Check | Meaning | Fix |
|---|---|---|
| `schema.<key>` | A setting has the wrong type, an unknown value or an unknown name (often a typo or wrong indentation in `cohort.yml`) | Compare with [Configuration](CONFIGURATION.md); `resolved_params.json` shows what was read |
| `required.<key>` | A setting this mode needs is missing, or its path does not exist | Add it to `cohort.yml`; check the path from a compute node |
| `required.restricted` | The restricted file is missing, unreadable, has a wrong checksum, or is a different version than this pipeline release expects | Use the file sent with this pipeline version; request it from MDACC if needed ([Setup](SETUP.md)) |
| `samplesheet.columns`, `.unique`, `.ids` | Column names, duplicate or empty `sample_id` | Header must be `sample_id,fastq_1,fastq_2` |
| `samplesheet.layout` | Rows disagree with `rnaseq.read_layout` | Empty `fastq_2` for single-end, filled for paired-end |
| `samplesheet.files` | FASTQ files not found | Relative paths are read from the samplesheet's folder |
| `metadata.columns`, `.unique` | Required columns missing or duplicate `rnaseq_id` | Required: `rnaseq_id`, `genotype_id`, `sex`, `GA` |
| `metadata.sex` / `metadata.sex_missing` | Sex codes not `M`/`F`; missing values (warning; imputed with the most common value) | Recode; fill in sex where known |
| `metadata.GA` | Gestational age missing, not numeric or outside 20–45 weeks | GA at delivery in weeks |
| `metadata.extra_columns` | Columns the pipeline ignores (warning) | Rename if an optional column was misspelled (`ppBMI`, `gdm`, `ogtt`, `maternal_age`, `parity`, `smoking`) |
| `link.samplesheet_metadata` | Samplesheet IDs not found as `rnaseq_id` | Every library needs a metadata row |
| `gxe.<exposure>` | An enabled exposure's column is missing, not numeric, not 0/1 (binary), or has few values | Fix the column or disable the exposure |
| `gxe.sensitivity.<name>` | An enabled sensitivity covariate's column is missing or not numeric | Fix the column or disable it |
| `mapping.modalities` | Isoform modalities requested without `expression` | Add `expression` |
| `site_refs` | `site_refs` missing or incomplete | Run `--mode prepare_references` |
| `site_refs.read_length` | Indices built for a different read length (warning) | Alignment still works; rebuild for best splice-junction detection |
| `site_refs.editing` | No editing sites in `site_refs` (warning) | `rna_editing` is skipped |
| `phenotypes_dir` | Missing folder or modality BEDs | See [Export specification](EXPORT_SPEC.md#phenotype-inputs-phenotypes_dir) |
| `imputation_prep.*` | Raw genotypes, GRCh38 FASTA or liftover chain not found; WES warning | Check `genotypes.raw` (a PLINK prefix has no extension) |
| `central.*` | Studies, GWAS catalog or annotation table problems | See [Central analysis](#central-analysis) below |

## Imputation preparation

| Message | Cause | Fix |
|---|---|---|
| `only N% of SNPs match the GRCh38 reference ... Wrong build?` | `genotypes.raw_build` does not match the data, or alleles are on the Illumina TOP/BOT strand | Set `raw_build: GRCh37` for hg19 data; convert TOP/BOT data to the forward strand first (`bcftools +fixref -m top`). No VCF is published |
| `duplicate sample IDs in the genotype input` | The same ID appears twice | Make IDs unique (and keep `genotype_id` in the metadata consistent) |
| Warning: chrX-inferred sex differs from the metadata | Possible sample swap or metadata error | Resolve before upload; swapped samples give wrong xQTL results |
| Warning: few SNPs in total | Low-density array or exome data | TOPMed imputation quality will be limited; discuss with MDACC |
| Samples above `imputation_prep.sample_missing_max` | Poor-quality samples (reported, not removed) | Consider removing them from the input |
| TOPMed QC report shows many strand flips or allele switches | Build or strand problem not caught here | Re-check the input build and strand |

## Site analysis

| Symptom | Cause | Fix |
|---|---|---|
| An ancestry group is missing from the export | Fewer than `ancestry.min_n` (50) analysed samples, or samples below `ancestry.assign_prob_min` | Expected; the group appears as `below_min_n` in `*_cohort_I_strata.tsv`. See `site/ancestry/cohort_ancestry.tsv` |
| Fewer samples than expected in a stratum | Genotype QC (sample missingness, genotype-PC outliers), expression outliers, unassigned ancestry, or IDs that do not match between metadata and VCFs | `site/<ANC>/<ANC>.samples.tsv`, `<ANC>.geno_outliers.tsv` and `*.outliers.tsv` |
| No imputed files found | `genotypes.imputed_pattern` does not match the file names | TOPMed names are `chr{chr}.dose.vcf.gz`; `{chr}` is 1–22 |
| Very few variants after QC | `qc.rsq_min` too strict, or the imputed files lack the TOPMed `INFO/R2` field | Use the TOPMed `.dose.vcf.gz` files unchanged |
| Modalities missing for a 3′ library | `rnaseq.library_type: 3prime` disables isoform, splicing, transcript-end and editing modalities | Expected |
| `gxc` or `gxe` is `stratum_below_min_n` or `no_cell_proportions` in `I_strata` | Too few samples, or no cell-type proportions | Expected; see [Interaction analyses](INTERACTIONS.md#run-status) |
| A GxE exposure has status `exposure_unavailable`, `insufficient_n` or `insufficient_cases` | Column missing, too few non-missing values, or too few cases or controls | Check the column and `gxe.min_n` / `gxe.min_cases` |
| Tasks killed for memory or time | Resource class too small for the cohort | [HPC configuration](HPC.md#retries-after-memory-or-time-kills) |
| The head process stopped | Login-node limits, wall time, maintenance | Run the same command with `-resume` from the same folder |

## Export

| Message | Cause | Fix |
|---|---|---|
| `LEAK SCAN FAILED ... the export was deleted` | A sample, RNA-seq or genotype identifier appears in an exported file name, column or value | Read `leak_report.tsv` in the export task's work directory (local only). Usually an ID that also occurs as a feature name or an exposure label; rename it, then `--mode export` |
| Export much larger than expected | File A covers every tested pair | `export.round_digits` and `export.min_mac_export` reduce the size; rebuild with `--mode export` |
| Need to change the format | | `export.format: tsv.gz`, then `--mode export` |

## Central analysis

Ingest writes `central/ingest/ingest_report.tsv` and stops on any `FAIL`.

| Check | Meaning | Fix |
|---|---|---|
| `manifest.<study>` FAIL: no `MANIFEST.tsv` | `export_dir` points at the wrong folder | Point it at `export/<cohort_id>/` |
| `manifest.<study>` FAIL: missing or checksum mismatch | The transfer was incomplete or a file was changed | Re-transfer the listed files; never edit received files |
| `versions.pipeline_version` / `versions.restricted_version` FAIL | Studies used different major pipeline versions or restricted files | Ask the study to rerun `--mode export` (or the site analysis) with the current release |
| `versions.image_sha256` WARN | Different image builds of compatible versions | Usually fine; record it |
| `strata.<study>` FAIL / INFO | No strata table; a stratum not mapped | FAIL: incomplete export. INFO: expected for small groups |
| `schema.<study>.<ANC>.<class>` FAIL | A file class is missing or lacks required columns | The export was made by an older or modified pipeline; ask for a rebuild |
| `central.studies` WARN: only one study | Meta-analysis reduces to within-study results | Expected while studies arrive |
| `central.gwas_catalog` / `central.annotations` FAIL | Missing columns, duplicate IDs or basenames, missing files, bad `trait_type` or `prop_cases` | Fix the table ([Central analysis](CENTRAL_ANALYSIS.md#running-it)) |

Other central results:

| Symptom | Cause | What to do |
|---|---|---|
| Many `drop_palindromic` or `mismatch` rows in `harmonization_summary.tsv` for one study | Strand or build problem at that site | Check the study's file D allele frequencies against the others; ask the site to check its imputation input |
| Many frequency outliers in one contribution | Wrong ancestry assignment, build, or allele coding | As above |
| Loci with status `too_few_variants_with_ld` | Little overlap between variants in file C across contributions | Lower `central.finemap.min_variants` with care; check `ld.*` settings used at the sites |
| Fine-mapping status `no_cs` | No credible set at the required purity | Expected for weak signals |
| Many secondary credible sets in strong loci | LD stored at low precision and averaged; very large effects | Review purity and PIPs; see [Caveats](CENTRAL_ANALYSIS.md#caveats) |
| Colocalization tasks `too_few_shared_variants` | GWAS lacks the locus variants, is on GRCh37, or uses other IDs | Lift the GWAS to GRCh38 and give `chrom`, `pos`, `ref`, `alt` |
| ColocBoost `no_cos` | No colocalization confidence set | Expected when the traits do not share a signal |
| FUSION finds no SNPs in the LD reference | The reference uses rsIDs | Rename the LD reference SNPs to `chr:pos:ref:alt` |

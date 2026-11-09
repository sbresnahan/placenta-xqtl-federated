# Walkthrough

This page follows one cohort from raw data to the files it shares, and then the MDACC side. Each step links to the page with the details.

| Step | Who | Mode | Result |
|---|---|---|---|
| 1. Setup and smoke test | site | `smoke` | `VERDICT GO` |
| 2. Input tables and `cohort.yml` | site | `--validate_only` | a clean validation report |
| 3. Genotypes for imputation | site | `imputation_prep` | per-chromosome GRCh38 VCFs and a TOPMed checklist |
| 4. TOPMed imputation | site | manual upload | decrypted `chr*.dose.vcf.gz` |
| 5. Site indices (once) | site | `prepare_references` | `site_refs/` |
| 6. Site analysis | site | `site` | `<outdir>/export/<cohort_id>/` |
| 7. Review the export | site | — | approval to share |
| 8. Transfer | site → MDACC | Globus or agreed method | received export |
| 9. Central analysis | MDACC | `central` | meta-analysis, fine-mapping, colocalization, TWAS inputs |

## 1. Setup and smoke test

Follow [Setup](SETUP.md): download and verify the image, run `bootstrap`, choose or write a site profile ([HPC configuration](HPC.md)), then run the smoke test until it reports `VERDICT GO`. Keep the restricted file sent by MDACC in a location only your team can read.

All later commands are run from a working folder on scratch, with the head process inside a batch job or `tmux` ([Running the head inside a batch job](HPC.md#running-the-head-inside-a-batch-job)):

```bash
cd /scratch/$USER/xqtl_COHORT1
XQTL=/project/xqtl/xqtl-1.0.0          # the bootstrap folder
```

## 2. Input tables and `cohort.yml`

Prepare the samplesheet, the metadata table and `cohort.yml` as described in [Configuration](CONFIGURATION.md). Then check them:

```bash
$XQTL/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode site --validate_only
```

Read `<outdir>/pipeline_info/validation_report.tsv`. Every `FAIL` must be fixed; read each `WARN` ([Troubleshooting](TROUBLESHOOTING.md#input-validation-messages)). Before imputation, `genotypes.imputed_dir` does not exist yet, so `required.genotypes.imputed_dir` fails; every other check should pass. The same command with `--mode imputation_prep` checks the raw genotypes.

## 3. Genotypes for imputation

```bash
$XQTL/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode imputation_prep
```

This reads `genotypes.raw` (VCF/BCF or PLINK, GRCh37 or GRCh38, array, WES or WGS) and writes to `<outdir>/imputation_prep/`:

- `chr1.vcf.gz` … `chr22.vcf.gz`: biallelic SNPs on GRCh38, REF/ALT aligned to the reference, filtered with the `imputation_prep` settings;
- `imputation_prep_report.tsv`: variant counts after each filter and the reference-match rate;
- `imputation_prep_samples.tsv`: per-sample missingness and the chrX sex check against the metadata;
- `liftover_summary.tsv` (GRCh37 input only);
- `TOPMED_CHECKLIST.md`: the upload settings for this cohort.

Check the reports before uploading. A reference-match rate below `imputation_prep.ref_match_min` stops the run with "Wrong build?" and publishes no VCF; see [Troubleshooting](TROUBLESHOOTING.md#imputation-preparation). Samples whose genotype sex disagrees with the metadata may be swapped; resolve them before imputation.

## 4. TOPMed imputation

The upload is manual and follows your institution's approval for sending genotypes to the TOPMed Imputation Server. Follow `TOPMED_CHECKLIST.md`:

- reference panel TOPMed R3 on GRCh38; array build GRCh38; rsq filter off; phasing with Eagle v2.4; mode Quality Control & Imputation;
- the one-time password for the results arrives by email and cannot be re-sent;
- results are deleted from the server after 7 days.

Unzip the results into one folder and set `genotypes.imputed_dir` (and `imputed_pattern` if your file names differ from `chr{chr}.dose.vcf.gz`).

## 5. Site indices (once per site)

```bash
$XQTL/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode prepare_references
```

This builds the STAR, Salmon and IRFinder-S indices into `site_refs` ([Phenotyping](PHENOTYPING.md#before-the-first-run-site-indices)). It is reused by later runs and pipeline versions unless the changelog says otherwise. Skip this step if you supply ready phenotypes with `phenotypes_dir`.

## 6. Site analysis

```bash
$XQTL/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode site
```

The run goes through:

1. phenotyping from FASTQs: nine modalities and cell-type proportions ([Phenotyping](PHENOTYPING.md));
2. post-imputation genotype QC and ancestry assignment;
3. per ancestry group with at least `ancestry.min_n` samples: normalization, hidden factors, cis-xQTL mapping, conditional analysis and LD;
4. G×C and GxE interaction tests ([Interaction analyses](INTERACTIONS.md));
5. the export with its leak scan ([Export specification](EXPORT_SPEC.md)).

If the run stops (a failed task, a killed head process, a maintenance window), fix the cause and run the same command with `-resume` from the same folder; finished tasks are not repeated.

Outputs:

| Folder | Content | Shared? |
|---|---|---|
| `<outdir>/phenotypes/` | Phenotype BEDs, QC metrics, cell-type proportions | No |
| `<outdir>/site/` | Genotype QC, ancestry, per-stratum mapping, interaction results and diagnostics | No |
| `<outdir>/pipeline_info/` | Validation report, settings, Nextflow trace and reports, `leak_report.tsv` | No |
| `<outdir>/export/<cohort_id>/` | Files A–J | Yes, after review |

## 7. Review the export

Before anything is shared, the site team reviews `export/<cohort_id>/`:

- `<cohort_id>_cohort_I_strata.tsv`: which ancestry groups were mapped and exported, and which were below `ancestry.min_n`;
- `DATA_DICTIONARY.md`: every column of every file;
- `MANIFEST.tsv`: every file with its size, row count and sha256, the pipeline and restricted-file versions and the image checksum;
- the leak scan passed: otherwise there is no export folder, and `pipeline_info/leak_report.tsv` has only its header line.

The [Export specification](EXPORT_SPEC.md#what-is-never-exported) lists what is never exported. If the cohort wants different rounding or a minimum minor-allele count in file A, change `export.round_digits` or `export.min_mac_export` and rebuild the export without recomputing anything:

```bash
$XQTL/run_xqtl.sh -profile <site profile> -params-file cohort.yml --mode export
```

## 8. Transfer

Send only the `export/<cohort_id>/` folder, by the method agreed with MDACC (for example Globus), keeping its structure. Expect roughly 100–150 GB per ancestry group, mostly files A and C. Tell MDACC the release tag and the sha256 of `MANIFEST.tsv`:

```bash
sha256sum export/COHORT1/MANIFEST.tsv
```

## 9. Central analysis (MDACC)

MDACC stores each received export in its own dated folder, verifies `MANIFEST.tsv`, adds the study to `studies.tsv` and runs `--mode central` ([Central analysis](CENTRAL_ANALYSIS.md)). Ingest checks every checksum again and stops if a file is missing or changed, or if the exports were made with incompatible pipeline versions. MDACC shares combined results with the cohorts; per-cohort files are not redistributed ([Collaboration brief](COLLABORATION_BRIEF.md#7-data-sharing-and-safeguards)).

# Central analysis (MDACC)

`--mode central` is run by MDACC. It combines the exports (files A–J, see [Export specification](EXPORT_SPEC.md)) received from the participating cohorts. It never needs individual-level data, and it runs offline from the same image as the site analysis.

The run has eight stages:

| Stage | What it does | Files used | Output folder |
|---|---|---|---|
| Ingest | Checks checksums, versions and file schemas; harmonizes alleles across contributions | J, I, A, B, D, E, G, H | `central/ingest/` |
| Nominal meta-analysis | Fixed- and random-effects inverse-variance weighting (IVW) of every variant–feature pair | A | `central/meta/` |
| Feature level | Combines per-cohort permutation p-values (ACAT); q-values; meta lead variant | B | `central/features.tsv.gz` |
| Fine-mapping | Multi-ancestry SuShiE in summary-statistics mode | A (meta), C | `central/finemap/` |
| Colocalization | SuSiE-coloc and ColocBoost against GWAS summary statistics | A (meta), C, GWAS | `central/coloc/` |
| Interactions | G×C and GxE synthesis | G, H | `central/interactions/` |
| TWAS inputs | OTTERS / mvOTTERS input files and FUSION weights | A (meta), E, fine-mapping | `central/twas/` |
| Summary | Counts per stage and enrichment of high-PIP variants in annotations | all of the above | `central/summary/` |

A **contribution** is one study × ancestry stratum: an `<ANC>/` folder whose status is `mapped` in the study's `*_cohort_I_strata.tsv`. Contributions are numbered from 0 in `ingest/contributions.tsv`. Every meta-analysis table has a `contrib_mask` column, a bit mask over those numbers: bit *i* is set when contribution *i* contributed. Meta-analyses are run within each ancestry (`scope` = the ancestry code) and across all contributions (`scope` = `ALL`).

## Running it

```bash
./xqtl/run_xqtl.sh -profile seadragon -params-file central.yml --mode central
```

`central.yml` needs only a cohort label, an output folder and the studies table. GWAS colocalization and annotation enrichment are optional:

```yaml
cohort_id: MDACC_central
outdir: /rsrch5/.../xqtl_central_2026_11
central:
  studies: studies.tsv            # required
  gwas_catalog: gwas_catalog.tsv  # optional: enables colocalization
  annotations: annotations.tsv    # optional: annotation BEDs for PIP enrichment
  fdr: 0.05
  finemap:
    L: 10
  coloc:
    pp_h4: 0.7
```

Every setting and its default is listed in [Configuration](CONFIGURATION.md#central-mdacc-only). Relative paths in the three tables are read relative to the table's own folder.

**`studies.tsv`** lists one received export per study. `export_dir` is the folder that holds `MANIFEST.tsv` (`export/<cohort_id>/` as written by the site):

```text
study_id	export_dir
COHORT1	/rsrch5/.../received/COHORT1_20261201/export/COHORT1
COHORT2	/rsrch5/.../received/COHORT2_20261215/export/COHORT2
```

**`gwas_catalog.tsv`** lists one GWAS per row. Each file needs GRCh38 coordinates, with either `variant_id` (`chr:pos:ref:alt`) or the columns `chrom`, `pos`, `ref`, `alt`, plus `beta` and `se` for the ALT allele; `eaf` and a per-variant `n` are optional. `ancestry` must be one of the ancestry codes used by the cohorts (`AFR`, `AMR`, `EAS`, `EUR`, `SAS`). For `trait_type: cc`, give `prop_cases`, the proportion of cases:

```text
trait_id	path	ancestry	n	trait_type	prop_cases
BW_EUR	gwas/birthweight_EUR.tsv.gz	EUR	298142	quantitative
GDM_EUR	gwas/gdm_EUR.tsv.gz	EUR	180000	cc	0.08
```

**`annotations.tsv`** lists BED files (GRCh38, 0-based) to test for enrichment of fine-mapped variants:

```text
name	bed
placenta_enhancers	annot/placenta_enhancers.bed.gz
```

The input validator checks all three tables before any task starts: required columns, duplicate IDs, folders and files that exist, `MANIFEST.tsv` and `*_cohort_I_strata.tsv` in each export folder, GWAS trait types and `prop_cases`. It warns when only one study is listed.

### Receiving exports

Keep every export exactly as the site sent it: the folder structure and file names are listed in `MANIFEST.tsv` with sha256 checksums, and ingest verifies every file. Put each transfer in its own dated folder, then point `export_dir` at the `export/<cohort_id>/` folder inside it. A study that sends a corrected export gets a new folder; do not overwrite the old one.

## Ingest and harmonization

Ingest stops the run (status `FAIL` in `ingest_report.tsv`) when:

- a folder has no `MANIFEST.tsv`, or any listed file is missing or has the wrong sha256 (`central.verify_checksums: false` skips the hashing, not the presence check);
- the exports do not share one export format: every study must have the same **major** pipeline version and the same restricted-file version;
- a study has no `*_cohort_I_strata.tsv`, or no stratum with status `mapped`;
- a required file class (A, B permutation, D, E) is missing, or a file lacks a required column.

Strata that were not mapped (for example `below_min_n`) are reported as `INFO` and skipped. Different image checksums across studies give a `WARN` only.

**Allele harmonization.** All exports use GRCh38 forward-strand IDs `chr:pos:ref:alt` with effects per ALT allele, so in a correct export no flips are needed; ingest still checks. For each site (chromosome, position and allele pair, strand-collapsed), the orientation used by most contributions becomes the canonical one. Other contributions are then:

- `swap` (REF and ALT exchanged): effect sign reversed, EAF = 1 − EAF;
- `flip` or `flip_swap` (opposite strand): complemented, and reversed when swapped;
- `mismatch` (alleles that cannot be reconciled): dropped.

Palindromic SNPs (A/T, C/G) cannot be told apart from a strand flip. A palindromic site that has the same orientation in every contribution is kept, because the IDs are forward-strand. When orientations differ, the site is dropped from every contribution if its MAF exceeds `central.palindromic_maf_max` (default 0.4) in any of them; below that it is read as a swap.

**Frequency outliers.** A contribution's variant is flagged when its harmonized EAF differs from the N-weighted mean EAF of the other contributions of the same ancestry by more than `central.freq_diff_max` (default 0.2). Flags are counted in the meta tables (`n_freq_outlier`); flagged variants are not removed. Many flags in one contribution suggest a build, strand or sample problem at that site.

Outputs: `ingest_report.tsv` (check, status, detail), `contributions.tsv`, `versions.tsv`, `harmonization_summary.tsv` (actions per contribution) and `harmonization/harm_<chrom>.parquet` (one row per contribution × variant).

## Nominal meta-analysis

For every modality and chromosome in file A, each variant–feature pair is meta-analysed by inverse-variance weighting:

- fixed effects: β = Σwᵢβᵢ / Σwᵢ with wᵢ = 1/SEᵢ², SE = (Σwᵢ)^−½;
- random effects: DerSimonian–Laird τ², weights 1/(SEᵢ² + τ²);
- heterogeneity: Cochran's Q with its p-value, and I² = max(0, (Q − df)/Q).

The meta EAF and N are N-weighted means and sums. Unit tests check these against a hand calculation. Splicing features in site-specific novel clusters (`nclu_`, `novel` = true in file E) are excluded, because their IDs differ between sites ([Phenotyping](PHENOTYPING.md#splicing-frozen-and-novel-clusters)).

Output: `meta/<modality>_<chrom>.parquet` with `scope, modality, chrom, phenotype_id, variant_id, tss_distance, n, eaf, k, beta_fe, se_fe, z_fe, p_fe, beta_re, se_re, p_re, tau2, q_het, p_het, i2, n_freq_outlier, contrib_mask`.

## Feature-level significance

Each contribution's beta-approximated permutation p-value (`pval_beta`, file B) measures the evidence that a feature has any cis-xQTL in that contribution. They are combined per arm, scope and feature with the aggregated Cauchy association test (ACAT), which is valid under arbitrary dependence [1, 2]. Storey q-values are computed per arm and scope, with Benjamini–Hochberg as the fallback when the q-value estimate fails (`qvalue_method`). A feature is `significant` when `q_acat` ≤ `central.fdr`.

For ungrouped arms, the meta lead variant is the smallest fixed-effect p-value among the feature's nominal meta rows (`lead_variant`, `lead_beta_fe`, `lead_p_fe`, `lead_p_re`, `lead_i2`, `lead_k`). `lead_concordant` is true when every contribution reported the same lead variant; `n_sig_contrib` counts contributions with site q ≤ `central.fdr`.

Output: `features.tsv.gz`.

## Fine-mapping

Every feature that is significant across all contributions (`scope` `ALL`, ungrouped arms) is fine-mapped, smallest p-value first (`central.finemap.max_loci` caps the number). For each ancestry:

- z-scores are the within-ancestry fixed-effect meta z, adjusted for the variance explained (z·√((n − 1)/(z² + n − 2)), the susieR convention);
- LD is the N-weighted mean of the contributions' signed in-sample correlation matrices (file C), sign-aligned to the canonical alleles. Because file C is stored as float16, the matrix is shrunk towards the identity by λ = max(`ld_ridge`, `ld_ridge` − 1.1 × its smallest eigenvalue).

The variant set is shared across ancestries and must have LD in every contribution of an ancestry. Ancestries are added largest first while at least `central.finemap.min_variants` variants remain; loci with fewer get status `too_few_variants_with_ld`. Only variants in file C (MAF ≥ `ld.maf_min`, imputation R² ≥ `ld.rsq_min` at the sites) can be fine-mapped.

SuShiE [3] then fits up to `central.finemap.L` single effects jointly across ancestries, allowing effect sizes to differ between ancestries, and reports 95% credible sets with purity ≥ `central.finemap.purity`.

Outputs in `finemap/`:

| File | Content |
|---|---|
| `loci.tsv` | One row per locus: variants, ancestries kept, region, preparation status |
| `status.tsv` | One row per locus: `ok`, `no_cs` or `error`, number of credible sets, run time |
| `pips.tsv` | Variant PIPs, PIP within credible sets, credible-set index |
| `cs.tsv` | Credible-set members with alpha and cumulative alpha |
| `corr.tsv` | Estimated cross-ancestry effect correlation per single effect |
| `weights.tsv` | Posterior effect per variant and ancestry (prediction weights, standardized scale) |

## Colocalization

When `central.gwas_catalog` is set, every fine-mapped locus is tested against every GWAS whose ancestry is among the locus's ancestries. GWAS variants are matched to the locus by ID, or by the allele-swapped ID with the sign reversed; strand flips are not attempted, so the GWAS must be on the GRCh38 forward strand. Tasks with fewer than `central.coloc.min_variants` shared variants are skipped.

Both traits use the same LD: the ancestry-matched N-weighted xQTL LD described above. No external LD panel is used. This is appropriate when the GWAS ancestry matches the cohorts, but GWAS signals in poorly matched LD can give spurious extra credible sets, particularly near the edge of an LD window.

- **SuSiE-coloc** [4]. `susieR::runsusie` is fitted to each trait with up to `central.coloc.L` effects (xQTL with sdY = 1; GWAS with its `n`, `trait_type` and `prop_cases`), then `coloc.susie` compares every pair of credible sets. When either trait has no credible set, `coloc.abf` (single causal variant) is used instead (`method` = `abf`). `coloc.tsv` keeps the pair with the highest PP.H4 per task; `coloc_pairs.tsv` keeps all pairs. A task colocalizes (`coloc` = TRUE) when PP.H4 ≥ `central.coloc.pp_h4`.
- **ColocBoost** [5, 6] (`central.coloc.colocboost`, default on). Multi-trait gradient boosting on the same z-scores and LD with `central.coloc.M` iterations; `colocboost.tsv` lists colocalization confidence sets with their purity, top variant and variant-level colocalization probabilities, or `no_cos`.

Treat SuSiE-coloc and ColocBoost as two readouts of the same data, not independent replication.

## Interactions

**G×C** (file G). For every lead pair × cell type, the interaction coefficients are harmonized and meta-analysed by IVW within each ancestry and across all contributions when at least `central.gxc_min_contributions` contributed (`method` = `ivw`); otherwise the single p-value is kept (`method` = `single`). BH q-values are computed per scope and modality. The central **primary cell type** of a pair is its cell type with the smallest q, if q ≤ `central.fdr`; otherwise `unassigned`. The maternal-contamination refits are meta-analysed the same way (`maternal_*` columns). Outputs: `interactions/gxc_meta.tsv.gz` and `gxc_annotation.tsv.gz`.

**GxE tier 1** (atlas script 54). Per modality × exposure × phenotype, the contributions' adaptive-permutation p-values are combined by ACAT (`p_feature_acat`), with BH q-values over all phenotypes of the modality × exposure (`q_feature_acat`). The **scope label** describes agreement between contributions:

| Scope | Meaning |
|---|---|
| `single_study` | one contribution tested the phenotype |
| `shared_concordant_lead` | several contributions, all with the same lead interaction variant; that variant's interaction effect is meta-analysed by IVW (`beta_meta`, `p_meta`, `i2`) |
| `shared_discordant_lead` | several contributions with different lead variants; only the ACAT p-value is given |

`significant_any` and `significant_all` show whether the site q ≤ `central.fdr` in any or all contributions. Output: `interactions/gxe_tier1_meta.tsv.gz`. The atlas computed this BH over features shared by several contributions only; here it covers all features, including `single_study` ones.

**GxE tier 2.** Sites estimate interaction effects per SD of the exposure in their own stratum. For continuous, untransformed exposures these are rescaled to a common unit, the pooled SD over all contributions (`gxe_exposure_scales.tsv`), using the raw-unit SDs in file H. Variant-level effects are then meta-analysed by IVW per modality × exposure. The top variant per phenotype is flagged and BH-adjusted across phenotypes (`q_top_bh`, descriptive only, as at the sites). Outputs: `interactions/tier2/gxe_tier2_meta_<modality>_<exposure>.parquet` and `gxe_tier2_top.tsv.gz`.

See [Interaction analyses](INTERACTIONS.md) for the site-side tests.

## TWAS inputs

The pipeline writes TWAS inputs; it does not run the TWAS. MDACC runs mvOTTERS separately.

- **OTTERS / mvOTTERS stage I** [7, 8]: per ancestry, `twas/otters/<ANC>/exp_anno.txt` (`CHROM GeneStart GeneEnd TargetID`, with the feature TSS from file E) and `eQTLSumStats.txt` (`CHROM POS A1 A2 Zscore TargetID N`; A1 is the effect allele, ALT). Z-scores are the within-ancestry fixed-effect meta z. `TargetID` is `<modality>:<phenotype_id>`. `central.twas.features` chooses the significant features (default) or all features.
- **FUSION weights** [9] from SuShiE: `twas/fusion/<ANC>/WEIGHTS/<modality>_<phenotype>.wgt.RDat` and `<ANC>.pos`. Each `.wgt.RDat` holds `wgt.matrix` (one column, `sushie`), `snps` (PLINK `.bim` layout, SNP IDs `chr:pos:ref:alt`, A1 = ALT), `cv.performance` (NA; no cross-validation), `hsq` (NA), `hsq.pv` (NA) and `N.tot`. Run FUSION with `--force_model sushie`, and use an LD reference whose SNP IDs are `chr:pos:ref:alt`. Weights are written for every locus SuShiE fitted, including loci without a credible set (`no_cs` in `finemap/status.tsv`); their posterior effects are small, and you can drop them using `status.tsv`.

## Summary and enrichment

`summary/central_summary.tsv` (section, item, value) counts contributions, individuals per ancestry, tested and significant features per arm and scope, fine-mapped loci, credible sets and variants with PIP ≥ `central.pip_threshold`, colocalized tasks, assigned G×C cell types and GxE scope labels. `ingest_report.tsv` and `harmonization_summary.tsv` are copied next to it.

`summary/pip_enrichment.tsv` tests whether high-PIP variants (PIP ≥ `central.pip_threshold`; best PIP per variant across loci) fall in each annotation class more often than the other fine-mapped variants: Fisher's exact test with BH q-values, and a PIP-weighted enrichment (share of total PIP in the class divided by the share of variants in the class). Classes are TSS-distance bins (≤ 1 kb, 1–10 kb, 10–100 kb, > 100 kb) and every BED in `central.annotations`. The background is the fine-mapped variant set, so these are enrichments within xQTL loci, not genome-wide.

## Caveats

- **LD comes from the cohorts.** In-sample LD is stored as float16 and averaged across contributions; with small strata and large effects SuShiE can report extra secondary credible sets. Review secondary sets with low purity or in loci with very strong primary signals before interpreting them.
- **Variant coverage.** Fine-mapping, colocalization and FUSION weights are limited to variants in file C at every site of an ancestry.
- **Unequal coverage.** A variant–feature pair can be meta-analysed with fewer contributions than its feature (`k` versus the number of contributions in `contrib_mask`), because sites filter variants independently.
- **Exposure scales.** Tier-2 rescaling assumes linear interaction effects in raw exposure units; log-transformed and binary exposures are not rescaled.
- **One ancestry per contribution.** Admixed individuals are assigned by the site's ancestry classifier; the central analysis cannot revisit that assignment.

## What the shared files cannot support

The files A–J are designed for the analyses above. Some analyses need data that do not leave the sites. They need a new pipeline release run by the sites, or, for pooled individual-level data, a separate data-transfer agreement (see [Collaboration brief](COLLABORATION_BRIEF.md#8-planned-extensions-beyond-the-first-release)).

| Analysis | Possible from A–J? | What is missing | Route |
|---|---|---|---|
| Within- and cross-ancestry meta-analysis, feature-level significance | Yes | — | Central |
| Multi-ancestry fine-mapping, GWAS colocalization | Yes, for variants in file C | LD for rare or poorly imputed variants | Central; an external LD panel for other variants |
| TWAS input files and weights | Yes | Out-of-sample prediction accuracy (`cv.performance`) | Central; sites could return aggregate prediction metrics in a later release |
| Heritability partitioning (S-LDSC, graphREML, MESC) | Yes, with public GWAS and LD resources | — | Central |
| G×C and GxE beyond the shared test sets (tier 2 genome-wide, new exposures) | No | Interaction statistics for untested pairs and exposures | New release run at the sites |
| Sex- or exposure-stratified mapping | No | Stratified association tables | New release run at the sites, subject to minimum cell sizes |
| Different normalization, covariates or new phenotype definitions | No | Per-sample phenotypes and covariates | New release run at the sites |
| Pooled (mega-) analysis, joint conditional analysis on individual genotypes | No | Individual genotypes and phenotypes | Data-transfer agreement |
| Ancestry reassignment, local-ancestry or admixture-aware mapping | No | Individual genotypes | Site analysis, or data-transfer agreement |
| Sample-level QC review (sample swaps, outliers, relatedness) | No | Individual data | Site only |
| varQTLs, maternal eQTLs, eQTMs and other site-directed analyses | No | Individual data | Site only (no files to share) |

## Test fixture

`tests/fixtures/central/make_gwas.py` writes a synthetic GWAS from one export's LD window: one causal variant shared with a planted cis effect (colocalizing) and one nearby independent causal variant (a distinct signal, H3). `tests/scripts/test_central.sh` runs `--mode site` for two independent test cohorts, then `--mode central`. Cohort 1 is the site fixture (`tests/fixtures/site`, parquet export). Cohort 2 (`tests/fixtures/site2`, built by `tests/fixtures/make_central_fixture.sh`) has 120 other 1000 Genomes samples (60 EUR, 60 EAS), the same planted effects with other random seeds, and is exported as `tsv.gz`. Splitting one fixture into two cohorts was not possible: each ancestry stratum needs at least `ancestry.min_n` = 50 samples, and that floor is a privacy safeguard. The test adds a BED of the planted variants (±1 kb) as an enrichment annotation and gives study 2 as a path relative to `studies.tsv`. `tests/scripts/check_central.py` checks ingest, the planted signals at meta, feature, fine-mapping, colocalization, G×C and GxE level, the tier-1 scope labels, the TWAS input files and the summary. `tests/scripts/check_central_ivw.R` recomputes the fixed- and random-effects IVW (Q, I², τ²) in R from the raw file A rows of both exports and loads every FUSION weight file. The test runs in CI once the development image contains the SuShiE environment.

## References

1. Liu Y, Xie J. Cauchy Combination Test: A Powerful Test With Analytic p -Value Calculation Under Arbitrary Dependency Structures. Journal of the American Statistical Association. 2019. doi:10.1080/01621459.2018.1554485.
2. Yao-Wu Liu, Si-Xing Chen, Zi-Lin Li, et al. ACAT: A Fast and Powerful P-value Combination Method for Rare-variant Analysis in Sequencing Studies. bioRxiv. 2018. doi:10.1101/482240.
3. mancusolab/sushie. [cited 2026 Oct 28]. Available from: https://github.com/mancusolab/sushie
4. C. Wallace. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genetics. 2021. doi:10.1101/2021.02.23.432421.
5. Xue-Wei Cao, Hao-Chen Sun, R. Feng, et al. Integrative multi-omics QTL colocalization maps regulatory architecture in aging human brain. medRxiv. 2025. doi:10.1101/2025.04.17.25326042.
6. Cao X, Sun H, Feng R, Nachun D, Dey K, Wang G. colocboost: Multi-Context Colocalization Analysis for QTL and GWAS Studies. CRAN: Contributed Packages. 2025. doi:10.32614/cran.package.colocboost.
7. Qile Dai, Geyu Zhou, Hong-Yu Zhao, et al. OTTERS: a powerful TWAS framework leveraging summary-level reference data. Nature Communications. 2023. doi:10.1038/s41467-023-36862-w.
8. daiqile96/OTTERS. [cited 2026 Nov 8]. Available from: https://github.com/daiqile96/OTTERS
9. Alexander Gusev, Arthur Ko, Huwenbo Shi, et al. Integrative approaches for large-scale transcriptome-wide association studies. Nature genetics. 2015. doi:10.1101/024083.

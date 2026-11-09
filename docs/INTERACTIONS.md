# Interaction analyses (G×C and GxE)

`--mode site` runs two interaction analyses in every mapped ancestry stratum after cis mapping:

- **G×C**: genotype × cell-type proportion at lead xQTL pairs. Exported as file G.
- **GxE**: genotype × maternal or pregnancy exposure (gestational age by default). Exported as file H.

Both reuse the stratum's normalized phenotypes, genotypes and final mapping covariates; nothing is re-normalized. They port the atlas code with the statistics unchanged: G×C from `06_colocalization_twas/49_aggregate_coloc_twas.py` (`_fit_gxc`), GxE from `07_gxe_mapping` (scripts 50, 52, 54 and 55 and `gxe_core.py`). The sections below list where the federated pipeline differs.

Settings are in the `gxc` and `gxe` blocks of `cohort.yml`. Both analyses are on by default. Turn one off with `enabled: false`. The run state of each stratum is in the `gxc` and `gxe` columns of `<cohort>_cohort_I_strata.tsv`:

- `tested`;
- `disabled`;
- `no_cell_proportions` (G×C only);
- `stratum_below_min_n`.

## Test sets

The test sets are built per stratum and modality from two sources:

- the site's own leads;
- the MDACC signal list in the restricted file (`mdacc_signals.tsv.gz`: atlas lead pairs per discovery ancestry).

**Site leads.** These are the top variants of phenotypes with q ≤ 0.05 in the stratum's feature-level permutation table: the `expression` table for expression and the `<modality>_ungrouped` table for the other modalities. The atlas used the same tables.

**G×C pairs (leads only).** The union of:

- site leads;
- MDACC pairs with `signal_rank = 1` from any discovery ancestry.

Every pair is tested in every stratum, whatever ancestry it was discovered in. `source` is `site`, `mdacc` or `both`; `both` means the same phenotype and variant. `mdacc_ancestry` records the discovery ancestry. A pair that cannot be tested is kept with `untested_reason`:

- `phenotype_absent`: the phenotype is not in the stratum after filtering;
- `variant_absent`: the variant did not pass the stratum's genotype QC.

**GxE tier-2 phenotypes.** Site-lead phenotypes plus every phenotype in the MDACC list (any rank). All cis variants of each one are tested. The site-lead and MDACC-lead variants are flagged (`site_lead_variant`, `mdacc_lead_variant`).

**Intron retention.** The atlas used MAJIQ segments; sites use IRFinder-S introns ([Phenotyping](PHENOTYPING.md#differences-from-the-atlas)). An MDACC segment is matched to the shortest site intron of the same gene, chromosome and strand that contains it. The original ID is kept in `mdacc_phenotype_id`.

## G×C: cell-type interaction

For each pair and target cell type, an OLS model with a t-test:

```
Y ~ 1 + g0 + c0 + g0:c0 + Z
```

| Term | Definition |
|---|---|
| Y | the stratum's normalized phenotype |
| g0 | ALT dosage, centred over the analysed samples |
| c0 | arcsinh(cell proportion), centred over the stratum, then re-centred over the analysed samples. Proportions are averaged over technical replicates |
| Z | the modality's final mapping covariates, without the target's own cell-type covariate. Covariates that are constant in the analysed samples are dropped |

**Targets.** Every cell type in `cell_proportions.tsv` except `gxc.exclude_cell_types` (default `Maternal`, the decidual fraction).

**Maternal sensitivity.** The same fit is repeated after removing samples whose maternal proportion is above `gxc.maternal_threshold` (default 0.10). These results are in the `maternal_sensitivity_*` columns.

**Fit status.** One of:

- `ok`;
- `no_common_samples`;
- `insufficient_n` (fewer than `gxc.min_n` complete samples, default 20);
- `monomorphic`;
- `no_celltype_variation`;
- `rank_deficient`;
- `insufficient_df`;
- `invalid_se`.

**Multiple testing.** BH q (`q_GxC`) over the `ok` tests, within stratum × modality.

**Annotation.** One row per pair. `best_cell_type` is the cell type with the smallest (q, p). `primary_cell_type` is the best cell type if its q ≤ `gxc.fdr` (default 0.05), otherwise `unassigned`. Annotation is per pair, not per phenotype: a phenotype with a site pair and a different MDACC pair has two rows.

Effects are also reported on an interpretable scale:

- `beta_GxC_iqr`: the interaction per interquartile range of the transformed proportion;
- `beta_G_at_cell_q25/q50/q75`: the genetic effect at the proportion quartiles;
- `cell_prop_raw_q25/q50/q75`: the quartiles of the raw proportion.

G×C needs cell-type proportions. Without them, G×C is skipped and I_strata records `no_cell_proportions`. Proportions come from MuSiC when phenotyping runs, or from `phenotypes_dir/cell_proportions.tsv`.

## GxE: exposure interaction

### Exposures

Each enabled exposure in `gxe.exposures` is read from its metadata column and processed in five steps:

1. averaged over technical replicates;
2. filtered to `min`–`max` on the raw scale;
3. transformed (`none`, `log` or `inverse_normal`);
4. continuous exposures are z-scored across all samples of the mapped strata;
5. binary exposures (0/1) are left unscaled.

Pooling the scale across strata, as the atlas did, puts every stratum's b_GxE on one scale for the central meta-analysis.

An exposure is used in a stratum only if both conditions hold:

- n ≥ `gxe.min_n` (default 50);
- for binary exposures, at least `gxe.min_cases` (default 20) cases and controls.

Otherwise its status is `insufficient_n` or `insufficient_cases`, and every component is `SKIPPED`. The exposure summary (H_exposures) gives n, the raw-unit mean and SD, the pooled scale, and counts and prevalence for binary exposures. The per-sample table stays at the site.

Defaults:

- `GA` is enabled;
- `ppBMI`, `gdm` and `ogtt` are defined but off.

You can add other exposures with the same fields. The validator fails if an enabled exposure has no metadata column.

### Model

```
Y = b_g G + b_e E + b_gi (G × E) + Z γ + ε
```

| Term | Definition |
|---|---|
| Y | the normalized phenotype, z-scored within the stratum |
| G | hard calls (atlas), or dosages with `gxe.dosages: true` |
| E | the exposure |
| Z | the final mapping covariates. A covariate named like the exposure is dropped |

The model is fitted per cis variant (±1 Mb). The algebra is tensorQTL's `calculate_interaction_nominal`:

- covariates are residualized out with a QR projection;
- dof = n − k − 4.

Variant filters:

- MAF ≥ 0.01 in the analysed samples;
- MAF ≥ `gxe.maf_interaction` (default 0.05) in both halves of the sorted exposure, so that the interaction is estimable.

The scans are sharded by chromosome. The scans use the process label `process_gpu_capable`, so with `--use_gpu` they run on GPU nodes ([HPC](HPC.md#gpus)).

### Tier 1: permutation scan of every phenotype

Every phenotype is tested; its statistic is the top cis variant for the interaction term. Significance comes from adaptive phenotype permutations:

- blocks are given by `gxe.tier1_permutation_blocks`, default 100, 400, 500 and 9,000 (cumulative 100, 500, 1,000 and 10,000);
- a further block runs only while the running empirical p ≤ 0.10;
- the null distribution is Beta-approximated (`pval_beta`), as in tensorQTL `map_cis`.

Phenotype *i* uses the seed `gxe.seed + i`, so results do not depend on sharding. Storey q-values (`qval`) are computed within stratum × modality × exposure, with a BH fallback recorded in `qvalue_method`. A non-finite `pval_beta` gives q = 1 (atlas).

### Tier 2: all cis variants of the test-set phenotypes

Nominal interaction statistics are reported for every cis variant of every tier-2 phenotype. The central workflow combines them across strata and cohorts by inverse-variance weighting. `af`, `b_gi`, `b_gi_se` and `tstat_gi` are stored as float32.

`gxe.tier2: pairs` keeps only the flagged lead variants. It is a fallback for sites that cannot store or share the full table. It is not the default, because the central analysis then cannot test variants the site did not flag.

Phenotypes or lead variants that cannot be tested are listed in `H_gxe_tier2_untested`, with one of these reasons:

- `phenotype_absent`;
- `variant_absent`;
- `outside_cis_window`;
- the exposure status.

`top_variant` marks the variant with the smallest nominal p for each phenotype. `q_top_bh` is the BH q over these per-phenotype minima within the stratum, as in the atlas finalization (54). **It is not corrected for the number of variants in each window, so it is anti-conservative.** Read it as a descriptive summary only; the same applies to the tier-2 `n_sig` in H_status. Tier-2 inference is made centrally.

### Sensitivity: SNP × covariate refits (atlas 55)

Tier-1 hits with q ≤ `gxe.fdr` are refitted with extra genotype × covariate terms:

```
Y = G + E + G:E + Σ_j G:C_j + covariates
```

Each C_j is an enabled entry of `gxe.sensitivity`:

- `maternal_age`, `parity`, `smoking` or `GA` read a metadata column, which also enters as a main effect;
- `ancestry_pcs` uses the genotype PCs already in the covariates, as interactions only.

All entries are off by default. The validator fails if an enabled covariate has no numeric metadata column. A hit is `stable` if the adjusted interaction changes by at most 20 % of the unadjusted one and stays on the same side of p = 0.05 (the atlas rule). `delta_b` and `rel_change` give the change.

### Run status

H_status has one row per stratum × modality × exposure × component (`tier1`, `tier2`, `sensitivity`). Each row has n, the number of covariates, the residual dof, the number tested and the number significant. Status values:

| Status | Meaning |
|---|---|
| `DONE` | the component ran |
| `SKIPPED` | the exposure failed `min_n` or `min_cases` in this stratum |
| `OFF` | sensitivity only: no covariate enabled |
| `BLOCKED` | sensitivity only: an enabled covariate's column is missing |

## What was not ported

| Atlas analysis | Reason |
|---|---|
| 56: transmitted vs non-transmitted alleles | needs parent–offspring trios, which most cohorts do not have |
| 06/50, 50a–50d: cell-type diagnostics (HC3 robust SEs, permutation calibration, EAS downsampling) | atlas-level diagnostics of the G×C model; they run centrally if needed, not at each site |
| HC3 standard errors in G×C | as above; the site fits use the atlas's classical OLS SEs |

## Caveats

- **Power.** Interaction tests need far more samples than main-effect tests. In strata near the minimum size, expect only very large interactions to reach significance. Their value is in the central meta-analysis, which is why every stratum reports full tier-2 statistics regardless of significance.
- **Phenotype permutation** in tier 1 is slightly conservative for the interaction null when a main genetic effect exists. tensorQTL's interaction mode uses the same scheme.
- **Exposure scale.** The z-scale is pooled within a cohort, not across cohorts. H_exposures reports the pooled mean and SD in raw units, so that b_GxE can be put on a common scale across cohorts centrally.
- **Splicing interactions.** Quantile normalization and the inverse normal transform across a cluster's introns attenuate a planted splicing interaction strongly. The test fixture therefore checks GxE recovery on an expression plant only, and checks the splicing plant only for its main effect and test-set flags.

## Test fixture

`tests/fixtures/site` (see `make_site_fixture.py`) plants:

- a G × Cytotrophoblast interaction on `ENSGFIX0000011`;
- a G × GA interaction on `ENSGFIX0000021` (planted variant chr22:20290948:A:G);
- a G × GA splicing plant on clu_31.

The restricted test file adds MDACC pairs that exercise every flag:

- a phenotype absent from the site;
- a variant absent from the site;
- a different MDACC lead for a site eGene;
- shared pairs;
- an MDACC-only pair.

`tests/scripts/test_site.sh` runs `-profile test` with GA and the maternal-age sensitivity enabled. `check_site_interactions.py` then checks, in EUR and EAS:

- the primary cell type of `ENSGFIX0000011`;
- the tier-1 q-value and sign for `ENSGFIX0000021`;
- the tier-2 top variant for `ENSGFIX0000021`;
- the test-set sources and untested reasons;
- H_status;
- the sensitivity refit.

It asserts only the planted pairs: with 60 samples per stratum, another pair can reach q ≤ 0.05 by chance.

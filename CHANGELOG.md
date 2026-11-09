# Changelog

All notable changes to this pipeline are recorded here. Versions follow [semantic versioning](https://semver.org/): a change in the export format (files A–J) or in the restricted-file version is a major version change.

## 1.0.0 (unreleased)

First release.

- **Run modes:** `smoke` (site validation with a GO / NO-GO verdict), `imputation_prep` (genotypes to per-chromosome GRCh38 VCFs for the TOPMed Imputation Server), `prepare_references` (site alignment and quantification indices), `site` (FASTQs and imputed genotypes to the shared export), `phenotyping_only`, `export` (rebuild the export from a finished site run) and `central` (MDACC meta-analysis).
- **Site analysis:** post-imputation genotype QC, ancestry assignment against the 1000 Genomes panel, nine RNA phenotype modalities against the HPLRv2 long-read placental annotation plus cell-type proportions, normalization and hidden-factor covariates, tensorQTL permutation, conditionally independent and nominal cis mapping per ancestry stratum, a combined-modality arm, in-sample LD, and G×C and GxE interaction tests.
- **Export:** files A–J with a mandatory leak scan ([Export specification](docs/EXPORT_SPEC.md)). Exports record the pipeline version (the export format is tied to the major version) and use restricted-file version 1.
- **Central analysis:** allele harmonization, fixed- and random-effects meta-analysis, feature-level ACAT, SuShiE multi-ancestry fine-mapping, colocalization (SuSiE-coloc, ColocBoost), G×C and GxE synthesis, TWAS inputs (OTTERS, FUSION) and PIP enrichment ([Central analysis](docs/CENTRAL_ANALYSIS.md)).
- **Distribution:** one offline image (software, public references, Java runtime, Nextflow, pipeline code and documentation), as a Docker image and an Apptainer `.sif`; `bootstrap` exports the launcher. Profiles for LSF, SLURM, PBS Pro, Torque, SGE and local runs, plus the built-in `seadragon` profile.
- **GPU:** the CUDA image uses PyTorch `cu132` wheels, so GPU nodes need a driver that supports CUDA 13.2.

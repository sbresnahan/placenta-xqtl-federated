# placenta-xqtl-federated

Federated multi-ancestry placental xQTL analysis. Each participating cohort runs this pipeline on its own computers and shares only summary-level files (A–J) with MD Anderson Cancer Center (MDACC). No individual-level genotype, expression or phenotype data leave the site.

The pipeline reproduces the analyses of the [placental xQTL atlas](https://github.com/sbresnahan/Placenta_xQTL_atlas) (modules 01–07) at each site:

- genotype preparation for the TOPMed Imputation Server, post-imputation QC and ancestry assignment;
- nine molecular phenotypes from placental RNA-seq against the long-read placental transcriptome (HPLRv2): gene and isoform expression, alternative TSS and polyadenylation, splicing, intron retention, RNA editing, RNA stability, and cell-type proportions;
- cis-xQTL mapping per ancestry group (tensorQTL);
- cell-type-interaction (G×C) and genotype-by-exposure (GxE) tests;
- LD export and assembly of the summary files.

MDACC then meta-analyses the exports across cohorts and ancestries (`--mode central`).

**Documentation:** <https://sbresnahan.github.io/placenta-xqtl-federated/> (also included inside the image for offline reading). The same pages are in [`docs/`](docs/).

> **Status: 1.0.0-dev.** Site validation (`--mode smoke`) and input validation work; the analysis modes are being added milestone by milestone (see [CHANGELOG.md](CHANGELOG.md)).

## Quick start

Everything runs from one self-contained image with no internet access needed. Sites need Linux, Apptainer (≥ 1.0) or Singularity (≥ 3.8), and optionally a batch scheduler (LSF, SLURM, PBS or SGE).

```bash
# 1. Get the image (once) and verify it
apptainer pull xqtl-federated_<version>.sif oras://ghcr.io/sbresnahan/xqtl-federated-sif:<version>
#    or download xqtl-federated_<version>.sif from the Zenodo record, then:
sha256sum -c xqtl-federated_<version>.sif.sha256

# 2. Export the launcher, pipeline code and documentation from the image
apptainer run xqtl-federated_<version>.sif bootstrap ./xqtl

# 3. Check that your cluster works with the pipeline (about 10 minutes)
./xqtl/run_xqtl.sh -profile seadragon,smoke          # or: -profile slurm,apptainer,smoke -c mysite.config

# 4. Build the alignment indices from the references inside the image (once)
./xqtl/run_xqtl.sh -profile seadragon -params-file cohort.yml --mode prepare_references

# 5. Run the site analysis
./xqtl/run_xqtl.sh -profile seadragon -params-file cohort.yml --mode site
```

> **Restricted file.** Site runs also need `xqtl-restricted-<version>.tar.gz`, a small file that holds the shared splicing clusters and the MDACC signal list. MDACC sends it by email after the data-sharing agreement is approved. If your site has not received it, request it from the scientific contact, **Sean T. Bresnahan (stbresnahan@mdanderson.org)**, stating your cohort, the pipeline version and the approved agreement. Set its path as `restricted:` in `cohort.yml`.

## Documentation map

| Page | For |
|---|---|
| [Collaboration brief](docs/COLLABORATION_BRIEF.md) | Cohort investigators: aims, analyses and the files requested |
| [HPC configuration](docs/HPC.md) | Running on your cluster: profiles, queues, memory, GPUs, the smoke test |
| [Seadragon validation](docs/SEADRAGON_VALIDATION.md) | Phase 0 results on MDACC's cluster |
| [Changelog](CHANGELOG.md) | Release history |

Setup, configuration, walkthrough, export specification, central analysis, troubleshooting and reference-data pages are added with their milestones.

## Repository layout

```text
main.nf, nextflow.config     entry point and configuration
conf/                        resource classes, executors, container engines, site profiles
workflows/, subworkflows/, modules/   pipeline code
bin/                         scripts run inside tasks (on PATH in the image)
assets/                      schema, head-job templates
containers/                  image recipe and software environment
docs/                        documentation (single source for the book)
book/                        bookdown builder for the documentation site
tests/                       unit tests and fixtures (synthetic or public data only)
scripts/                     developer and MDACC maintenance scripts
```

## Development

```bash
scripts/dev_env.sh /path/to/env        # conda/micromamba environment matching the image
export PATH=/path/to/env/align/bin:/path/to/env/xqtl/bin:/path/to/env/irfinder/bin:$PWD/bin:$PATH
python -m pytest tests/unit            # unit tests
scripts/check_configs.sh               # every profile resolves
nextflow run main.nf -profile local,dev,smoke --refs_dir tests/fixtures/refs_mini \
    --restricted tests/fixtures/xqtl-restricted-test.tar.gz
Rscript book/render_book.R             # documentation site into book/_book
```

## Licence and citation

GPL-3.0 (see [LICENSE](LICENSE)). Embedded third-party software and reference data keep their own licences. Please cite the placental xQTL atlas and the tools listed in the collaboration brief.

Contact: Sean T. Bresnahan, MD Anderson Cancer Center, Bhattacharya Lab (stbresnahan@mdanderson.org).

# Reference data and software provenance

Everything a site run reads, apart from the cohort's own data and the restricted file, is inside the image. This page records where each embedded file comes from, how it is checked, and how MDACC builds and releases it.

## Embedded reference data

The image holds the public references at `/opt/xqtl/refs`. The build downloads each file, verifies it against a sha256 pinned in [`assets/refs_sources.tsv`](../assets/refs_sources.tsv), processes it with the image's own tools, and writes a checksum list, `/opt/xqtl/refs/MANIFEST.tsv`. That list must equal [`assets/refs_manifest.tsv`](../assets/refs_manifest.tsv) for a release build to succeed. `apptainer run <sif> refs-check` repeats the comparison at any time, and `bootstrap` exports the list as `refs_MANIFEST.tsv`.

| Path in `/opt/xqtl/refs` | Contents | Source | Licence |
|---|---|---|---|
| `genome/GRCh38_no_alt_analysis_set.fa.gz` (+ `.fai`, `.gzi`) | GRCh38 no-alt analysis set, UCSC-style contig names; recompressed with bgzip and indexed | NCBI GenBank `GCA_000001405.15`, `seqs_for_alignment_pipelines.ucsc_ids` | Public domain (NCBI/GRC) |
| `liftover/hg19ToHg38.over.chain.gz` | hg19 → hg38 liftOver chain, for genotype arrays on GRCh37 (`--mode imputation_prep`) | UCSC Genome Browser downloads | Free use (UCSC downloads) |
| `annotation/HPLRv2.0.annotated.gtf.gz` | HPLRv2 long-read placental transcriptome annotation (SQANTI3-annotated) | MDACC public reference record (Zenodo) | CC-BY-4.0 (to be confirmed at upload) |
| `annotation/HPLRv2.0.annotated.PANTRY.gtf.gz` | The same annotation normalized for the phenotyping steps (as used by the atlas) | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `annotation/HPLRv2_classification.annotated.txt.gz` | SQANTI3 isoform classification table | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `txrevise/txrevise.{grp_1,grp_2}.{upstream,contained,downstream}.gff3.gz` | txrevise event annotations for alternative TSS and polyadenylation | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `music/placenta_music_reference.rds` | Placental single-cell reference for MuSiC deconvolution (built from public GEO data, GSE182381) | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `ancestry/1kGP_hg38_pca_panel.{pgen,pvar.zst,psam}` | 1000 Genomes high-coverage reference panel for ancestry PCA (see below) | MDACC record, derived from 1000 Genomes (IGSR) | IGSR: free, no restrictions |
| `ancestry/1kg_superpop.tsv` | 1000 Genomes sample → super-population labels | MDACC record, from IGSR | IGSR: free, no restrictions |
| `editing/rediportal_hg38.bed.gz` | Known RNA-editing sites (optional; see below) | MDACC record, from REDIportal | To be confirmed before inclusion |
| `genome/atlas_genome_contigs.tsv` | Contig names, lengths and MD5s of the genome FASTA the atlas used (provenance check only) | MDACC record | CC-BY-4.0 |
| `MANIFEST.tsv`, `README.txt` | sha256, size, description, source and licence of every file above | Written by the build | — |

The MDACC files come from one public Zenodo record, built from the atlas by `scripts/make_public_reference_record.sh`; its record number is pinned on the `# ZENODO_RECORD=` line of `assets/refs_sources.tsv`. Primary-source files (genome, chain) are downloaded from their original servers and pinned by sha256; their published MD5s were checked when they were pinned (genome `a08035b6a6e31780e96a34008ff21bd6`, chain `35887f73fe5e2231656504d1f6430900`).

### Genome

The atlas aligned to a cleaned copy of the GRCh38 no-alt analysis set (`..._cleaned_ready_for_salmon.fasta`). The image ships the original NCBI file. During the build, `fetch_references.sh` compares the contig names, lengths and sequence MD5s of the shipped genome with `atlas_genome_contigs.tsv` from the MDACC record and prints either `genome: identical to the atlas genome` or the contigs that differ. Any difference is reviewed before a release and recorded in the [changelog](../CHANGELOG.md).

### Ancestry PCA panel (deviation from the atlas)

The atlas pruned each cohort's genotypes (MAF ≥ 0.05, missingness ≤ 0.05, `--indep-pairwise 1000 100 0.1`) and intersected them with the full 1000 Genomes high-coverage panel (`chr:pos:ref:alt` IDs). Shipping the full panel would add many gigabytes to the image, so the federated panel is pre-filtered: biallelic autosomal A/C/G/T SNPs of the 1000 Genomes samples with super-population labels, duplicate IDs removed, MAF ≥ 0.05 in 1000 Genomes, then `--indep-pairwise 1000 100 0.2`. Sites still prune their own genotypes at the atlas setting (`1000 100 0.1`) and intersect with this panel. The looser panel threshold keeps more candidate SNPs than the cohort-side prune, so the intersection should be close to what the atlas used. Ancestry assignments from this panel are to be compared with the atlas assignments before the first release.

### RNA-editing sites (REDIportal)

The atlas used REDIportal sites to define the RNA-editing phenotype. REDIportal's terms for redistributing the site list have not been confirmed, so this file is **optional**: it is left out of the MDACC record (`INCLUDE_REDIPORTAL=no`) until the terms are confirmed. Without it, the RNA-editing phenotype will take a site-supplied BED of editing sites (a parameter added with the phenotyping release).

### Built at each site

STAR, Salmon and IRFinder-S indices are tens of gigabytes and depend on the aligner versions, so they are not shipped. `--mode prepare_references` builds them once per site from the genome and annotation above, into `site_refs:`.

### Not public: the restricted file

`xqtl-restricted-<version>.tar.gz` holds the frozen leafCutter cluster set (atlas clusters merged across cohorts and ancestries, junction coordinates only) and the MDACC signal list (atlas lead and conditionally independent cis-signals at q ≤ 0.05 for nine modalities and both ancestry groups, as variant and phenotype IDs). It contains no per-sample data, but it is not published. MDACC sends it by email to approved sites; request it from **Sean T. Bresnahan (stbresnahan@mdanderson.org)**, stating your cohort, the pipeline version and the approved agreement. The pipeline checks its version and internal checksums before use.

## Software

All software versions are pinned (see `containers/env/*.yml`, `containers/requirements-pip.txt`, `containers/install_r_extra.R`, `containers/build_irfinder.sh` and the `ARG` lines of `containers/Dockerfile`). The exact resolved versions are in the image at `/opt/xqtl/VERSIONS.tsv` and `/opt/xqtl/env_lock/`, printed by `apptainer run <sif> versions`, exported by `bootstrap`, and attached to each GitHub release.

Notes on how the image differs from a straight copy of the atlas environment:

- **Two conda environments.** Salmon 1.10.2 and the R/ICU stack need conflicting libraries, so the aligners (STAR, Salmon) live in a separate `align` environment. Both are on `PATH`; nothing changes for users.
- **IRFinder-S is compiled from source** (release v2.0.1, commit-verified) in a separate build stage.
- **Java and Nextflow** are a pinned Temurin 17 runtime and the self-contained Nextflow distribution, both checksum-verified at build time and exported to the host by `bootstrap`.
- **GPU variant.** `-cuda` images use PyTorch wheels built for CUDA 13.2 (`cu132`), the oldest CUDA build published for the pinned PyTorch release. GPU nodes need an NVIDIA driver that supports CUDA 13.2 (check the CUDA version shown by `nvidia-smi`). The CPU image runs everywhere and is the default.

## MDACC maintainer procedures

### Building the public reference record (once per reference release)

On seadragon, after confirming the paths marked `CHECK` in `scripts/packaging/seadragon.env`:

```bash
module load plink samtools
scripts/make_public_reference_record.sh --env scripts/packaging/seadragon.env --out /rsrch5/.../xqtl_public_refs --dry-run   # lists every input, or MISSING
scripts/make_public_reference_record.sh --env scripts/packaging/seadragon.env --out /rsrch5/.../xqtl_public_refs
```

Upload every file in the output directory to a new public Zenodo record (CC-BY-4.0 for the MDACC-derived files), publish it, then pin it and commit:

```bash
scripts/pin_reference_record.py --record <record id> --files /rsrch5/.../xqtl_public_refs/files.tsv
git commit -m "Pin public reference record <id>" assets/refs_sources.tsv
```

### Building the restricted file

```bash
scripts/make_restricted_bundle.sh --env scripts/packaging/seadragon.env --out /rsrch5/.../restricted --dry-run
scripts/make_restricted_bundle.sh --env scripts/packaging/seadragon.env --out /rsrch5/.../restricted
```

This writes `xqtl-restricted-<version>.tar.gz` and its `.sha256`. The archive is byte-reproducible: rebuilding from the same inputs gives the same checksum. Send it by institutional secure email, or as an encrypted archive with the password sent separately, to approved sites only.

### Building an image

`containers/build_image.sh` downloads and verifies the references, renders this documentation if needed, and runs `docker buildx build` with the right build contexts and arguments:

```bash
containers/build_image.sh                       # development build (unpinned rows allowed)
containers/build_image.sh --variant cuda        # GPU variant
containers/build_image.sh --release             # every row pinned; references must equal assets/refs_manifest.tsv
```

### Releasing a version

1. Pin the public reference record (above) and push to `main`. The `image` workflow builds `ghcr.io/sbresnahan/placenta-xqtl-federated:dev` and uploads `refs_MANIFEST.tsv` as an artifact.
2. Commit that file as `assets/refs_manifest.tsv` (first release, or whenever a reference changes).
3. Set `manifest.version` in `nextflow.config`, update the changelog, and push the tag `v<version>`. The workflow then builds the CPU and `-cuda` images in strict mode, pushes them, builds the `.sif`, pushes it to `ghcr.io/sbresnahan/xqtl-federated-sif:<version>` (pulled with `apptainer pull oras://...`), archives it on Zenodo, and attaches the `.sif` checksum, reference manifest and version list to the GitHub release.

One-time repository settings:

- **Secret `ZENODO_TOKEN`**: a Zenodo personal access token with `deposit:write` and `deposit:actions`. Optional repository variables: `ZENODO_CONCEPT_ID` (an earlier deposition of the image record, so every release becomes a new version under one concept DOI) and `ZENODO_PUBLISH=true` (publish automatically; otherwise the deposition is left as a draft to review on zenodo.org).
- **Package visibility**: after the first push, set both GHCR packages (`placenta-xqtl-federated` and `xqtl-federated-sif`) to public under the repository's Packages settings, so sites can pull without logging in.
- **Documentation site**: Settings → Pages → Deploy from a branch → `bookdown` / (root). The `docs` workflow renders the book and pushes it to that branch.

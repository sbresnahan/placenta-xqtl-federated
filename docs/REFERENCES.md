# Reference data and software provenance

Everything a site run reads, apart from the cohort's own data and the restricted file, is inside the image. This page records where each embedded file comes from, how it is checked, and how MDACC builds and releases it.

## Embedded reference data

The image holds the public references at `/opt/xqtl/refs`. The build downloads each file, verifies it against a sha256 pinned in [`assets/refs_sources.tsv`](../assets/refs_sources.tsv), processes it with the image's own tools, and writes a checksum list, `/opt/xqtl/refs/MANIFEST.tsv`. That list must equal [`assets/refs_manifest.tsv`](../assets/refs_manifest.tsv) for a release build to succeed. `apptainer run <sif> refs-check` repeats the comparison at any time, and `bootstrap` exports the list as `refs_MANIFEST.tsv`.

| Path in `/opt/xqtl/refs` | Contents | Source | Licence |
|---|---|---|---|
| `genome/GRCh38_no_alt_analysis_set.fa.gz` (+ `.fai`, `.gzi`) | GRCh38 no-alt analysis set, UCSC-style contig names; recompressed with bgzip and indexed | NCBI GenBank `GCA_000001405.15`, `seqs_for_alignment_pipelines.ucsc_ids` | Public domain (NCBI/GRC) |
| `liftover/hg19ToHg38.over.chain.gz` | hg19 → hg38 liftOver chain, for genotype arrays on GRCh37 (`--mode imputation_prep`) | UCSC Genome Browser downloads | Free use (UCSC downloads) |
| `annotation/HPLRv2.0.annotated.gtf.gz` | HPLRv2 long-read placental transcriptome annotation (SQANTI3-annotated) | MDACC public reference record (Zenodo) | CC-BY-4.0 (to be confirmed at upload) |
| `annotation/HPLRv2.0.annotated.PANTRY.gtf.gz` | The same annotation normalized for the phenotyping steps (Pantry GTF form) | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `annotation/HPLRv2_classification.annotated.txt.gz` | SQANTI3 isoform classification table | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `txrevise/txrevise.{grp_1,grp_2}.{upstream,contained,downstream}.gff3.gz` | txrevise event annotations for alternative TSS and polyadenylation | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `music/placenta_music_reference.rds` | Placental single-cell reference for MuSiC deconvolution (built from public GEO data, GSE182381) | MDACC record | CC-BY-4.0 (to be confirmed at upload) |
| `ancestry/1kGP_hg38_pca_panel.{pgen,pvar.zst,psam}` | 1000 Genomes high-coverage reference panel for ancestry PCA (see below) | MDACC record, derived from 1000 Genomes (IGSR) | IGSR: free, no restrictions |
| `ancestry/1kg_superpop.tsv` | 1000 Genomes sample → super-population labels | MDACC record, from IGSR | IGSR: free, no restrictions |
| `editing/rediportal_hg38.bed.gz` | Known RNA-editing sites (optional; see below) | MDACC record, from REDIportal | Usable without restriction |
| `genome/mdacc_genome_contigs.tsv` | Contig names, lengths and MD5s of the genome FASTA used for the MDACC signal list (provenance check only) | MDACC record | CC-BY-4.0 |
| `MANIFEST.tsv`, `README.txt` | sha256, size, description, source and licence of every file above | Written by the build | — |

The MDACC files come from one public Zenodo record, built by `scripts/make_public_reference_record.sh`; its record number is pinned on the `# ZENODO_RECORD=` line of `assets/refs_sources.tsv`. Primary-source files (genome, chain) are downloaded from their original servers and pinned by sha256; their published MD5s were checked when they were pinned (genome `a08035b6a6e31780e96a34008ff21bd6`, chain `35887f73fe5e2231656504d1f6430900`).

### Genome

The MDACC signal list and splicing clusters in the restricted file were derived on a cleaned copy of the GRCh38 no-alt analysis set. The image ships the original NCBI file. During the build, `fetch_references.sh` compares the contig names, lengths and sequence MD5s of the shipped genome with `mdacc_genome_contigs.tsv` from the MDACC record and prints either `genome: identical to the MDACC genome` or the contigs that differ. Any difference is reviewed before a release and recorded in the [changelog](../CHANGELOG.md).

### Ancestry PCA panel

Each site prunes its own genotypes (MAF ≥ 0.05, missingness ≤ 0.05, `--indep-pairwise 1000 100 0.1`) and intersects them with the 1000 Genomes high-coverage panel (`chr:pos:ref:alt` IDs). Shipping the full panel would add many gigabytes to the image, so the shipped panel is pre-filtered: biallelic autosomal A/C/G/T SNPs of the 1000 Genomes samples with super-population labels, duplicate IDs removed, MAF ≥ 0.05 in 1000 Genomes, then `--indep-pairwise 1000 100 0.2`. The looser panel threshold keeps more candidate SNPs than the cohort-side prune, so the intersection stays close to an intersection with the full panel.

### RNA-editing sites (REDIportal)

REDIportal sites define the RNA-editing phenotype. The REDIportal paper states that its tables are downloadable and usable without restriction ([Picardi et al. 2017, Nucleic Acids Research](https://pmc.ncbi.nlm.nih.gov/articles/PMC5210607/)), so the MDACC record includes the preprocessed hg38 site list (`INCLUDE_REDIPORTAL=yes`). The file stays optional in the pipeline: if it is absent, `rna_editing` is skipped with a warning. A site can also pass its own BED6 of editing sites with `--editing_sites` to `--mode prepare_references`.

### Built at each site

STAR, Salmon and IRFinder-S indices are tens of gigabytes and depend on the aligner versions, so they are not shipped. `--mode prepare_references` builds them once per site from the genome and annotation above, into `site_refs:` ([Phenotyping](PHENOTYPING.md#before-the-first-run-site-indices)).

### Not public: the restricted file

`xqtl-restricted-<version>.tar.gz` holds the frozen leafCutter cluster set (clusters pooled across the MDACC discovery cohorts and merged across ancestries, junction coordinates only) and the MDACC signal list (lead and conditionally independent cis-signals of the MDACC discovery analysis at q ≤ 0.05 for nine modalities and both ancestry groups, as variant and phenotype IDs). It contains no per-sample data, but it is not published. MDACC sends it by email to approved sites; request it from **Sean T. Bresnahan (stbresnahan@mdanderson.org)**, stating your cohort, the pipeline version and the approved agreement. The pipeline checks its version and internal checksums before use.

## Software

All software versions are pinned (see `containers/env/*.yml`, `containers/requirements-pip.txt`, `containers/install_r_extra.R`, `containers/build_irfinder.sh` and the `ARG` lines of `containers/Dockerfile`). The exact resolved versions are in the image at `/opt/xqtl/VERSIONS.tsv` and `/opt/xqtl/env_lock/`, printed by `apptainer run <sif> versions`, exported by `bootstrap`, and attached to each GitHub release.

Notes on the software layout:

- **Two conda environments.** Salmon 1.10.2 and the R/ICU stack need conflicting libraries, so the aligners (STAR, Salmon) live in a separate `align` environment. Both are on `PATH`; nothing changes for users.
- **IRFinder-S is compiled from source** (release v2.0.1, commit-verified) in a separate build stage.
- **Java and Nextflow** are a pinned Temurin 17 runtime and the self-contained Nextflow distribution, both checksum-verified at build time and exported to the host by `bootstrap`.
- **GPU variant.** `-cuda` images use PyTorch wheels built for CUDA 13.2 (`cu132`), the oldest CUDA build published for the pinned PyTorch release. GPU nodes need an NVIDIA driver that supports CUDA 13.2 (check the CUDA version shown by `nvidia-smi`). The CPU image runs everywhere and is the default.

## MDACC maintainer procedures

### Building the public reference record (once per reference release)

Copy `scripts/packaging/mdacc.env.template` to `mdacc.env` (ignored by git) and fill in the input paths. Then, on seadragon:

```bash
module load plink samtools
scripts/make_public_reference_record.sh --env mdacc.env --out /rsrch5/.../xqtl_public_refs --dry-run   # lists every input, or MISSING
scripts/make_public_reference_record.sh --env mdacc.env --out /rsrch5/.../xqtl_public_refs
```

Upload every file in the output directory to a new public Zenodo record (CC-BY-4.0 for the MDACC-derived files), publish it, then pin it and commit:

```bash
scripts/pin_reference_record.py --record <record id> --files /rsrch5/.../xqtl_public_refs/files.tsv
git commit -m "Pin public reference record <id>" assets/refs_sources.tsv
```

### Building the restricted file

```bash
scripts/make_restricted_bundle.sh --env mdacc.env --out /rsrch5/.../restricted --dry-run
scripts/make_restricted_bundle.sh --env mdacc.env --out /rsrch5/.../restricted
```

This writes `xqtl-restricted-<version>.tar.gz` and its `.sha256`. The archive is byte-reproducible: rebuilding from the same inputs gives the same checksum. Send it by institutional secure email, or as an encrypted archive with the password sent separately, to approved sites only.

Building the image and releasing a version are described in [Testing and validation](TESTING.md#building-and-releasing).

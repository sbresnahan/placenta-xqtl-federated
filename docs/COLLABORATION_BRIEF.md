# Collaboration brief: proposed analyses and summary files requested from participating cohorts

**Project:** Multi-ancestry placental xQTL atlas

**Prepared for:** Participating placental cohort investigators, including PACEplacenta collaborators

**Coordinating institution:** MD Anderson Cancer Center (MDACC), Bhattacharya Lab

**Scientific contact:** Sean T. Bresnahan (<stbresnahan@mdanderson.org>)

**Analysis software:** MDACC-provided self-contained pipeline (first release in preparation)

## 1. Purpose of the collaboration

We are studying how fetal genetic differences affect RNA produced by the placenta and whether conditions during pregnancy, such as maternal glucose levels or BMI, change those effects. We will combine results from cohorts with different genetic ancestries to identify placental genes and RNA forms that may be relevant to birth outcomes and childhood metabolic health. \[1, 2\]

For cohorts that cannot transfer participant-level data to MDACC, each cohort runs the analyses on its own computing system. It reviews the summary files and sends MDACC only those approved under its data-sharing rules. No participant-level genotypes, RNA-seq files, molecular measurements, or clinical records are sent to MDACC. Separately, each site must upload quality-controlled genotype files to the external TOPMed Imputation Server, as described in Section 5. Cohorts able to share de-identified individual-level data may discuss a separate arrangement with MDACC under an appropriate agreement.

This document describes the planned analyses, what each cohort would run locally, what MDACC would do with the shared results, and which files sites would send. The final section lists possible later analyses that would require discussion with collaborators before any additional work or data sharing.

## 2. Research questions and planned analyses

### 2.1 Identify genetic effects on different RNA measurements

We aim to identify genetic variants that affect placental RNA measurements, including changes that gene-level expression alone may miss. Using a placenta-specific reference built from long-read RNA sequencing (HPLRv2), participating sites will measure nine types of RNA features: **gene expression, isoform expression, isoform usage, splicing, intron retention, alternative transcription start sites, alternative polyadenylation, RNA editing, and RNA stability**. For each RNA feature, sites will test which nearby fetal genetic variants are associated with its measurement (cis-xQTL analysis), separately for each eligible ancestry group. The analyses will also identify independent genetic signals within each region. \[3, 4\]

Sites will estimate the proportions of major placental cell types using a public single-cell RNA reference. These estimates will be included in the statistical models to account for differences among samples. At significant xQTL signals from the site's own analysis, and at a list of xQTL signals supplied by MDACC, sites will also test whether the genetic association changes with the proportion of a cell type, such as trophoblast, endothelial, or stromal cells. \[5\]

MDACC will combine association results across cohorts and ancestry groups to identify consistent and varying genetic effects. We will also use each cohort's linkage disequilibrium (LD) results, which describe correlations among nearby genetic variants, to narrow the list of candidate causal variants across ancestries. \[6\]

### 2.2 Compare placental findings with birth and childhood health studies

We will integrate placental xQTL results with publicly available genome-wide association study (GWAS) summary statistics for birth weight, gestational duration, childhood BMI and anthropometric traits, glycemic phenotypes, and youth-onset type 2 diabetes. Planned analyses include statistical fine-mapping, colocalization to assess whether placental xQTL and GWAS signals may share causal variants, and transcriptome-wide association analyses (TWAS), which use genetic evidence to prioritize genes and RNA isoforms. Functional annotation will further characterize credible variants and candidate regulatory mechanisms. \[7–12\]

These analyses use association summaries and LD matrices shared by sites, together with published GWAS results. They do not require sites to release participant-level childhood outcomes or other participant-level clinical data to MDACC. More advanced methods, including multivariate isoform-level TWAS, will be added after we confirm that they work with summary statistics from separate cohorts.

### 2.3 Test whether pregnancy conditions modify genetic effects

Genotype-by-exposure (GxE) analysis will test whether the association between fetal genotype and placental RNA phenotypes changes with gestational exposures. The initial local workflow includes genotype-by-gestational-age interaction scans and, where sufficiently measured and permitted, interactions with **pre-pregnancy BMI, gestational diabetes mellitus (GDM), and oral glucose tolerance test (OGTT) glucose**. Exposure variables remain at the cohort site; only interaction statistics and aggregate exposure descriptors are requested. \[13\]

In the first release, sites will run interaction tests separately for each eligible ancestry group. MDACC will combine results from cohorts with comparable exposures and models. Planned extensions will examine more extensive interaction statistics, colocalization of exposure-modified signals with GWAS loci, and whether exposure-dependent regulatory annotations are enriched for childhood metabolic trait heritability.

### 2.4 Interpret the findings

We plan to report placental xQTL associations across cohorts and ancestry groups, likely causal variants, comparisons across types of RNA measurements, genes and isoforms prioritized at GWAS loci, and genetic effects that vary with gestational exposures. Public results will report combined findings. Cohort-specific results will be shared only as agreed with each participating cohort.

## 3. Analyses performed by participating sites

Each cohort runs analyses separately for each genetic ancestry group with at least 50 samples. Optional analyses depend on the available data and exposures. The letters in Table 1 identify the types of summary files requested from sites; Table 3 describes each file type.

**Table 1.** Analyses performed at each cohort and summary files shared with MDACC

| Task at the cohort | Analyses performed at the participating site | Summary files requested by MDACC |
|----|----|----|
| Genotype preparation and ancestry | Quality control (QC) before and after site-initiated imputation on the TOPMed Imputation Server; ancestry assignment using a public reference; sample and variant summaries for each ancestry group. \[14, 15\] | Aggregate QC and variant summaries (D, I). |
| RNA measurement | Measurement of nine types of RNA features using common HPLRv2 references; matching molecular features to the common reference; cell-type deconvolution. | Feature definitions and aggregate QC (E, I). |
| Data processing and statistical models | Cohort- and ancestry-specific normalization, sample/feature filtering, genotype principal components, latent factors, and covariate selection. | Model settings without participant-level covariate values (F). |
| Genetic association testing | Per-modality cis-association testing, permutation-based feature significance, and conditionally independent signals; nominal variant-by-feature statistics. | Association and signal tables (A, B). |
| Effects by cell type or pregnancy exposure | Genotype-by-cell-type proportion at the site's significant xQTL signals and an MDACC-supplied signal list; ancestry-stratified, transcriptome-wide genotype-by-exposure scans and follow-up interaction tests at the same signals where relevant covariates exist. | Interaction statistics and aggregate exposure descriptors (G, H). |
| LD, quality checks, and file preparation | Cohort-derived LD matrices; allele-frequency and imputation summaries; local quality checks; file list and checksums, and automated checks for sample identifiers. | LD, QC, and file documentation (C, D, I, J). |

All sites use the same set of reference RNA features so their results can be compared. Splicing analyses use a shared set of reference junction clusters; locally identified clusters are matched using genomic coordinates.

## 4. Analyses performed by MD Anderson

Before combining results, MDACC will check file format and integrity, RNA feature coordinates, allele coding, and allele frequencies across cohorts. Table 2 lists the planned analyses and the files they use. The letters refer to the file types defined in Table 3.

**Table 2.** Analyses at MDACC and the files used from each cohort

| Central analysis | Scientific purpose | Files used (Table 3) |
|----|----|----|
| Within- and cross-ancestry meta-analysis | Combine cohort estimates; assess heterogeneity and consistency; establish feature-level significance. | A, B, D, E, F |
| Cross-ancestry fine-mapping | Prioritize putative causal variants and credible sets using association statistics and LD estimated by the cohorts. | A, C |
| GWAS colocalization | Identify loci with statistical evidence of shared placental regulatory and birth or metabolic trait signals. | A, C plus public GWAS |
| Summary-statistic TWAS and functional interpretation | Prioritize genes and isoforms, associate regulatory predictions with GWAS, and evaluate functional annotation enrichment. | A, C plus public annotations and GWAS |
| Combined interaction analyses | Meta-analyze genotype-by-cell-type and genotype-by-exposure statistics; assess consistency across ancestry groups and cohorts. | G, H and relevant summary descriptors |

## 5. Data and computing needed at each cohort

The pipeline is intended for cohorts with fetal genotype data (SNP arrays, whole-exome sequencing, or whole-genome sequencing where suitable for downstream imputation and QC) and placental bulk RNA-seq. Broad transcript coverage, such as total RNA-seq, is important for isoform, splicing, transcript-end, and RNA-editing analyses; 3′-capture assays may not support the complete modality set. Data availability and sequencing design will be reviewed with each site before analysis.

Each site runs the analyses on its own research computing system using a supplied, self-contained software container that includes the reference data and runs without internet access, which keeps software and reference versions consistent across sites. After approval, MDACC also supplies by email a small restricted file containing the shared splicing reference clusters and the list of MDACC xQTL signals that every site tests for cell-type and exposure interactions. It uses the site's fetal genotypes, RNA sequencing data, and a limited sample metadata table. **Required local metadata** are an RNA-seq identifier, genotype identifier, fetal sex, and gestational age at delivery. **Optional local metadata** include pre-pregnancy BMI, GDM, OGTT glucose, maternal age, parity, and smoking; interaction modules are enabled only when relevant variables and sufficient samples are available.

Genotype imputation is required and must be performed through the **TOPMed Imputation Server**. Each participating site is responsible for uploading quality-controlled genotype files directly to the TOPMed Imputation Server in accordance with its terms of use and the site's data-governance requirements. \[15\]

## 6. Files requested from each cohort

Each site prepares summary files separately for each eligible ancestry group. Variants are identified by GRCh38 chromosome, position, and reference and alternate alleles. RNA features use the common annotation IDs and coordinates. The files contain no sample identifiers or participant-level records. The main requested results are variant-by-feature effect estimates (betas), standard errors (SEs), p-values, minor-allele frequencies (MAFs), minor-allele counts (MACs), sample sizes, and LD matrices. Table 3 gives the complete list of requested files and estimated sizes. Letters A–J match the file types in Tables 1 and 2.

**Table 3.** Summary files requested from each cohort

| File ID | File type | Contents and purpose | Level of detail | Approx. size per ancestry group |
|:--:|----|----|----|----|
| A | Variant-RNA association results | For each tested variant and RNA feature: effect allele, regression coefficient (beta), standard error (SE), p-value, effect-allele frequency, MAF, MAC, and analyzed sample size (N); used for meta-analysis and GWAS integration. | Variant x feature | 30–60 GB |
| B | Lead variants and independent signals | Best variant, permutation significance, q-value, and independent signals; supports feature-level evidence synthesis. | Feature or signal | \< 1 GB |
| C | Genetic correlation (LD) matrices | Pairwise variant correlations in genomic windows for quality-filtered common variants; supports fine-mapping, colocalization, and summary-based TWAS. | Variant x variant matrix per window | ~70 GB |
| D | Variant quality and frequency summaries | Allele frequency, imputation quality, Hardy–Weinberg statistics, and sample counts; supports harmonization and QC. | Variant | \< 1 GB |
| E | RNA feature descriptions | Molecular feature coordinates, grouping, and aggregate detection fraction; supports cross-site feature matching. | Feature | \< 1 GB |
| F | Model settings | Normalization method, covariate names, and latent-factor count; **no individual covariate values**. | Model | KB |
| G | Genetic effects by cell type | Interaction coefficients (betas), SEs, and p-values at the site's significant xQTL signals and the MDACC-supplied signal list. | Variant x feature x cell type | \< 1 GB |
| H | Genetic effects by pregnancy exposure | Genotype-by-exposure interaction betas, SEs, p-values, permutation results, and aggregate exposure summaries (n, mean, SD, or prevalence). | Feature or variant x feature x exposure | 1–5 GB |
| I | Group-level quality summaries | Ancestry-group counts, distributions of sequencing metrics, and aggregate cell-type proportions. | Summary metric | KB |
| J | File list and data dictionary | File paths, checksums, file sizes, row counts, data types, descriptions, and software and reference versions. | Shared file | KB |

We estimate approximately 100–150 GB of files per ancestry group, mostly association statistics (A) and LD matrices (C). Actual sizes depend on sample size, the number of RNA features, quality filters, and optional analyses. Each site can review the files, file list, and data dictionary before sharing them. Transfers require the cohort's approval and an agreed secure transfer method, such as Globus.

## 7. Data sharing and safeguards

**Files sent to MDACC exclude** individual genotypes or sequencing reads (FASTQ/BAM), participant identifiers, per-sample molecular or clinical data, exposure values, latent factors, or cell-type proportions. Summary results are prepared separately for ancestry groups with at least 50 analyzed samples.

MDACC will store approved shared files on access-controlled research computing infrastructure and use them for the collaborative placental xQTL analyses described here. Per-cohort files will not be redistributed or publicly released without the contributing cohort's agreement. Publicly released findings will combine results across cohorts. Authorship, acknowledgement, access, and any additional uses should be established with collaborators in advance.

## 8. Planned extensions beyond the first release

The analyses in Table 4 may be added after the first pipeline release. Participating cohorts would discuss and approve any additional work, including sample-size requirements, covariates, files to share, and data-governance requirements. MDACC would use summary results from each cohort where possible. Analyses requiring individual-level data from multiple cohorts would need separate approval and an appropriate data-transfer agreement.

**Table 4.** Possible analyses after the first pipeline release

| Planned analysis | Purpose | Possible additional site work or files |
|----|----|----|
| Sex-stratified and exposure-stratified xQTL mapping | Compare regulatory effects by fetal sex and metabolic exposure strata. | Stratified association and signal tables, subject to minimum cell sizes. |
| Expanded genotype-by-exposure association testing and colocalization | Identify interaction signals beyond the signals tested in the first release and assess overlap with childhood metabolic GWAS. | More extensive interaction summary statistics and matching LD. |
| Heritability partitioning and mediation from summary data | Test whether placental modality-specific and exposure-dependent annotations are enriched for GWAS heritability, using stratified LD score regression (S-LDSC), graphREML, or mediated expression score regression (MESC). \[16–18\] | MDACC analyses using existing summary results and public GWAS and LD resources. |
| Validation of summary-based TWAS models | Test whether centrally derived isoform models predict phenotypes in held-out local samples. | Aggregate out-of-sample prediction metrics, if jointly agreed. |
| Site-directed analyses using the placenta annotation | Support cohort interests such as variance QTLs (varQTLs), maternal eQTLs, and expression-methylation QTLs (eQTMs). | Can be performed entirely at the cohort site; no additional files need to be shared with MDACC. |

## References

1\. Peng, S. et al. Genetic regulation of the placental transcriptome underlies birth weight and risk of childhood obesity. PLoS Genetics 14, e1007799 (2018).

2\. Bhattacharya, A. et al. Placental genomics mediates genetic associations with complex health traits and disease. Nature Communications 13, 706 (2022).

3\. Bresnahan, S. T. et al. Long-read assembly reveals vast transcriptional complexity in the placenta associated with metabolic and endocrine function. Nature Communications 17, 4795 (2026). doi:10.1038/s41467-026-71303-4.

4\. Munro, D. et al. Multimodal analysis of RNA sequencing data powers discovery of complex trait genetics. Nature Communications 15, 10387 (2024).

5\. Wang, X., Park, J., Susztak, K., Zhang, N. R. & Li, M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nature Communications 10, 380 (2019).

6\. Lu, Z. et al. Improved multi-ancestry fine-mapping identifies cis-regulatory variants underlying molecular traits and disease risk. Preprint (2024). doi:10.1101/2024.04.15.24305836.

7\. EGG Consortium et al. Maternal and fetal genetic effects on birth weight and their relevance to cardio-metabolic risk factors. Nature Genetics 51, 804–814 (2019).

8\. Solé-Navais, P. et al. Genetic effects on the timing of parturition and links to fetal birth weight. Nature Genetics 55, 559–567 (2023).

9\. Kumasaka, N. et al. Genome-wide association study on longitudinal and cross-sectional traits of child health and development in a Japanese population. Preprint (2025). doi:10.1101/2025.02.18.25322458.

10\. Srinivasan, S. et al. The first genome-wide association study for type 2 diabetes in youth: the ProDiGY Consortium. Diabetes 70, 996–1005 (2021).

11\. Wallace, C. A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genetics 17, e1009440 (2021).

12\. Bhattacharya, A. et al. Isoform-level transcriptome-wide association uncovers genetic risk mechanisms for neuropsychiatric disorders in the human brain. Nature Genetics 55, 2117–2128 (2023).

13\. BIOS Consortium et al. A linear mixed-model approach to study multivariate gene–environment interactions. Nature Genetics 51, 180–186 (2019).

14\. Conomos, M. P., Miller, M. B. & Thornton, T. A. Robust inference of population structure for ancestry prediction and correction of stratification in the presence of relatedness. Genetic Epidemiology 39, 276–293 (2015).

15\. Taliun, D. et al. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program. Nature 590, 290–299 (2021).

16\. Gazal, S. et al. Linkage disequilibrium-dependent architecture of human complex traits shows action of negative selection. Nature Genetics 49, 1421–1427 (2017).

17\. Li, H., Kamath, T., Mazumder, R. et al. Improved heritability partitioning and enrichment analyses using summary statistics with graphREML. Nature Genetics (2026). doi:10.1038/s41588-026-02649-0.

18\. Yao, D. W., O’Connor, L. J., Price, A. L. & Gusev, A. Quantifying genetic effects on disease mediated by assayed gene expression levels. Nature Genetics 52, 626–633 (2020).

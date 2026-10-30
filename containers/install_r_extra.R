# R packages that are not on conda-forge/bioconda, installed at pinned versions/commits
# into the xqtl environment's R library. Used by containers/Dockerfile and scripts/dev_env.sh.
# Requires network at build time only; the image itself never needs it.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L),
        timeout = 600)
# bioconda installs GenomeInfoDbData through a post-link download that does not fail the
# environment build when bioconductor.org is unreachable; stop here with a clear message.
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE))
  stop("GenomeInfoDbData is missing: the bioconda post-link download failed; rebuild the environment")
pins <- list(
  # HCP hidden covariates (atlas 19/25 steps)
  Rhcpp = "mvaniterson/Rhcpp@a285e2baf379ef4ca21d14dc64733dd8dfc4f0d8",
  # MuSiC deconvolution (atlas 18 steps)
  MuSiC = "xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c"
)
for (p in names(pins)) {
  if (!requireNamespace(p, quietly = TRUE))
    remotes::install_github(pins[[p]], upgrade = "never", dependencies = FALSE, build_vignettes = FALSE)
}
# colocboost (central colocalization), CRAN release pinned
if (!requireNamespace("colocboost", quietly = TRUE) ||
    as.character(packageVersion("colocboost")) != "1.0.10")
  remotes::install_version("colocboost", version = "1.0.10", upgrade = "never")
for (p in c(names(pins), "colocboost")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat(sprintf("%s\t%s\n", p, as.character(packageVersion(p))))
}

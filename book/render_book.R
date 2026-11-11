# render_book.R -- build the documentation site into book/_book and check every link.
#   Rscript book/render_book.R
# Requires R with bookdown (>= 0.40) and pandoc. Used locally, by .github/workflows/docs.yml
# and by the image build (rendered book copied to /opt/xqtl/docs_html).
invisible(Sys.setlocale("LC_ALL", "C.UTF-8"))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grepl("^--file=", args)])))
repo <- normalizePath(file.path(here, ".."))
source(file.path(here, "prepare_book.R"))
res <- prepare_book(repo = repo, build = file.path(here, "_build"))

cfg <- readLines(file.path(repo, "nextflow.config"))
ver <- sub(".*'([^']+)'.*", "\\1", grep("^\\s*version\\s*=", cfg, value = TRUE)[1])
env <- new.env()
assign("PIPELINE_VERSION", ver, envir = env)

out <- file.path(here, "_book")
unlink(out, recursive = TRUE)
owd <- setwd(res$build)
bookdown::render_book("index.Rmd", output_format = "bookdown::gitbook", envir = env, quiet = TRUE)
setwd(owd)
file.create(file.path(out, ".nojekyll"))

# ---- link-preview image ----------------------------------------------------------------
# gitbook writes og:image / twitter:image from `url` + `cover-image` in index.Rmd but does not copy
# the file, and always uses the small "summary" X card. Copy the image, use the large card, and add
# its size and alt text so previews render on the first share.
img <- "social-preview.png"
if (!file.copy(file.path(here, img), file.path(out, img), overwrite = TRUE))
  stop("could not copy book/", img, " into the rendered book", call. = FALSE)
for (f in list.files(out, pattern = "\\.html$", full.names = TRUE)) {
  x <- readLines(f, warn = FALSE)
  x <- sub('<meta name="twitter:card" content="summary" />',
           '<meta name="twitter:card" content="summary_large_image" />', x, fixed = TRUE)
  i <- grep('<meta property="og:image"', x, fixed = TRUE)
  if (length(i) == 1)
    x <- append(x, c('  <meta property="og:image:width" content="1200" />',
                     '  <meta property="og:image:height" content="630" />',
                     '  <meta property="og:image:alt" content="placenta-xqtl-federated: federated multi-ancestry placental xQTL pipeline documentation" />'),
                after = i)
  writeLines(x, f)
}
if (length(grep("og:image:width", readLines(file.path(out, "index.html")))) != 1)
  stop("og:image meta tags missing from index.html (check url and cover-image in index.Rmd)", call. = FALSE)

# ---- link check on the rendered HTML ------------------------------------------------
html <- list.files(out, pattern = "\\.html$", full.names = TRUE)
ids <- lapply(html, function(f) {
  x <- paste(readLines(f, warn = FALSE), collapse = "\n")
  unique(regmatches(x, gregexpr('(?<=\\sid=")[^"]+', x, perl = TRUE))[[1]])
})
names(ids) <- basename(html)
bad <- character(0)
for (f in html) {
  x <- paste(readLines(f, warn = FALSE), collapse = "\n")
  hrefs <- unique(regmatches(x, gregexpr('(?<=href=")[^"]+', x, perl = TRUE))[[1]])
  hrefs <- hrefs[!grepl("^([a-zA-Z][a-zA-Z0-9+.-]*:|//)", hrefs)]          # external
  hrefs <- hrefs[!grepl("^(libs/|style\\.css|search_index\\.json|\\./?$)", hrefs)]
  for (h in hrefs) {
    page <- sub("#.*$", "", h); anchor <- if (grepl("#", h)) sub("^[^#]*#", "", h) else ""
    if (page == "") page <- basename(f)
    if (!(page %in% names(ids))) {
      if (!file.exists(file.path(out, page))) bad <- c(bad, sprintf("%s -> %s (missing page)", basename(f), h))
      next
    }
    if (nzchar(anchor) && !(utils::URLdecode(anchor) %in% ids[[page]]))
      bad <- c(bad, sprintf("%s -> %s (missing anchor)", basename(f), h))
  }
}
if (length(bad)) stop("broken links in rendered book:\n  ", paste(unique(bad), collapse = "\n  "), call. = FALSE)
message(sprintf("book rendered: %d pages in %s; all internal links resolve", length(html), out))

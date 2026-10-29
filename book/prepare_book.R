# prepare_book.R -- assemble the documentation site sources from docs/*.md.
#
# docs/*.md (and CHANGELOG.md) are the single source: they stay readable on GitHub,
# and this script copies them into book/_build/ in a fixed chapter order, adding what
# bookdown needs:
#   * chapter and section IDs: each heading gets {#<chapter>--<github-slug>} so that
#     anchors are unique across the book and match GitHub's anchors;
#   * cross-document links rewritten: HPC.md#queues-and-partitions -> #hpc--queues-and-partitions,
#     same-page #anchor -> #<chapter>--anchor, links to repository files -> GitHub URLs;
#   * part and appendix headers.
# It stops with an error listing every link or anchor that does not resolve.
#
# Usage: Rscript book/prepare_book.R            (normally called by render_book.R)

repo_url <- "https://github.com/sbresnahan/placenta-xqtl-federated"

# Chapter order. Chapters whose file does not exist yet are skipped (with a note);
# a link to a skipped chapter is reported as broken.
chapters <- data.frame(
  file  = c("docs/COLLABORATION_BRIEF.md",
            "docs/SETUP.md", "docs/HPC.md", "docs/CONFIGURATION.md", "docs/WALKTHROUGH.md",
            "docs/TROUBLESHOOTING.md",
            "docs/EXPORT_SPEC.md", "docs/CENTRAL_ANALYSIS.md",
            "docs/REFERENCES.md", "docs/SEADRAGON_VALIDATION.md", "CHANGELOG.md"),
  part  = c("Part I: For cohort investigators",
            rep("Part II: Running the pipeline", 5),
            rep("Part III: Outputs and analysis", 2),
            rep("APPENDIX", 3)),
  title = c("Collaboration brief", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

chapter_id <- function(path) gsub("_", "-", tolower(sub("\\.md$", "", basename(path))))

# GitHub heading slug: strip markdown, lowercase, keep letters/digits/space/-/_, spaces -> -
gh_slug <- function(text) {
  t <- gsub("\\{#[^}]*\\}\\s*$", "", text)
  t <- gsub("!?\\[([^]]*)\\]\\([^)]*\\)", "\\1", t, perl = TRUE)   # links -> text
  t <- gsub("[`*]", "", t)
  t <- gsub("<[^>]+>", "", t)
  t <- tolower(trimws(t))
  t <- gsub("[^\\p{L}\\p{N} _-]", "", t, perl = TRUE)
  gsub(" ", "-", t, fixed = TRUE)
}

is_fence <- function(line) grepl("^\\s*(```|~~~)", line)

# Pass 1: collect heading slugs per document
scan_headings <- function(lines) {
  in_code <- FALSE
  seen <- character(0)
  out <- data.frame(line = integer(0), level = integer(0), text = character(0), slug = character(0))
  for (i in seq_along(lines)) {
    if (is_fence(lines[i])) { in_code <- !in_code; next }
    if (in_code) next
    m <- regmatches(lines[i], regexec("^(#{1,6})\\s+(.*?)\\s*#*\\s*$", lines[i], perl = TRUE))[[1]]
    if (length(m) == 0) next
    base <- gh_slug(m[3]); s <- base; k <- 0
    while (s %in% seen) { k <- k + 1; s <- paste0(base, "-", k) }
    seen <- c(seen, s)
    out[nrow(out) + 1, ] <- list(i, nchar(m[2]), m[3], s)
  }
  out
}

prepare_book <- function(repo = normalizePath(file.path(dirname(sys_script()), "..")),
                         build = file.path(repo, "book", "_build")) {
  unlink(build, recursive = TRUE)
  dir.create(build, recursive = TRUE)
  present <- file.exists(file.path(repo, chapters$file))
  for (f in chapters$file[!present]) message("note: ", f, " not present yet; chapter skipped")
  ch <- chapters[present, ]
  ch$id <- vapply(ch$file, chapter_id, "")
  docs <- lapply(ch$file, function(f) readLines(file.path(repo, f), warn = FALSE, encoding = "UTF-8"))
  names(docs) <- ch$file
  heads <- lapply(docs, scan_headings)
  for (f in ch$file) if (!any(heads[[f]]$level == 1)) stop(f, ": needs one level-1 heading (# Title)")

  problems <- character(0)
  resolve_anchor <- function(target_file, anchor, src) {
    cid <- ch$id[ch$file == target_file]
    h <- heads[[target_file]]
    if (is.na(anchor) || anchor == "") return(paste0("#", cid))
    if (!(anchor %in% h$slug)) {
      problems <<- c(problems, sprintf("%s: anchor #%s not found in %s", src, anchor, target_file))
      return(paste0("#", cid))
    }
    if (anchor == h$slug[h$level == 1][1]) return(paste0("#", cid))
    paste0("#", cid, "--", anchor)
  }
  rewrite_link <- function(target, src) {
    if (grepl("^[a-zA-Z][a-zA-Z0-9+.-]*:", target)) return(target)          # http:, mailto:, ...
    if (startsWith(target, "#")) return(resolve_anchor(src, substring(target, 2), src))
    parts <- strsplit(target, "#", fixed = TRUE)[[1]]
    path <- parts[1]; anchor <- if (length(parts) > 1) parts[2] else NA
    rel <- normalizePath(file.path(repo, dirname(src), path), mustWork = FALSE)
    rel <- sub(paste0("^", normalizePath(repo), "/?"), "", rel)
    if (grepl("\\.md$", rel)) {
      if (rel %in% ch$file) return(resolve_anchor(rel, anchor, src))
      problems <<- c(problems, sprintf("%s: link to %s, which is not a chapter of the book", src, target))
      return(target)
    }
    if (!file.exists(file.path(repo, rel))) {
      problems <<- c(problems, sprintf("%s: link to missing file %s", src, target))
      return(target)
    }
    kind <- if (dir.exists(file.path(repo, rel))) "tree" else "blob"
    paste0(repo_url, "/", kind, "/main/", rel, if (!is.na(anchor)) paste0("#", anchor) else "")
  }

  rmd_files <- "index.Rmd"
  last_part <- ""
  for (j in seq_len(nrow(ch))) {
    f <- ch$file[j]; lines <- docs[[f]]; h <- heads[[f]]; cid <- ch$id[j]
    in_code <- FALSE
    for (i in seq_along(lines)) {
      if (is_fence(lines[i])) { in_code <- !in_code; next }
      if (in_code) next
      k <- match(i, h$line)
      if (!is.na(k)) {
        txt <- sub("\\s*\\{#[^}]*\\}\\s*$", "", h$text[k])
        if (h$level[k] == 1) {
          title <- if (!is.na(ch$title[j])) ch$title[j] else txt
          lines[i] <- sprintf("# %s {#%s}", title, cid)
          if (!is.na(ch$title[j]) && ch$title[j] != txt) lines[i] <- paste0(lines[i], "\n\n*", txt, "*")
        } else {
          lines[i] <- sprintf("%s %s {#%s--%s}", strrep("#", h$level[k]), txt, cid, h$slug[k])
        }
        next
      }
      # inline links [text](target "title"), skipping inline code spans
      line <- lines[i]
      codes <- regmatches(line, gregexpr("`+[^`]*`+", line))[[1]]
      for (q in seq_along(codes)) line <- sub(codes[q], sprintf("\001%d\001", q), line, fixed = TRUE)
      mm <- gregexpr("\\]\\(([^)\\s]+)((\\s+\"[^\"]*\")?)\\)", line, perl = TRUE)[[1]]
      if (mm[1] != -1) {
        st <- attr(mm, "capture.start")[, 1]; ln <- attr(mm, "capture.length")[, 1]
        for (q in rev(seq_along(st))) {
          old <- substr(line, st[q], st[q] + ln[q] - 1)
          line <- paste0(substr(line, 1, st[q] - 1), rewrite_link(old, f), substr(line, st[q] + ln[q], nchar(line)))
        }
      }
      for (q in seq_along(codes)) line <- sub(sprintf("\001%d\001", q), codes[q], line, fixed = TRUE)
      lines[i] <- line
    }
    part <- ch$part[j]
    if (part != last_part) {
      hdr <- if (part == "APPENDIX") "# (APPENDIX) Appendices {-}" else sprintf("# (PART\\*) %s {-}", part)
      lines <- c(hdr, "", lines)
      last_part <- part
    }
    out <- sprintf("%02d-%s.md", j, cid)
    writeLines(lines, file.path(build, out), useBytes = TRUE)
    rmd_files <- c(rmd_files, out)
  }

  if (length(problems)) {
    stop("broken documentation links:\n  ", paste(unique(problems), collapse = "\n  "), call. = FALSE)
  }
  for (f in c("index.Rmd", "_output.yml", "style.css"))
    file.copy(file.path(repo, "book", f), file.path(build, f), overwrite = TRUE)
  yml <- readLines(file.path(repo, "book", "_bookdown.yml"))
  yml <- c(yml, "rmd_files:", paste0("  - \"", rmd_files, "\""))
  writeLines(yml, file.path(build, "_bookdown.yml"))
  message(sprintf("prepared %d chapters in %s", nrow(ch), build))
  invisible(list(build = build, chapters = ch))
}

sys_script <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(f) else normalizePath("book/prepare_book.R", mustWork = FALSE)
}

if (sys.nframe() == 0L) prepare_book()

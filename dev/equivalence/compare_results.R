# Compare the outputs of two scenario runs (reference vs. new) and report, per
# scenario: whether numeric values are bit-identical, whether all.equal() holds
# (including attributes), the largest absolute/relative difference, and timings.
#
# Usage: Rscript compare_results.R <ref_dir> <new_dir> [scenario_regex]

args <- commandArgs(trailingOnly = TRUE)
ref_dir <- args[1]
new_dir <- args[2]
pattern <- if (length(args) >= 3) args[3] else NULL

# Flatten any result (data frame, list of matrices, nested list) to a named list of atomic leaves
flatten <- function(x, prefix = "") {
  if (is.data.frame(x)) {
    out <- list()
    for (n in names(x)) out[[paste0(prefix, "$", n)]] <- x[[n]]
    return(out)
  }
  if (is.list(x)) {
    out <- list()
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    for (i in seq_along(x)) out <- c(out, flatten(x[[i]], paste0(prefix, "$", nms[i])))
    return(out)
  }
  stats::setNames(list(x), prefix)
}

strip_selfref <- function(x) {
  if (is.data.frame(x)) { attr(x, ".internal.selfref") <- NULL; return(x) }
  if (is.list(x)) return(lapply(x, strip_selfref))
  x
}

compare_one <- function(a, b) {
  res <- list(same_class = identical(class(a), class(b)), identical_values = TRUE,
              all_equal = TRUE, max_abs = 0, max_rel = 0, note = "")
  if (inherits(a, "scenario_error") || inherits(b, "scenario_error")) {
    res$note <- "error in a run"; res$identical_values <- res$all_equal <- FALSE; return(res)
  }
  ae <- all.equal(strip_selfref(a), strip_selfref(b), check.attributes = TRUE)
  res$all_equal <- isTRUE(ae)
  if (!res$all_equal) res$note <- paste(head(ae, 3), collapse = " | ")
  fa <- flatten(a); fb <- flatten(b)
  if (!identical(names(fa), names(fb))) {
    res$identical_values <- FALSE
    res$note <- paste("structure differs:", paste(setdiff(union(names(fa), names(fb)), intersect(names(fa), names(fb))), collapse = ","))
    return(res)
  }
  for (n in names(fa)) {
    va <- fa[[n]]; vb <- fb[[n]]
    if (is.factor(va)) va <- as.character(va)
    if (is.factor(vb)) vb <- as.character(vb)
    if (!identical(unname(va), unname(vb))) {
      res$identical_values <- FALSE
      if (is.numeric(va) && is.numeric(vb) && length(va) == length(vb)) {
        d <- abs(va - vb); ok <- is.finite(d)
        if (any(ok)) {
          res$max_abs <- max(res$max_abs, max(d[ok]))
          rel <- d[ok] / pmax(abs(va[ok]), .Machine$double.xmin)
          res$max_rel <- max(res$max_rel, max(rel))
        }
        if (!identical(is.na(va), is.na(vb))) res$note <- paste(res$note, "NA pattern differs in", n)
      } else {
        res$note <- paste(res$note, "non-numeric difference in", n)
      }
    }
  }
  res
}

ref_files <- list.files(ref_dir, pattern = "\\.rds$")
if (!is.null(pattern)) ref_files <- ref_files[grepl(pattern, ref_files)]
t_ref <- read.csv(file.path(ref_dir, "_timings.csv"))
t_new <- if (file.exists(file.path(new_dir, "_timings.csv"))) read.csv(file.path(new_dir, "_timings.csv")) else NULL

rows <- list()
for (f in ref_files) {
  nm <- sub("\\.rds$", "", f)
  if (!file.exists(file.path(new_dir, f))) next
  a <- readRDS(file.path(ref_dir, f)); b <- readRDS(file.path(new_dir, f))
  r <- compare_one(a, b)
  rows[[nm]] <- data.frame(scenario = nm, identical = r$identical_values, all_equal = r$all_equal,
                           same_class = r$same_class, max_abs = r$max_abs, max_rel = r$max_rel,
                           t_ref = t_ref$seconds[match(nm, t_ref$scenario)],
                           t_new = if (is.null(t_new)) NA else t_new$seconds[match(nm, t_new$scenario)],
                           note = trimws(r$note), stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
tab$speedup <- round(tab$t_ref / tab$t_new, 1)
print(tab[, c("scenario", "identical", "all_equal", "same_class", "max_abs", "max_rel", "t_ref", "t_new", "speedup")], row.names = FALSE)
if (any(nzchar(tab$note))) { cat("\nNotes:\n"); for (i in which(nzchar(tab$note))) cat(" -", tab$scenario[i], ":", tab$note[i], "\n") }
write.csv(tab, file.path(new_dir, "_comparison.csv"), row.names = FALSE)

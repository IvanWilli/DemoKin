# Internal helpers for the fast (C++) path, `fast = TRUE`. Not used when fast = FALSE.
# Summary indicators of kin_full (count_living, mean_age, sd_age, count_dead,
# count_cum_dead, mean_age_lost) by Focal's age, kin (and sex) and period/cohort.
#
# This reproduces exactly what the original dplyr/tidyr pipeline returned:
# groups sorted by `agrupar`, sums accumulated in long double over the rows of
# kin_full in their original order, cumulative deaths across Focal's ages
# within kin (x sex) x period/cohort, and a tibble with the same columns.
kin_summary_exact <- function(kin_full, agrupar, agrupar_no_age_focal) {
  kin_full <- as.data.frame(kin_full)
  keys <- kin_full[agrupar]
  # stable radix order == C-locale sort of the group keys, rows keep their order within a group
  ord <- do.call(order, c(unname(keys), list(method = "radix")))
  n <- length(ord)
  if (n == 0) {
    return(dplyr::tibble(keys, count_living = numeric(0), mean_age = numeric(0), sd_age = numeric(0),
                         count_dead = numeric(0), count_cum_dead = numeric(0), mean_age_lost = numeric(0)))
  }
  # first row of each group (keys equal, or both missing), found in C++
  # without copying the sorted key columns
  first <- group_starts_cpp(unname(as.list(keys)), ord)
  groups <- as.data.frame(lapply(keys, function(k) k[ord[first]]), stringsAsFactors = FALSE,
                          optional = TRUE)
  # cumulative groups: same keys without age_focal, in order of first appearance
  cum_key <- do.call(paste, c(groups[agrupar_no_age_focal], sep = "\r"))
  cg <- match(cum_key, unique(cum_key))
  age_kin <- kin_full$age_kin
  if (!is.integer(age_kin)) age_kin <- as.numeric(age_kin)
  res <- kin_summary_cpp(as.numeric(kin_full$living), as.numeric(kin_full$dead), age_kin,
                         as.integer(ord), as.integer(first), length(first), as.numeric(groups$age_focal),
                         as.integer(cg - 1L), max(cg))
  dplyr::as_tibble(c(as.list(groups), res))
}

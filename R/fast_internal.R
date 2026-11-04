# Internal helpers for the fast (C++) path, `fast = TRUE`. Not used when fast = FALSE.
# Internal helpers shared by the kin functions: thin wrappers around the C++
# projection kernel and the builders that turn kin matrices into the long
# data frames the package returns.
#
# The data frame builders reproduce exactly (values, column order, row order,
# column types and class) what the original reshape (pivot / melt / dcast /
# merge) pipeline produced, without materialising the intermediate long tables.

kin_names_14 <- c("d", "gd", "ggd", "m", "gm", "ggm", "os", "ys",
                  "nos", "nys", "oa", "ya", "coa", "cya")

# Focal's ageing operator in age-only models: a shift within the living rows,
# so that phi[, i] is the unit vector of age i-1 and F %*% phi[, i] = F[, i].
# Passed to the kernel as a triplet list (the same non-zero entries as the
# dense matrix it replaces) and memoised: it is the same every year.
focal_shift <- local({
  cache <- list()
  function(n, nlive) {
    key <- paste(n, nlive)
    if (is.null(cache[[key]])) {
      r <- if (nlive > 1) 2:nlive else integer(0)
      cache[[key]] <<- list(i = as.integer(r), j = as.integer(r - 1L), x = rep(1, length(r)),
                            nrow = as.integer(n), ncol = as.integer(n))
    }
    cache[[key]]
  }
})

# One-sex projection: Focal's whole trip (prev = NULL) or one time step.
project_kin_1sex <- function(Ut, ft, pit, ages, prev = NULL) {
  n <- 2 * ages
  pit <- pit[1:ages]
  k <- kin_project_cpp(Ut, ft, ft, focal_shift(n, ages),
                       phi0 = c(1, numeric(n - 1)), m0 = pit,
                       pi_all = pit, pi_os = pit,
                       ncol = ages, nlive = ages, prev = prev)
  k[kin_names_14]
}

# Layouts: which matrix rows are living kin (in the row order of the output)
# and which rows hold their deaths, plus the id columns that label them.
layout_1sex <- function(age) {
  ages <- length(age)
  list(ids = data.frame(age_kin = age),
       live = seq_len(ages),
       dead = list(dead = ages + seq_len(ages)))
}

# Sort kin names like merge()/data.table did: plain ASCII, locale-independent.
sort_kin <- function(x) sort(x, method = "radix")

# Values of one kin matrix block in output row order: id rows (lay$live or a
# death row set) by Focal age, Focal age fastest.
block_values <- function(x, rows, cols) as.vector(t(x[rows, cols, drop = FALSE]))

# Time-invariant long data frame: one block per kin (sorted), rows ordered by
# the id columns then Focal's age; deaths relocated one Focal age earlier
# (NA at the last age), as merge(all.x = TRUE) did.
build_kin_df_ti <- function(kin_list, age, lay) {
  ages <- length(age)
  nid <- nrow(lay$ids)
  kin_sorted <- sort_kin(names(kin_list))
  nk <- length(kin_sorted)
  nblock <- nid * ages
  out <- data.frame(kin = rep(kin_sorted, each = nblock), stringsAsFactors = FALSE)
  # id columns indexed directly: indexing the ids data frame by repeated rows
  # (then dropping the row names) gave the same columns but spent its time in
  # make.unique() on the row names
  idx <- rep(rep(seq_len(nid), each = ages), nk)
  for (nm in names(lay$ids)) out[[nm]] <- lay$ids[[nm]][idx]
  out$age_focal <- rep(as.integer(age), times = nid * nk)
  out$living <- unlist(lapply(kin_sorted, function(k) block_values(kin_list[[k]], lay$live, seq_len(ages))), use.names = FALSE)
  for (dn in names(lay$dead)) {
    out[[dn]] <- unlist(lapply(kin_sorted, function(k) {
      x <- kin_list[[k]]
      as.vector(t(cbind(x[lay$dead[[dn]], -1, drop = FALSE], NA_real_)))
    }), use.names = FALSE)
  }
  out
}

# Time-varying long data table: `store` is a list by year (character names,
# ascending) of list(ages = selected Focal ages, kin = list of kin matrices
# restricted to those ages). Deaths at (year + step, age + step) are relocated
# to (year, age) when that cell was selected, NA otherwise.
# keep_years: if given, only rows of those years are built (the same rows, in
# the same order, as building all years and then subsetting by year); the
# other stored years still supply the deaths.
build_kin_df_tv <- function(store, year_step, lay, keep_years = NULL) {
  years <- as.integer(names(store))
  nid <- nrow(lay$ids)
  ny <- length(store)
  year_col <- kin_col <- af_col <- liv_col <- id_idx <- vector("list", ny)
  dead_cols <- lapply(lay$dead, function(d) vector("list", ny))
  for (yi in seq_len(ny)) {
    Y <- years[yi]
    if (!is.null(keep_years) && !(Y %in% keep_years)) next
    A <- store[[yi]]$ages
    nA <- length(A)
    nxt <- match(as.character(Y + year_step), names(store))
    # column (in next year's stored matrices) holding deaths for each selected age
    dcols <- if (is.na(nxt)) rep(NA_integer_, nA) else match(A + year_step, store[[nxt]]$ages)
    ok <- !is.na(dcols)
    kin_sorted <- sort_kin(names(store[[yi]]$kin))
    nk <- length(kin_sorted)
    nblock <- nid * nA
    year_col[[yi]] <- rep(Y, nblock * nk)
    kin_col[[yi]] <- rep(kin_sorted, each = nblock)
    id_idx[[yi]] <- rep(rep(seq_len(nid), each = nA), nk)
    af_col[[yi]] <- rep(A, times = nid * nk)
    liv_col[[yi]] <- unlist(lapply(kin_sorted, function(k) block_values(store[[yi]]$kin[[k]], lay$live, seq_len(nA))), use.names = FALSE)
    for (dn in names(lay$dead)) {
      dead_cols[[dn]][[yi]] <- unlist(lapply(kin_sorted, function(k) {
        dcol <- matrix(NA_real_, length(lay$dead[[dn]]), nA)
        if (any(ok)) dcol[, ok] <- store[[nxt]]$kin[[k]][lay$dead[[dn]], dcols[ok], drop = FALSE]
        as.vector(t(dcol))
      }), use.names = FALSE)
    }
  }
  out <- data.frame(year = unlist(year_col, use.names = FALSE), kin = unlist(kin_col, use.names = FALSE),
                    stringsAsFactors = FALSE)
  idx <- unlist(id_idx, use.names = FALSE)
  for (nm in names(lay$ids)) out[[nm]] <- lay$ids[[nm]][idx]
  out$age_focal <- unlist(af_col, use.names = FALSE)
  out$cohort <- out$year - out$age_focal
  out$living <- unlist(liv_col, use.names = FALSE)
  for (dn in names(lay$dead)) out[[dn]] <- unlist(dead_cols[[dn]], use.names = FALSE)
  data.table::setDT(out)
  out
}

# Selected Focal ages (ascending) per year, from output_period_cohort_combination()
selected_ages_by_year <- function(out_selected) {
  sp <- split(out_selected$age, out_selected$year)
  lapply(sp, function(a) sort(unique(a)))
}

# Re-group kin types (os + ys -> s, ...) by summing within groups.
#
# by = "sorted" (kin2sex()) reproduces exactly what kin2sex() did with
#   group_by(<keys>) %>% summarise_at(vars(c("living", starts_with("dea"))),
#                                     funs(sum)) %>% ungroup()
# * rows sorted by the keys, ascending, character keys in the C locale and NA
#   last (dplyr >= 1.1 group ordering); order(method = "radix") is stable, so
#   rows keep their original order within a group, as in dplyr;
# * each value column summed with R's sum() arithmetic (group_sum_cpp());
# * a tibble with the keys, then the summed columns.
# by = "first" (kin()) reproduces
#   summarise(living = sum(living), dead = sum(dead), .by = <keys>)
# * groups in order of first appearance;
# * a tibble if the input is one, a plain data.frame otherwise (as.data.frame()
#   in summarise.data.frame());
# * vals: the summed columns, by default living and every dea* column.
# dplyr evaluated sum() once per group in R, which is where kin() and kin2sex()
# spent most of their time on large outputs.
regroup_kin_sum <- function(df, keys, by = c("sorted", "first"), vals = NULL) {
  by <- match.arg(by)
  if (is.null(vals)) {
    vals <- names(df)[names(df) == "living" | startsWith(tolower(names(df)), "dea")]
    vals <- c("living", setdiff(vals, "living"))
  }
  n <- nrow(df)
  kc <- lapply(keys, function(k) df[[k]])
  if (n == 0L) {
    start <- integer(0); o <- integer(0)
  } else {
    o <- do.call(order, c(unname(kc), list(na.last = TRUE, method = "radix")))
    start <- group_starts_cpp(kc, o)
  }
  if (by == "first" && length(start)) {
    # same groups, reordered by the row where each first appears; within a
    # group, o keeps the original row order (the radix order is stable)
    size <- diff(c(start, n + 1L))
    g_ord <- order(o[start], method = "radix")
    o <- o[unlist(lapply(g_ord, function(g) start[g] + seq_len(size[g]) - 1L), use.names = FALSE)]
    start <- c(1L, cumsum(size[g_ord])[-length(g_ord)] + 1L)
  }
  first <- o[start]
  grp <- NULL
  sum_col <- function(x) {
    if (is.double(x) && is.null(attributes(x))) return(group_sum_cpp(x, o, start))
    # anything else (e.g. an all-NA logical column): sum() itself, per group
    if (is.null(grp)) grp <<- split(o, rep.int(seq_along(start), diff(c(start, n + 1L))))
    unname(vapply(grp, function(i) sum(x[i]), sum(x[0])))
  }
  out <- c(stats::setNames(lapply(kc, function(x) x[first]), keys),
           stats::setNames(lapply(vals, function(v) sum_col(df[[v]])), vals))
  out <- dplyr::as_tibble(out)
  if (by == "first" && !inherits(df, "tbl_df")) out <- as.data.frame(out)
  out
}

# MASS::ginv(diag(s)) without the SVD, bit for bit.
#
# ginv() takes La.svd() (LAPACK dgesdd) of the diagonal matrix. For a
# diagonal input every Householder reflection of the bidiagonal reduction is
# the identity (tau = 0), so the bidiagonal matrix is diag(s) with a zero
# superdiagonal, and the bidiagonal SVD (dbdsdc) only rescales: for n > 25
# (SMLSIZ) it divides by orgnrm = max|s| and multiplies back, computing
# (s * (1/orgnrm)) * orgnrm, which can differ from s in the last bit; for
# n <= 25 it leaves s unchanged. The singular values are the absolute values
# of that, the singular vectors carry only signs and a permutation, and
# ginv() keeps the values above sqrt(eps) times the largest. The result is
# diagonal: sign(s) / d where kept, 0 elsewhere. Verified against MASS::ginv()
# on random, wide-range, zero, negative, tied and near-tolerance inputs from
# n = 1 to 1010 (dev/regroup/test_ginv_diag.R); outside the range where that
# derivation holds (a scale beyond 1e-290 .. 1e290, or non-finite values) the
# SVD is used.
ginv_diag <- function(s) {
  s <- as.vector(s)
  n <- length(s)
  m <- if (n) max(abs(s)) else 0
  if (!n || !all(is.finite(s)) || (m != 0 && (m < 1e-290 || m > 1e290)))
    return(MASS::ginv(Matrix::diag(s, n)))
  d <- if (n > 25L && m != 0) abs((s * (1 / m)) * m) else abs(s)
  keep <- d > max(sqrt(.Machine$double.eps) * max(d), 0)
  out <- matrix(0, n, n)
  if (any(keep)) {
    i <- which(keep)
    out[cbind(i, i)] <- (1 / d[i]) * sign(s[i])
  }
  out
}

# MASS::ginv(make_diag(v)) for the call sites: ginv_diag() when v has at
# least two entries; otherwise the original call, since diag() of a single
# number is an identity matrix of that size, not a 1 x 1 matrix.
ginv_of_diag <- function(v, make_diag = diag) {
  v <- as.vector(v)
  if (length(v) >= 2L) ginv_diag(v) else MASS::ginv(make_diag(v))
}

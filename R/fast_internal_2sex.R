# Internal helpers for the fast (C++) path, `fast = TRUE`. Not used when fast = FALSE.
# Two-sex helpers (with or without causes of death): projection wrapper and
# row layouts. The matrix block structure follows Caswell (2022) and
# Caswell, Margolis & Verdery (2023): living rows are [females, males], death
# rows follow, by sex, and (for causes of death) by age with cause fastest.

# Focal's ageing operator for two sexes, as triplets for kin_project_cpp():
# ones on the subdiagonal of each sex block of the living rows. Memoised.
focal_shift_2sex <- local({
  cache <- list()
  function(ages, n) {
    key <- paste(ages, n)
    if (is.null(cache[[key]])) {
      r <- c(2:ages, ages + 2:ages)
      cache[[key]] <<- list(i = as.integer(r), j = as.integer(r - 1L), x = rep(1, length(r)),
                            nrow = as.integer(n), ncol = as.integer(n))
    }
    cache[[key]]
  }
})

project_kin_2sex <- function(Ut, Ft, Ft_star, pit, sex_focal, ages, prev = NULL) {
  agess <- 2 * ages
  n <- if (is.matrix(Ut)) nrow(Ut) else Ut$nrow   # dense matrix or triplet list
  pif <- pit[1:ages]
  pim <- pit[(ages + 1):agess]
  # G matrix moves focal by age: bdiag(G, G) with G the subdiagonal shift,
  # written directly, as a sparse triplet list (it is the same every year)
  Gt <- focal_shift_2sex(ages, n)
  # locate focal at age 0 depending sex
  phi0 <- numeric(n)
  phi0[ifelse(sex_focal == "f", 1, ages + 1)] <- 1
  # initial conditions on os and nos depend on focal sex
  pios <- if (sex_focal == "f") pif else pim
  k <- kin_project_cpp(Ut, Ft, Ft_star, Gt, phi0, m0 = pit,
                       pi_all = pif + pim, pi_os = pios,
                       ncol = ages, nlive = agess, prev = prev)
  k[kin_names_14]
}

# time invariant: rows sorted by (age_kin, sex_kin)
layout_2sex_ti <- function(age, causes = 0) {
  ages <- length(age)
  agess <- 2 * ages
  live <- as.vector(rbind(seq_len(ages), ages + seq_len(ages)))
  ids <- data.frame(age_kin = rep(age, each = 2), sex_kin = rep(c("f", "m"), ages), stringsAsFactors = FALSE)
  if (causes == 0) {
    dead <- list(dead = agess + live)
  } else {
    sex0 <- rep(0:1, ages)      # 0 female, 1 male, in `live` order
    a0 <- rep(age, each = 2)
    dead <- lapply(seq_len(causes), function(cc) agess + sex0 * ages * causes + a0 * causes + cc)
    names(dead) <- paste0("deadcause", seq_len(causes))
  }
  list(ids = ids, live = live, dead = dead)
}

# time varying: rows sorted by (sex_kin, age_kin); death columns in C-locale order (dcast)
layout_2sex_tv <- function(age, causes = 0) {
  ages <- length(age)
  agess <- 2 * ages
  live <- seq_len(agess)
  ids <- data.frame(sex_kin = rep(c("f", "m"), each = ages), age_kin = rep(age, 2), stringsAsFactors = FALSE)
  if (causes == 0) {
    dead <- list(dead = agess + live)
  } else {
    sex0 <- rep(0:1, each = ages)
    a0 <- rep(age, 2)
    dead <- lapply(seq_len(causes), function(cc) agess + sex0 * ages * causes + a0 * causes + cc)
    names(dead) <- paste0("deadcause", seq_len(causes))
    dead <- dead[sort(names(dead), method = "radix")]
  }
  list(ids = ids, live = live, dead = dead)
}

# as.matrix(Matrix::bdiag(...)) for dense blocks, without the sparse round
# trip: the blocks placed on the diagonal of a zero matrix. A vector is a
# one-column block, and one list argument is a list of blocks, as in bdiag().
# bdiag() goes through a sparse matrix, which drops zero entries, so -0 comes
# back as +0; the same is done here, so the result is identical.
bdiag_dense <- function(...) {
  x <- list(...)
  if (length(x) == 1L && is.list(x[[1]])) x <- x[[1]]
  x <- lapply(x, function(b) if (is.matrix(b)) b else matrix(b, ncol = 1L))
  nr <- vapply(x, nrow, integer(1)); nc <- vapply(x, ncol, integer(1))
  out <- matrix(0, sum(nr), sum(nc))
  r0 <- c(0L, cumsum(nr)); c0 <- c(0L, cumsum(nc))
  for (k in seq_along(x)) {
    if (nr[k] && nc[k]) out[r0[k] + seq_len(nr[k]), c0[k] + seq_len(nc[k])] <- x[[k]]
  }
  out[which(out == 0)] <- 0
  out
}

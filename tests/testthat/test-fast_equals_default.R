# `fast = TRUE` (C++ projection, direct table builders) must return the same
# object as `fast = FALSE` (the original R code): same values, attributes and
# classes. No fixtures: both modes run here, on small inputs.
#
# Bit-identity relies on the C++ kernel accumulating in the same order as R's
# `%*%` does with R's internal (reference) BLAS. With an external BLAS
# (OpenBLAS, MKL, Accelerate) R's own products are not reproducible to the last
# bit either, so the comparison falls back to a tight tolerance.

strip_selfref <- function(x) {
  if (is.data.frame(x)) { attr(x, ".internal.selfref") <- NULL; return(x) }
  if (is.list(x)) return(lapply(x, strip_selfref))
  x
}

internal_blas <- function() {
  b <- extSoftVersion()[["BLAS"]]
  !nzchar(b) || grepl("Rblas|libRblas", b)
}

expect_fast_equals_default <- function(f) {
  fast <- suppressWarnings(suppressMessages(f(TRUE)))
  slow <- suppressWarnings(suppressMessages(f(FALSE)))
  fast <- strip_selfref(fast); slow <- strip_selfref(slow)
  if (internal_blas()) {
    expect_identical(fast, slow)
  } else {
    expect_equal(fast, slow, tolerance = 1e-12)
  }
}

# ---- one-sex -----------------------------------------------------------------

test_that("kin(): time invariant, grouped kin types, selected ages, summary", {
  expect_fast_equals_default(function(fast)
    kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"],
        output_kin = c("m", "d", "s", "a", "gm"), output_age_focal = c(0, 30, 65),
        birth_female = 0.49, fast = fast))
})

test_that("kin(): time invariant without summary, and with pi given", {
  pi <- swe_pop[, "2015"] * swe_asfr[, "2015"]; pi <- pi / sum(pi)
  expect_fast_equals_default(function(fast)
    kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"], pi = pi, summary_kin = FALSE, fast = fast))
})

test_that("kin_time_invariant(): data frame and list outputs", {
  expect_fast_equals_default(function(fast)
    kin_time_invariant(p = swe_px[, "2015"], f = swe_asfr[, "2015"], output_kin = c("m", "os"), fast = fast))
  expect_fast_equals_default(function(fast)
    kin_time_invariant(p = swe_px[, "2015"], f = swe_asfr[, "2015"], list_output = TRUE, fast = fast))
})

test_that("kin_time_variant(): period output, cohort output and list output", {
  yr <- as.character(2005:2015)
  expect_fast_equals_default(function(fast)
    kin_time_variant(p = swe_px[, yr], f = swe_asfr[, yr], output_period = 2010,
                     output_kin = c("d", "m"), fast = fast))
  expect_fast_equals_default(function(fast)
    kin_time_variant(p = swe_px[, yr], f = swe_asfr[, yr], output_cohort = 2005,
                     output_kin = c("d", "gm"), list_output = TRUE, fast = fast))
})

test_that("kin(): time varying with cohort output, output_year_last and grouped kin", {
  yr <- as.character(2010:2018)
  expect_fast_equals_default(function(fast)
    kin(p = swe_px[, yr], f = swe_asfr[, yr], time_invariant = FALSE, output_cohort = 2010,
        output_year_last = TRUE, output_kin = c("s", "n", "m", "d"), fast = fast))
  expect_fast_equals_default(function(fast)
    kin(p = swe_px[, yr], f = swe_asfr[, yr], n = swe_pop[, yr], time_invariant = FALSE,
        output_period = c(2012, 2018), output_age_focal = c(0, 40, 80), fast = fast))
})

test_that("timevarying_kin(): one projection step", {
  p <- swe_px[, "2015"]; f <- swe_asfr[, "2015"]
  ages <- length(p)
  pkin <- suppressMessages(kin_time_invariant(p = p, f = f, list_output = TRUE))
  Ut <- Mt <- zeros <- matrix(0, ages, ages)
  Ut[row(Ut) - 1 == col(Ut)] <- p[-ages]; Ut[ages, ages] <- p[ages]; diag(Mt) <- 1 - p
  Ut <- rbind(cbind(Ut, zeros), cbind(Mt, zeros))
  ft <- matrix(0, 2 * ages, 2 * ages); ft[1, 1:ages] <- f * 0.49
  pit <- swe_pop[, "2015"] * f; pit <- pit / sum(pit)
  expect_fast_equals_default(function(fast) timevarying_kin(Ut, ft, pit, ages, pkin, fast = fast))
})

# ---- two-sex -----------------------------------------------------------------

test_that("kin2sex(): time invariant, male Focal, grouped and ungrouped kin", {
  expect_fast_equals_default(function(fast)
    kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
            sex_focal = "m", output_kin = c("d", "s", "c"), output_age_focal = c(20, 50), fast = fast))
  expect_fast_equals_default(function(fast)
    kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
            output_kin = c("s", "a", "n", "m", "ys"), output_age_focal = c(5, 35, 65), fast = fast))
})

test_that("kin_time_invariant_2sex(): list output with pif/pim given", {
  pif <- fra_asfr_sex[, "ff"] * (0.3 + seq_len(101) / 101); pif <- pif / sum(pif)
  pim <- fra_asfr_sex[, "fm"] * (0.5 + seq_len(101) / 101); pim <- pim / sum(pim)
  expect_fast_equals_default(function(fast)
    kin_time_invariant_2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
                            pif = pif, pim = pim, list_output = TRUE, fast = fast))
})

test_that("kin2sex(): time invariant with three causes of death", {
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101)
  Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  expect_fast_equals_default(function(fast)
    kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
            Hf = Hf3, Hm = Hm3, sex_focal = "m", output_kin = c("m", "d", "s"),
            output_age_focal = c(30, 70), fast = fast))
  expect_fast_equals_default(function(fast)
    kin_time_invariant_2sex_cod(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
                                Hf = Hf3, Hm = Hm3, output_kin = c("m", "os"), list_output = TRUE, fast = fast))
})

test_that("timevarying_kin_2sex(): one projection step", {
  pf <- fra_surv_sex[, "pf"]; pm <- fra_surv_sex[, "pm"]; ff <- fra_asfr_sex[, "ff"]; fm <- fra_asfr_sex[, "fm"]
  ages <- length(pf); agess <- 2 * ages; bf <- 0.49
  pkin <- suppressMessages(kin_time_invariant_2sex(pf, pm, ff, fm, list_output = TRUE))
  Uf <- Um <- Ff <- Fm <- zeros <- matrix(0, ages, ages)
  Uf[row(Uf) - 1 == col(Uf)] <- pf[-ages]; Uf[ages, ages] <- pf[ages]
  Um[row(Um) - 1 == col(Um)] <- pm[-ages]; Um[ages, ages] <- pm[ages]
  Ut <- as.matrix(rbind(cbind(Matrix::bdiag(Uf, Um), Matrix::bdiag(zeros, zeros)),
                        cbind(Matrix::bdiag(diag(1 - pf), diag(1 - pm)), Matrix::bdiag(zeros, zeros))))
  Ff[1, ] <- ff; Fm[1, ] <- fm
  Ft <- Ft_star <- matrix(0, 2 * agess, 2 * agess)
  Ft[1:agess, 1:agess] <- rbind(cbind(bf * Ff, bf * Fm), cbind((1 - bf) * Ff, (1 - bf) * Fm))
  Ft_star[1:agess, 1:ages] <- rbind(bf * Ff, (1 - bf) * Ff)
  pit <- c(ff / sum(ff), fm / sum(fm))
  expect_fast_equals_default(function(fast)
    timevarying_kin_2sex(Ut, Ft, Ft_star, pit, "f", ages, pkin, fast = fast))
})

test_that("kin2sex(): time varying, with and without causes of death", {
  skip_on_cran()
  yrs <- as.character(2000:2004)
  swe_pm <- swe_px[, yrs]^1.5
  swe_fm <- rbind(matrix(0, 5, length(yrs)), swe_asfr[-(97:101), yrs]) * 1.05
  dimnames(swe_fm) <- dimnames(swe_asfr[, yrs])
  expect_fast_equals_default(function(fast)
    kin2sex(pf = swe_px[, yrs], pm = swe_pm, ff = swe_asfr[, yrs], fm = swe_fm,
            nf = swe_pop[, yrs], nm = swe_pop[, yrs] * 0.98, time_invariant = FALSE, sex_focal = "m",
            output_kin = c("c", "d", "gm", "n", "s"), output_period = c(2002, 2004),
            output_age_focal = c(0, 10, 40, 70), fast = fast))
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101)
  Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  Hf_list <- stats::setNames(lapply(seq_along(yrs), function(i) Hf3 * (1 + i / 50)), yrs)
  Hm_list <- stats::setNames(lapply(seq_along(yrs), function(i) Hm3 * (1 + i / 60)), yrs)
  expect_fast_equals_default(function(fast)
    kin2sex(pf = swe_px[, yrs], pm = swe_pm, ff = swe_asfr[, yrs], fm = swe_fm,
            Hf = Hf_list, Hm = Hm_list, time_invariant = FALSE,
            output_cohort = 2000, output_kin = c("m", "s", "gm"), fast = fast))
  expect_fast_equals_default(function(fast)
    kin_time_variant_2sex_cod(pf = swe_px[, yrs], pm = swe_pm, ff = swe_asfr[, yrs], fm = swe_fm,
                              Hf = Hf_list, Hm = Hm_list, output_period = 2003,
                              output_kin = c("d", "m"), list_output = TRUE, fast = fast))
})

# ---- multistage --------------------------------------------------------------

test_that("kin_multi_stage(): parity model, data frame and list outputs", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  loadNamespace("Matrix")   # the svk_* data are Matrix objects
  expect_fast_equals_default(function(fast)
    kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs, birth_female = 1,
                    parity = TRUE, list_output = TRUE, fast = fast))
  expect_fast_equals_default(function(fast)
    kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs, birth_female = 0.49,
                    parity = FALSE, output_kin = c("m", "d", "os"), fast = fast))
})

test_that("kin_multi_stage_time_variant_2sex(): education (male Focal) and UK parity data", {
  skip_on_cran()
  skip_if_not_installed("Matrix")
  loadNamespace("Matrix")
  expect_fast_equals_default(function(fast)
    kin_multi_stage_time_variant_2sex(
      U_list_females = U_mat_fem_edu[1:3], U_list_males = U_mat_male_edu[1:3],
      F_list_females = F_mat_fem_edu[1:3], F_list_males = F_mat_male_edu[1:3],
      T_list_females = T_mat_fem_edu[1:3], T_list_males = T_mat_male_edu[1:3],
      birth_female = c(.49, .5), parity = FALSE, summary_kin = TRUE,
      sex_Focal = "Male", initial_stage_Focal = 2, output_kin = c("Focal", "d", "m", "coa"),
      output_years = 1:3, fast = fast))
  expect_fast_equals_default(function(fast)
    kin_multi_stage_time_variant_2sex(
      U_list_females = Female_parity_mortality_list_UK[1:2], U_list_males = Male_parity_mortality_list_UK[1:2],
      F_list_females = Female_parity_fert_list_UK[1:2], F_list_males = Female_parity_fert_list_UK[1:2],
      T_list_females = Parity_transfers_by_age_list_UK[1:2], T_list_males = Parity_transfers_by_age_list_UK[1:2],
      birth_female = 1 - 0.51, parity = TRUE, output_kin = c("d", "oa", "ys", "os"), summary_kin = TRUE,
      sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1965:1966, fast = fast))
})

test_that("commutation_matrix() equals matrixcalc::commutation.matrix()", {
  expect_identical(DemoKin:::commutation_matrix(6, 110), matrixcalc::commutation.matrix(6, 110))
  expect_identical(DemoKin:::commutation_matrix(3, 5), matrixcalc::commutation.matrix(3, 5))
})

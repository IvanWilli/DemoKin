# Scenario grid shared by the reference run (original DemoKin) and the
# comparison run (DemoKin_fast). Each scenario is a zero-argument closure that
# uses only the package's exported functions and bundled data, so the same
# code runs unchanged against either implementation.
#
# Inputs are built once (build_inputs()) and captured by the closures.

build_inputs <- function() {
  inp <- list()

  # one-sex (Sweden, 1900-2018)
  inp$p_2015 <- swe_px[, "2015"]
  inp$f_2015 <- swe_asfr[, "2015"]
  pi_2015 <- swe_pop[, "2015"] * swe_asfr[, "2015"]
  inp$pi_2015 <- pi_2015 / sum(pi_2015)
  inp$pi_mat <- t(t(swe_pop * swe_asfr) / colSums(swe_pop * swe_asfr))

  # two-sex (France 2012, time invariant)
  inp$fra_pf <- fra_surv_sex[, "pf"]
  inp$fra_pm <- fra_surv_sex[, "pm"]
  inp$fra_ff <- fra_asfr_sex[, "ff"]
  inp$fra_fm <- fra_asfr_sex[, "fm"]
  pif <- inp$fra_ff * c(0.3 + seq_len(101) / 101)
  inp$fra_pif <- pif / sum(pif)
  pim <- inp$fra_fm * c(0.5 + seq_len(101) / 101)
  inp$fra_pim <- pim / sum(pim)

  # two-sex time varying (vignette 1_4 construction)
  years <- ncol(swe_px)
  ages <- nrow(swe_px)
  inp$swe_pf <- swe_px
  inp$swe_ff <- swe_asfr
  inp$swe_pm <- swe_px^1.5
  swe_fm <- rbind(matrix(0, 5, years), swe_asfr[-((ages - 4):ages), ]) * 1.05
  dimnames(swe_fm) <- dimnames(swe_asfr)
  inp$swe_fm <- swe_fm
  inp$swe_nm <- swe_pop * 0.98

  # causes of death (vignette 1_4 construction, plus a 3-cause variant)
  Hf <- matrix(c(.5, 1), nrow = 2, ncol = 101)
  Hm <- matrix(c(.33, 1), nrow = 2, ncol = 101)
  Hf[, 1:50] <- Hm[, 1:50] <- 1
  inp$Hf2 <- Hf
  inp$Hm2 <- Hm
  Hf3 <- matrix(c(.5, 1, .25), nrow = 3, ncol = 101)
  Hm3 <- matrix(c(.33, 1, .4), nrow = 3, ncol = 101)
  Hf3[, 1:40] <- Hm3[, 1:40] <- 1
  inp$Hf3 <- Hf3
  inp$Hm3 <- Hm3
  yrs_cod <- as.character(1990:2018)
  inp$yrs_cod <- yrs_cod
  inp$Hf2_list <- stats::setNames(lapply(seq_along(yrs_cod), function(i) Hf * (1 + i / 100)), yrs_cod)
  inp$Hm2_list <- stats::setNames(lapply(seq_along(yrs_cod), function(i) Hm * (1 + i / 200)), yrs_cod)
  yrs_cod3 <- as.character(2000:2008)
  inp$yrs_cod3 <- yrs_cod3
  inp$Hf3_list <- stats::setNames(lapply(seq_along(yrs_cod3), function(i) Hf3 * (1 + i / 50)), yrs_cod3)
  inp$Hm3_list <- stats::setNames(lapply(seq_along(yrs_cod3), function(i) Hm3 * (1 + i / 60)), yrs_cod3)

  # multistage, two-sex time varying (test file construction)
  inp$joe_Ff <- (1 / 0.49) * svk_fxs
  inp
}

build_scenarios <- function(inp) {
  s <- list()
  q <- function(expr) suppressMessages(suppressWarnings(expr))

  # ---- one-sex, time invariant -------------------------------------------
  s$ti_1sex_default <- function() q(kin(p = inp$p_2015, f = inp$f_2015))
  s$ti_1sex_kinsel <- function() q(kin(p = inp$p_2015, f = inp$f_2015,
                                        output_kin = c("m", "d", "s", "a", "gm"),
                                        output_age_focal = c(0, 30, 65), birth_female = 0.49))
  s$ti_1sex_pi <- function() q(kin(p = inp$p_2015, f = inp$f_2015, pi = inp$pi_2015))
  s$ti_1sex_nosum <- function() q(kin(p = inp$p_2015, f = inp$f_2015, summary_kin = FALSE))
  s$ti_1sex_matrix_input <- function() q(kin(p = swe_px[, c("2010", "2011")], f = swe_asfr[, c("2010", "2011")]))
  s$ti_1sex_direct_df <- function() q(kin_time_invariant(p = inp$p_2015, f = inp$f_2015, output_kin = c("m", "os")))
  s$ti_1sex_direct_list <- function() q(kin_time_invariant(p = inp$p_2015, f = inp$f_2015, list_output = TRUE))

  # ---- one-sex, time varying ---------------------------------------------
  s$tv_1sex_cohort <- function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE,
                                        output_cohort = c(1910, 1960)))
  s$tv_1sex_period <- function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE,
                                        output_period = c(1950, 2000, 2018)))
  s$tv_1sex_period_n <- function() q(kin(p = swe_px, f = swe_asfr, n = swe_pop, time_invariant = FALSE,
                                          output_period = 2000, output_kin = c("m", "gm", "s")))
  s$tv_1sex_all_small <- function() q(kin(p = swe_px[, as.character(2000:2018)], f = swe_asfr[, as.character(2000:2018)],
                                           time_invariant = FALSE))
  s$tv_1sex_yearlast <- function() q(kin(p = swe_px[, as.character(2000:2018)], f = swe_asfr[, as.character(2000:2018)],
                                          time_invariant = FALSE, output_period = c(2017, 2018), output_year_last = TRUE))
  s$tv_1sex_bf_vec <- function() q(kin(p = swe_px[, as.character(1950:2018)], f = swe_asfr[, as.character(1950:2018)],
                                        time_invariant = FALSE, birth_female = c(0.49, 0.485), output_cohort = 1950))
  s$tv_1sex_pi_matrix <- function() q(kin(p = swe_px, f = swe_asfr, pi = inp$pi_mat, time_invariant = FALSE,
                                           output_cohort = 1980))
  s$tv_1sex_age_focal <- function() q(kin(p = swe_px, f = swe_asfr, time_invariant = FALSE,
                                           output_cohort = 1960, output_age_focal = c(10, 40), summary_kin = FALSE))
  s$tv_1sex_direct_list <- function() q(kin_time_variant(p = swe_px[, as.character(2000:2018)], f = swe_asfr[, as.character(2000:2018)],
                                                          output_cohort = 2000, list_output = TRUE))
  s$tv_1sex_direct_df <- function() q(kin_time_variant(p = swe_px[, as.character(2005:2015)], f = swe_asfr[, as.character(2005:2015)],
                                                        output_period = 2010, output_kin = c("d", "m")))

  # ---- two-sex, time invariant -------------------------------------------
  s$ti_2sex_default <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm, birth_female = .5))
  s$ti_2sex_male <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm, sex_focal = "m",
                                          output_kin = c("d", "s", "c"), output_age_focal = c(20, 50)))
  s$ti_2sex_pi <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm,
                                        pif = inp$fra_pif, pim = inp$fra_pim))
  s$ti_2sex_nosum <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm, summary_kin = FALSE))
  s$ti_2sex_direct_list <- function() q(kin_time_invariant_2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm,
                                                                 sex_focal = "m", list_output = TRUE))

  # ---- two-sex, time varying ---------------------------------------------
  s$tv_2sex_cohort <- function() q(kin2sex(pf = inp$swe_pf, pm = inp$swe_pm, ff = inp$swe_ff, fm = inp$swe_fm,
                                            time_invariant = FALSE, birth_female = .5, output_cohort = 1900))
  s$tv_2sex_period <- function() q(kin2sex(pf = inp$swe_pf, pm = inp$swe_pm, ff = inp$swe_ff, fm = inp$swe_fm,
                                            nf = swe_pop, nm = inp$swe_nm, sex_focal = "m",
                                            time_invariant = FALSE, output_period = c(1960, 2000)))
  # grouped kin types (s, c, a, n): the regrouping step of kin2sex(), with
  # all eight kin types and many output years, as applied work often asks
  s$tv_2sex_grouped_period <- function() q(kin2sex(pf = inp$swe_pf, pm = inp$swe_pm, ff = inp$swe_ff, fm = inp$swe_fm,
                                                    nf = swe_pop, nm = inp$swe_nm, time_invariant = FALSE,
                                                    output_kin = c("c", "d", "gd", "gm", "m", "n", "a", "s"),
                                                    output_period = 1990:2018))
  s$tv_2sex_grouped_cohort <- function() q(kin2sex(pf = inp$swe_pf, pm = inp$swe_pm, ff = inp$swe_ff, fm = inp$swe_fm,
                                                    nf = swe_pop, nm = inp$swe_nm, time_invariant = FALSE,
                                                    sex_focal = "m", output_kin = c("s", "a", "m", "os"),
                                                    output_cohort = c(1950, 1980)))
  s$ti_2sex_grouped <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm,
                                            output_kin = c("c", "d", "gd", "gm", "m", "n", "a", "s")))
  s$tv_2sex_all_small <- function() {
    yr <- as.character(2000:2018)
    q(kin2sex(pf = inp$swe_pf[, yr], pm = inp$swe_pm[, yr], ff = inp$swe_ff[, yr], fm = inp$swe_fm[, yr],
              time_invariant = FALSE))
  }
  s$tv_2sex_bf_vec_yearlast <- function() {
    yr <- as.character(2000:2018)
    q(kin2sex(pf = inp$swe_pf[, yr], pm = inp$swe_pm[, yr], ff = inp$swe_ff[, yr], fm = inp$swe_fm[, yr],
              time_invariant = FALSE, birth_female = c(.5, .49, .48), output_period = 2018, output_year_last = TRUE,
              output_kin = c("m", "a", "n")))
  }
  s$tv_2sex_direct_list <- function() {
    yr <- as.character(2010:2018)
    q(kin_time_variant_2sex(pf = inp$swe_pf[, yr], pm = inp$swe_pm[, yr], ff = inp$swe_ff[, yr], fm = inp$swe_fm[, yr],
                            output_cohort = 2010, list_output = TRUE))
  }

  # ---- causes of death ---------------------------------------------------
  s$ti_cod_default <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm, Hf = inp$Hf2, Hm = inp$Hm2))
  s$ti_cod_3causes_male <- function() q(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm, Hf = inp$Hf3, Hm = inp$Hm3,
                                                 sex_focal = "m", output_kin = c("m", "d", "s"), output_age_focal = c(30, 70)))
  s$ti_cod_direct_list <- function() q(kin_time_invariant_2sex_cod(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm,
                                                                    Hf = inp$Hf2, Hm = inp$Hm2, list_output = TRUE))
  s$tv_cod_cohort <- function() {
    yr <- inp$yrs_cod
    q(kin2sex(pf = inp$swe_pf[, yr], pm = inp$swe_pm[, yr], ff = inp$swe_ff[, yr], fm = inp$swe_fm[, yr],
              Hf = inp$Hf2_list, Hm = inp$Hm2_list, time_invariant = FALSE, output_cohort = 1990))
  }
  s$tv_cod_period <- function() {
    yr <- inp$yrs_cod3
    q(kin2sex(pf = inp$swe_pf[, yr], pm = inp$swe_pm[, yr], ff = inp$swe_ff[, yr], fm = inp$swe_fm[, yr],
              Hf = inp$Hf3_list, Hm = inp$Hm3_list, time_invariant = FALSE, output_period = c(2003, 2005),
              output_kin = c("m", "s", "gm")))
  }

  # ---- multistage (Caswell 2020) -----------------------------------------
  s$ms_parity_list <- function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                                    birth_female = 1, list_output = TRUE, parity = TRUE))
  s$ms_df_kinsel <- function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                                  birth_female = 1, parity = FALSE, output_kin = c("m", "d", "os")))
  s$ms_df_all <- function() q(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs,
                                               birth_female = 0.49, parity = TRUE))

  # ---- multistage, two-sex, time varying (Butterick) ----------------------
  s$joe_svk_1yr <- function() q(kin_multi_stage_time_variant_2sex(
    list(svk_pxs), list(svk_pxs), list(inp$joe_Ff), list(inp$joe_Ff), list(svk_Uxs), list(svk_Uxs),
    birth_female = 0.49, output_kin = NULL, parity = TRUE, summary_kin = FALSE,
    sex_Focal = "Female", initial_stage_Focal = 1, output_years = seq(1990, 1990)))
  s$joe_edu <- function() q(kin_multi_stage_time_variant_2sex(
    U_list_females = U_mat_fem_edu, U_list_males = U_mat_male_edu,
    F_list_females = F_mat_fem_edu, F_list_males = F_mat_male_edu,
    T_list_females = T_mat_fem_edu, T_list_males = T_mat_fem_edu,
    birth_female = 1 / (1.06 + 1), parity = FALSE, summary_kin = TRUE,
    sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1:15))
  s$joe_edu_male <- function() q(kin_multi_stage_time_variant_2sex(
    U_list_females = U_mat_fem_edu[1:4], U_list_males = U_mat_male_edu[1:4],
    F_list_females = F_mat_fem_edu[1:4], F_list_males = F_mat_male_edu[1:4],
    T_list_females = T_mat_fem_edu[1:4], T_list_males = T_mat_male_edu[1:4],
    birth_female = c(.49, .5), parity = FALSE, summary_kin = FALSE,
    sex_Focal = "Male", initial_stage_Focal = 2, output_kin = c("Focal", "d", "m", "coa"), output_years = 1:4))
  s$joe_parity_uk_3yrs <- function() q(kin_multi_stage_time_variant_2sex(
    U_list_females = Female_parity_mortality_list_UK[1:3], U_list_males = Male_parity_mortality_list_UK[1:3],
    F_list_females = Female_parity_fert_list_UK[1:3], F_list_males = Female_parity_fert_list_UK[1:3],
    T_list_females = Parity_transfers_by_age_list_UK[1:3], T_list_males = Parity_transfers_by_age_list_UK[1:3],
    birth_female = 1 - 0.51, parity = TRUE, output_kin = c("d", "oa", "ys", "os"), summary_kin = TRUE,
    sex_Focal = "Female", initial_stage_Focal = 1, output_years = 1965:1967))

  # ---- helpers -----------------------------------------------------------
  s$rename_kin_f <- function() q(rename_kin(kin(p = inp$p_2015, f = inp$f_2015)$kin_summary, sex = "f"))
  s$rename_kin_2sex <- function() q(rename_kin(kin2sex(inp$fra_pf, inp$fra_pm, inp$fra_ff, inp$fra_fm)$kin_summary, sex = "2sex"))
  s$apc_combination <- function() q(output_period_cohort_combination(output_cohort = c(1950, 1980), output_period = NULL,
                                                                     age = 0:100, years_data = 1900:2018))
  s$apc_combination_period <- function() q(output_period_cohort_combination(output_period = c(1950, 1980),
                                                                            age = 0:100, years_data = seq(1900, 2015, 5)))
  s
}

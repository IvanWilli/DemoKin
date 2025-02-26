
# Here I test the two-sex time-variant multi-stage function against caswell 2020. To do so I restrict the input years
# to only 1990 (so the model becomes invariant, and restrict the output results to females only)

# Technical note: to ensure that the two sex c(pi_f, pi_m) has a pi_f equal to the pi used in Caswell 2020 we need to guarantee
# that the spectral radius of the block structured A = U_proj + F_star comes from the upper-left block:

# A = [U_fem , 0 ; 0 , U_male] + [(1-alpha)*F_fem , 0 ; alpha*F_male , 0].
# If lambda spectral radius of A[1:n,1:n] then lambda*w[1:n] = A[1,1] %*% w[1:n] = (U_fem + (1-alpha)*F_fem) %*% w[1:n] then w[1:n]*(1/(1-alpha)) %*% F_fem[1,] is the same pi used in Caswell.
# However, if lambda comes from bottom right block of A then w[1:n]*F_fem[1,] is not the same as pi in Caswell 2020.
# Numerically check!


test_that("same output in multi_stage (caswell 2020)", {
  Tf <- svk_Uxs
  Tm <- svk_Uxs
  Ff <- svk_fxs
  Fm <- svk_fxs
  Ff <- (1/0.49)*Ff
  Fm <- (1/0.49)*Fm
  Uf <- svk_pxs
  Um <- svk_pxs
  H <- svk_Hxs

  joe_output <- kin_multi_stage_time_variant_2sex(list(Uf),
                                                  list(Um),
                                                  list(Ff),
                                                  list(Fm),
                                                  list(Tf),
                                                  list(Tm),
                                                  list(H),
                                                  birth_female = 0.49, ## svk_fxs already divided
                                                  output_kin = NULL,
                                                  parity = TRUE,
                                                  summary_kin = FALSE,
                                                  sex_Focal = "Female", ##  define Focal's sex at birth
                                                  initial_stage_Focal = 1, ## Define Focal's stage at birth
                                                  output_years = seq(1990, (1990)))


  ## Younger sisters
  jcmp_ys <- joe_output %>% dplyr::filter(sex_kin == "Female", group == "ys") %>%
    dplyr::select(age_focal, age_kin, stage_kin, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count) # Joe's
  hals_output_ys <- kin_svk1990_caswell2020$ys # Hal's
  hals_output_ys <- as.data.frame(hals_output_ys)
  colnames(hals_output_ys) <- seq(0,109,1)
  hals_output_ys$age_kin <- rep(seq(0, (110-1), 1), each = 6)
  hals_output_ys$stage_kin <- rep(seq(1, 6), 110)
  hcmp_ys <- hals_output_ys %>% reshape2::melt(id = c("age_kin","stage_kin")) %>%
    dplyr::mutate(age_focal = variable,
                  count = value) %>%
    dplyr::select(age_kin, stage_kin, age_focal, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count)

  ## Older sisters
  jcmp_os <- joe_output %>% dplyr::filter(sex_kin == "Female", group == "os") %>%
    dplyr::select(age_focal, age_kin, stage_kin, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count) # Joe's
  hals_output_os <- kin_svk1990_caswell2020$os # Hal's
  hals_output_os <- as.data.frame(hals_output_os)
  colnames(hals_output_os) <- seq(0,109,1)
  hals_output_os$age_kin <- rep(seq(0, (110-1), 1), each = 6)
  hals_output_os$stage_kin <- rep(seq(1, 6), 110)
  hcmp_os <- hals_output_os %>% reshape2::melt(id = c("age_kin","stage_kin")) %>%
    dplyr::mutate(age_focal = variable,
                  count = value) %>%
    dplyr::select(age_kin, stage_kin, age_focal, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count)

  ## Check equivalence
  expect_equal(jcmp_ys$count, hcmp_ys$count)
  expect_equal(jcmp_os$count, hcmp_os$count)


})

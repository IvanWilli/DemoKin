




test_that("same output in multi_stage (caswell 2020)", {
  load("data/svk_Uxs.rda")
  Tf <- unname(svk_Uxs)
  Tm <- unname(svk_Uxs)
  load("data/svk_fxs.rda")
  Ff <- unname(svk_fxs)
  Fm <- unname(svk_fxs)
  Ff <- (1/0.49)*Ff
  Fm <- (1/0.49)*Fm
  load("data/svk_pxs.rda")
  Uf <- unname(svk_pxs)
  Um <- unname(svk_pxs)
  load("data/svk_Hxs.rda")
  H <- unname(svk_Hxs)

  load("data/kin_svk1990_caswell2020.rda")

  joe_output <- kin_multi_stage_time_variant_2sex(list(Uf),
                                                  list(Um),
                                                  list(Ff),
                                                  list(Fm),
                                                  list(Tf),
                                                  list(Tm),
                                                  list(H),
                                                  birth_female = 0.49, ## svk_fxs already divided
                                                  output_kin = FALSE,
                                                  parity = TRUE,
                                                  summary_kin = FALSE,
                                                  sex_Focal = "Female", ##  define Focal's sex at birth
                                                  initial_stage_Focal = 1, ## Define Focal's stage at birth
                                                  seq(1990, (1990)))


  ## Younger sisters
  jcmp_ys <- joe_output %>% dplyr::filter(sex_kin == "Female", group == "ys") %>%
    dplyr::select(age_focal, age_kin, stage_kin, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count)
  hals_output_ys <- kin_svk1990_caswell2020$ys
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
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count)
  hals_output_os <- kin_svk1990_caswell2020$os
  hals_output_os <- as.data.frame(hals_output_os)
  colnames(hals_output_os) <- seq(0,109,1)
  hals_output_os$age_kin <- rep(seq(0, (110-1), 1), each = 6)
  hals_output_os$stage_kin <- rep(seq(1, 6), 110)
  hcmp_os <- hals_output_os %>% reshape2::melt(id = c("age_kin","stage_kin")) %>%
    dplyr::mutate(age_focal = variable,
                  count = value) %>%
    dplyr::select(age_kin, stage_kin, age_focal, count) %>%
    dplyr::transmute(age_focal = age_focal, age_kin = age_kin, stage_kin = stage_kin, count = count)




  expect_equal(jcmp_ys$count, hcmp_ys$count)
  expect_equal(jcmp_os$count, hcmp_os$count)


})

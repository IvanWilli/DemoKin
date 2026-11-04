# The `fast` argument: default FALSE, validated, present on every modelling function.

test_that("the default is the R implementation and equals fast = FALSE", {
  res_default <- suppressMessages(kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"], output_kin = "m"))
  res_slow <- suppressMessages(kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"], output_kin = "m", fast = FALSE))
  expect_identical(res_default, res_slow)
})

test_that("fast must be a single TRUE or FALSE", {
  expect_error(kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"], fast = NA), "`fast` must be TRUE or FALSE")
  expect_error(kin(p = swe_px[, "2015"], f = swe_asfr[, "2015"], fast = "yes"), "`fast` must be TRUE or FALSE")
  expect_error(kin2sex(fra_surv_sex[, "pf"], fra_surv_sex[, "pm"], fra_asfr_sex[, "ff"], fra_asfr_sex[, "fm"],
                       fast = c(TRUE, FALSE)), "`fast` must be TRUE or FALSE")
  expect_error(kin_multi_stage(U = svk_Uxs, f = svk_fxs, D = svk_pxs, H = svk_Hxs, fast = 1), "`fast` must be TRUE or FALSE")
})

test_that("every modelling function has fast = FALSE as its last argument", {
  fns <- c("kin", "kin2sex", "kin_multi_stage", "kin_multi_stage_time_variant_2sex",
           "kin_time_invariant", "kin_time_variant", "kin_time_invariant_2sex", "kin_time_variant_2sex",
           "kin_time_invariant_2sex_cod", "kin_time_variant_2sex_cod",
           "timevarying_kin", "timevarying_kin_2sex", "timevarying_kin_2sex_cod")
  for (fn in fns) {
    fm <- formals(get(fn, asNamespace("DemoKin")))
    expect_identical(names(fm)[length(fm)], "fast", info = fn)
    expect_false(fm$fast, info = fn)
  }
})

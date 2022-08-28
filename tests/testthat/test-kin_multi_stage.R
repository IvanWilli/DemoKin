# test multi-stage function against caswell 2020

test_that("same output in multi_stage (caswell 2020)", {
  demokin_svk1990_caswell2020 <- kin_multi_stage(U = svk_Uxs,
                                                 f = svk_fxs,
                                                 D = svk_pxs,
                                                 H = svk_Hxs, birth_female=1, list_output = TRUE)
  expect_equal(demokin_svk1990_caswell2020$d[1:(110*6),], kin_svk1990_caswell2020$d)
  expect_equal(demokin_svk1990_caswell2020$gd[1:(110*6),], kin_svk1990_caswell2020$gd)
  expect_equal(demokin_svk1990_caswell2020$ggd[1:(110*6),], kin_svk1990_caswell2020$ggd)
  expect_equal(demokin_svk1990_caswell2020$m[1:(110*6),], kin_svk1990_caswell2020$m)
  expect_equal(demokin_svk1990_caswell2020$gm[1:(110*6),], kin_svk1990_caswell2020$gm)
  expect_equal(demokin_svk1990_caswell2020$ggm[1:(110*6),], kin_svk1990_caswell2020$ggm)
  expect_equal(demokin_svk1990_caswell2020$os[1:(110*6),], kin_svk1990_caswell2020$os)
  expect_equal(demokin_svk1990_caswell2020$ys[1:(110*6),], kin_svk1990_caswell2020$ys)
})





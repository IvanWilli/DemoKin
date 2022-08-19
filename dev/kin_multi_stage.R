
# replicate multi-satage paper ---------------------------------------------------------------

# input and output for testing dowloaded from http://www.demographic-research.org/Volumes/Vol42/38/

# # libraries
# library(matrixcalc)
# library(Matrix)
# library(MASS)
# library(tidyverse)
# library(R.matlab)

# select years for testing
# years_run <- c(seq(1960,2014,4),2014)
#
# # apply demokin function (source function starting line ~ 60)
# kin_out_1960_2014 <- map(years_run, function(year){
#   print(year)
#   input_matrix <- readMat(paste0("tests/SVK_kinmats/SVKmats",year-10,".mat")) # Hal has en error in matlab code, is always reading 10 years before
#   U <- input_matrix$U
#   D <- input_matrix$D
#   f <- input_matrix$F
#   H <- input_matrix$H
#   kin_multi_stage(U, f, D, H, list_output = TRUE)
# })
#
# # read paper output for same years
# kin_paper_1960_2014 <- map(years_run, function(year){
#   print(year)
#   readMat(paste0("tests/SVK_kinout/SVKkinout",year,".mat"))
# })
#
# # give names both lists
# names(kin_out_1960_2014) <- years_run
# names(kin_paper_1960_2014) <- years_run
#
# # build data frame cheking same output
# check_results <- map_df(names(kin_paper_1960_2014), function(year){
#
#   print(year)
#   # how to read matlab files
#   d_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,1]
#   gd_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,2]
#   ggd_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,3]
#   m_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,4]
#   gm_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,5]
#   ggm_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,6]
#   os_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,7]
#   ys_paper <- kin_paper_1960_2014[[year]]$kinout[,,1]$allkin[,,8]
#
#   # output form demokin funtcion
#   d_demokin <- kin_out_1960_2014[[as.character(year)]]$d
#   gd_demokin <- kin_out_1960_2014[[as.character(year)]]$gd
#   ggd_demokin <- kin_out_1960_2014[[as.character(year)]]$ggd
#   m_demokin <- kin_out_1960_2014[[as.character(year)]]$m
#   gm_demokin <- kin_out_1960_2014[[as.character(year)]]$gm
#   ggm_demokin <- kin_out_1960_2014[[as.character(year)]]$ggm
#   os_demokin <- kin_out_1960_2014[[as.character(year)]]$os
#   ys_demokin <- kin_out_1960_2014[[as.character(year)]]$ys
#
#   # logical data frame
#   data.frame(
#     year = year,
#     d = all.equal(d_paper,d_demokin),
#     gd = all.equal(gd_paper,gd_demokin),
#     ggd = all.equal(ggd_paper,ggd_demokin),
#     m = all.equal(m_paper,m_demokin),
#     gm = all.equal(gm_paper,gm_demokin),
#     ggm = all.equal(ggm_paper,ggm_demokin),
#     os = all.equal(os_paper,os_demokin),
#     ys = all.equal(ys_paper,ys_demokin))
# })
#
# # see same results
# check_results

# end OK

# function ----------------------------------------------------------------

#' Estimate kin counts by age and stage in a time invariant framework

#' @description Implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019). Stages are implied in length of input lists.

#' @param U numeric. stage transition prob: list of w elements with s x s matrix.
#' @param f numeric. state-specific fertility: list of w elements with s x s matrix.
#' @param D numeric. age advance (depends if U includes survival or not): list of s elements with w x w matrix.
#' @param H numeric. assigns the offspring of individuals in stage j to the appropriate age class: list of s elements with w x w matrix.
#' @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See the `vignette` for all kin types.
#'
#' @return A data frame or a list with ...
#' @export
#'

kin_multi_stage <- function(U = NULL, f = NULL, D = NULL, H = NULL,
                               birth_female = 1,
                               output_kin = NULL,
                               list_output = FALSE){

  # stages and ages
  s <- length(D)
  ages <- length(U)
  age <- (1:ages)-1

  # build block matrix
  bbU <- as.matrix(bdiag(lapply(U, function(x) as.matrix(x[[1]]))))
  bbF <- as.matrix(bdiag(lapply(f, function(x) as.matrix(x[[1]])))) * birth_female
  bbD <- as.matrix(bdiag(lapply(D, function(x) as.matrix(x[[1]]))))
  bbH <- as.matrix(bdiag(lapply(H, function(x) as.matrix(x[[1]]))))

  # rearrange amtrix
  K <- commutation.matrix(s, ages)
  Ut <- t(K) %*% bbD %*% K %*% bbU
  ft <-  t(K) %*% bbH %*% K %*% bbF
  Gt <-  Ut%*% ginv(diag(colSums(Ut)))

  # stable distribution mothers: age x stage
  At <- Ut+ft
  A_decomp <- eigen(At)
  lambda <- as.double(A_decomp$values[1])
  wt <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
  pi <- wt*At[1,]/sum(wt*At[1,])

  # marginal mothers age
  Iom <- diag(1,ages, ages);
  ones <- t(rep(1,s))
  piage <- kronecker(Iom,ones) %*% pi

  # initial vectors
  phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,ages*s,ages)
  phi[1,1] = 1

  # momarray is an array with pit in each column
  momarray <- pi %*% matrix(1,1,ages)
  Iom = diag(1, ages)
  Is = diag(1, s)
  Isom = diag(1, s*ages)
  Z=Is;
  Z[1,1]=0;
  for(i in 1:ages){
    # imom = 1
    E <- Iom[,i] %*% t(Iom[i,]); # al cuadrado?
    momarray[,i] <- kronecker(E,Z) %*% momarray[,i]
  }
  # re-scale
  momarray <- momarray %*% ginv(diag(colSums(momarray)))

  # focal´s trip
  m[,1] = momarray %*% piage;
  for(i in 1:(ages-1)){
    phi[,i+1] = Gt %*% phi[,i]
    d[,i+1]   = Ut %*% d[,i] + ft %*% phi[,i]
    gd[,i+1]  = Ut %*% gd[,i] + ft %*% d[,i]
    ggd[,i+1]  = Ut %*% ggd[,i] + ft %*% gd[,i]
    m[,i+1]   = Ut %*% m[,i]
    ys[,i+1]  = Ut %*% ys[,i] + ft %*% m[,i]
    nys[,i+1] = Ut %*% nys[,i] + ft %*% ys[,i]
  }

  gm[,1] = m %*% piage
  for(i in 1:(ages-1)){
    gm[,i+1]  = Ut %*% gm[,i]
  }

  ggm[,1] = gm %*% piage
  for(i in 1:(ages-1)){
    ggm[,i+1]  = Ut %*% ggm[,i]
  }

  os[,1]  = d %*% piage
  nos[1:ages,1] = gd[1:ages,] %*% piage
  for(i in 1:(ages-1)){
    os[,i+1]  = Ut %*% os[,i]
    nos[,i+1] = Ut %*% nos[,i] + ft %*% os[,i]
  }

  oa[,1]  = os %*% piage
  ya[,1]  = ys %*% piage
  coa[,1] = nos %*% piage
  cya[,1] = nys %*% piage
  for(i in 1:(ages-1)){
    oa[,i+1]  = Ut %*% oa[,i]
    ya[,i+1]  = Ut %*% ya[,i]  + ft %*% gm[,i]
    coa[,i+1] = Ut %*% coa[,i] + ft %*% oa[,i]
    cya[,i+1] = Ut %*% cya[,i] + ft %*% ya[,i]
  }

  # get results
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                   nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list %>% keep(names(.) %in% output_kin)
  }

  # as data.frame
  kin <- map2(kin_list, names(kin_list),
              function(x,y){
                # y <- "m"
                # x <- kin_list[["m"]]
                out <- as.data.frame(x)
                colnames(out) <- age
                out %>%
                  mutate(kin = y,
                         age_kin = sort(rep(age,s)),
                         stage_kin = rep(1:s,ages),
                         alive = rep("yes",s*ages)) %>%
                  pivot_longer(c(-age_kin, -stage_kin, -kin, -alive), names_to = "age_focal", values_to = "count") %>%
                  mutate(age_focal = as.integer(age_focal))
              }
  ) %>%
    reduce(rbind)

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- kin
  }

  return(out)
}

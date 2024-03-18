#' Estimate kin counts in a time invariant framework for two-sex model.

#' @description Two-sex matrix framework for kin count and death estimates.This produces kin counts grouped by kin, age and sex of
#' each relatives at each Focal´s age. For example, male cousins from aunts and uncles from different sibling's parents
#' are grouped in one male count of cousins. This also produces kin deaths grouped by kin, age, sex of
#' each relatives at each Focal´s age, and cause of death.
#' @details See Caswell (2022) for details on formulas.
#' @param pf numeric. A vector of survival probabilities for females with same length as ages.
#' @param Hf numeric. A matrix with cause-specific hazards for females with rows as causes and columns as ages, being the name of each col the age.
#' @param ff numeric. A vector of age-specific fertility rates for females with same length as ages.
#' @param pm numeric. A vector of survival probabilities for males with same length as ages.
#' @param fm numeric. A vector of age-specific fertility rates for males with same length as ages.
#' @param Hm numeric. A matrix with cause-specific hazards for males with rows as causes and columns as ages, being the name of each col the age.
#' @param sex_focal character. "f" for female or "m" for male.
#' @param birth_female numeric. Female portion at birth.
#' @param pif numeric. For using some specific non-stable age distribution of childbearing for mothers (same length as ages). Default `NULL`.
#' @param pim numeric. For using some specific non-stable age distribution of childbearing for fathers (same length as ages). Default `NULL`.
#' @param output_kin character. kin to return, considering matrilineal names. For example "m" for parents, "d" for children, etc. See the `vignette` for all kin types.
#' @param list_output logical. Results as a list with `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
#'
#' @return A data frame with focal´s age, related ages and type of kin
#' (for example `d` is children, `oa` is older aunts/uncles, etc.), sex, alive and death. If `list_output = TRUE` then this is a list.
#' @export

# BEN: Added hazard matrices as inputs.
#      Assume that input of cause-specific mortality will be in terms of
#      matrices of cause-specific hazards for the two sexes (causes * ages).
#      Alternative: a matrix (causes * ages) containing the ratio mxi/mx.
kin_time_invariant_2sex_cod <- function(pf = NULL,
                                        pm = NULL,
                                        ff = NULL,
                                        fm = NULL,
                                        Hf = NULL,
                                        Hm = NULL,
                                        sex_focal = "f",
                                        birth_female = 1 / 2.04,
                                        pif = NULL,
                                        pim = NULL,
                                        output_kin = NULL,
                                        list_output = FALSE) {


  # global vars
  .<-sex_kin<-alive<-count<-living<-dead<-age_kin<-age_focal<-cohort<-year<-total<-mean_age<-count_living<-sd_age<-count_dead<-mean_age_lost<-indicator<-value<-NULL

  # same input length

  # BEN: Now we should also check the dimensions of the cause-specific hazard
  #      matrices.
  if(!all(length(pf)==length(pm), length(pf)==length(ff), length(pf)==length(fm),
          nrow(Hf)==nrow(Hm), ncol(Hf)==ncol(Hm), ncol(Hf)==length(pf))) stop("Number of age groups of p's, h's, and f's should match")

  # make matrix transition from vectors. Include death counts with matrix M
  age = 0:(length(pf)-1)
  ages = length(age)
  agess = ages * 2
  Uf = Um = Ff = Fm = Gt = matrix(0, nrow=ages, ncol=ages)

  # BEN: The zero matrix was deleted from line above and has
  #      to be made specific according to living/dead kin
  #      part of the block matrix Ut.
  causes <- nrow(Hf) # number of causes of death
  zeros_l <- matrix(0, nrow = ages, ncol = (causes*ages)) # zero matrix for living kin part
  zeros_d = matrix(0, nrow = (causes*ages), ncol = (causes*ages)) # zero matrix for death kin part

  Uf[row(Uf)-1 == col(Uf)] <- pf[-ages]

  # BEN: What is the purpose of the following line? By default it is zero due to
  #      how the matrix is created
  Uf[ages, ages] = Uf[ages]

  Um[row(Um)-1 == col(Um)] <- pm[-ages]
  Um[ages, ages] = Um[ages]

  # BEN: Building of M, matrix of cause-specific prob. of dying.
  #      Hence, M = H D(h_tilde)^{-1} D(q)
  #      where h_tilde are the summed hazards for each age, and
  #      q = 1 - p
  sum_hf <- t(rep(1, causes)) %*% Hf # h_tilde female
  sum_hm <- t(rep(1, causes)) %*% Hm # h_tilde male
  Mf <- Hf %*% solve(diag(c(sum_hf))) %*% diag(1-pf)
  Mm <- Hm %*% solve(diag(c(sum_hm))) %*% diag(1-pm)
  # Mm <- diag(1-pm)
  # Mf <- diag(1-pf)

  # BEN: In order to classify kin death by both cause and age at death,
  #      we need a mortality matrices M_hat of dimension
  #      ((causes*ages) * ages). See eq.12 in Caswell et al. (2024).
  # Store columns of M as a list of vectors
  Mf.cols <- lapply(1:ncol(Mf), function(j) return(Mf[,j]))
  Mm.cols <- lapply(1:ncol(Mm), function(j) return(Mm[,j]))
  # Create M_hat using the vectors as elements of the block diagonal
  Ut <- as.matrix(rbind(
    cbind(Matrix::bdiag(Uf, Um), Matrix::bdiag(zeros_l, zeros_l)),
    cbind(Matrix::bdiag(Matrix::bdiag(Mf.cols), Matrix::bdiag(Mm.cols)), Matrix::bdiag(zeros_d, zeros_d))))

  Ff[1,] = ff
  Fm[1,] = fm

  # BEN: Accounting for causes of death leads to have different dimensions
  #      in Ft and Ft_star.
  Ft <- Ft_star <- matrix(0, (agess + agess*causes), (agess + agess*causes))
  Ft[1:agess,1:agess] <- rbind(cbind(birth_female * Ff, birth_female * Fm),
                               cbind((1-birth_female) * Ff, (1-birth_female) * Fm))

  # mother and father do not reproduce independently to produce focal´s siblings. Assign to mother
  Ft_star[1:agess,1:ages] <- rbind(birth_female * Ff, (1-birth_female) * Ff)

  # parents age distribution under stable assumption in case no input
  if(is.null(pim) | is.null(pif)){
    A = Matrix::bdiag(Uf, Um) + Ft_star[1:agess,1:agess]
    A_decomp = eigen(A)
    lambda = as.double(A_decomp$values[1])
    w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
    wf = w[1:ages]
    wm = w[(ages+1):(2*ages)]
    pif = wf * ff / sum(wf * ff)
    pim = wm * fm / sum(wm * fm)
  }

  # initial count matrix (kin ages in rows and focal age in column)
  # BEN: Changed dimensions of lower part (dead kin) to account for death from causes.
  phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0, (agess + agess*causes), ages)

  # locate focal at age 0 depending sex
  sex_index <- ifelse(sex_focal == "f", 1, ages+1)
  phi[sex_index, 1] <- 1

  # G matrix moves focal by age
  G <- matrix(0, nrow=ages, ncol=ages)
  G[row(G)-1 == col(G)] <- 1

  # BEN: Changed dimensions
  Gt <- matrix(0, (agess + agess*causes), (agess + agess*causes))

  Gt[1:(agess), 1:(agess)] <- as.matrix(Matrix::bdiag(G, G))

  # focal´s trip
  # names of matrix count by kin refers to matrilineal as general reference
  m[1:(agess),1] = c(pif, pim)
  for(i in 1:(ages-1)){
    # i = 1
    phi[,i+1] = Gt %*% phi[,i]
    d[,i+1]   = Ut %*% d[,i]   + Ft %*% phi[,i]
    gd[,i+1]  = Ut %*% gd[,i]  + Ft %*% d[,i]
    ggd[,i+1] = Ut %*% ggd[,i] + Ft %*% gd[,i]
    m[,i+1]   = Ut %*% m[,i]
    ys[,i+1]  = Ut %*% ys[,i]  + Ft_star %*% m[,i]
    nys[,i+1] = Ut %*% nys[,i] + Ft %*% ys[,i]
  }

  gm[1:(agess),1] = m[1:(agess),] %*% (pif + pim)
  for(i in 1:(ages-1)){
    gm[,i+1]  = Ut %*% gm[,i]
  }

  ggm[1:(agess),1] = gm[1:(agess),] %*% (pif + pim)
  for(i in 1:(ages-1)){
    ggm[,i+1]  = Ut %*% ggm[,i]
  }

  os[1:(agess),1]  = d[1:(agess),] %*% pif
  nos[1:(agess),1] = gd[1:(agess),] %*% pif
  for(i in 1:(ages-1)){
    os[,i+1]  = Ut %*% os[,i]
    nos[,i+1] = Ut %*% nos[,i] + Ft %*% os[,i]
  }

  oa[1:(agess),1]  = os[1:(agess),] %*% (pif + pim)
  ya[1:(agess),1]  = ys[1:(agess),] %*% (pif + pim)
  coa[1:(agess),1] = nos[1:(agess),] %*% (pif + pim)
  cya[1:(agess),1] = nys[1:(agess),] %*% (pif + pim)
  for(i in 1:(ages-1)){
    oa[,i+1]  = Ut %*% oa[,i]
    ya[,i+1]  = Ut %*% ya[,i]  + Ft_star %*% gm[,i]
    coa[,i+1] = Ut %*% coa[,i] + Ft %*% oa[,i]
    cya[,i+1] = Ut %*% cya[,i] + Ft %*% ya[,i]
  }

  # get results
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                   nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list %>% purrr::keep(names(.) %in% output_kin)
  }

  # as data.frame
  kin <- purrr::map2(kin_list, names(kin_list),
                     function(x,y){

                       # BEN: Death take place in the same year and age!
                       #      I adapted the code
                       #      below such that it works with the new dimensions.

                       # reassign deaths to Focal experienced age
                       x[(agess+1):(agess + agess*causes),1:(ages-1)] <- x[(agess+1):(agess + agess*causes),2:ages]
                       x[(agess+1):(agess + agess*causes),ages] <- 0
                       out <- as.data.frame(x)
                       colnames(out) <- age
                       out %>%
                         # BEN: the matrices have different dimensions when
                         #      we accounf for causes of death so what follows
                         #      has been substantially changed.
                         dplyr::mutate(kin = y,
                                       age_kin = c(rep(age,2), rep(rep(age,each=causes),2)),
                                       sex_kin = c(rep(c("f", "m"),each=ages), rep(c("f", "m"),each=ages*causes)),
                                       alive = c(rep("living",2*ages), rep(paste0("deadcause",1:causes),2*ages))) %>%
                         tidyr::pivot_longer(c(-age_kin, -kin, -sex_kin, -alive), names_to = "age_focal", values_to = "count") %>%
                         dplyr::mutate(age_focal = as.integer(age_focal)) %>%
                         tidyr::pivot_wider(names_from = alive, values_from = count)
                     }
  ) %>%
    purrr::reduce(rbind)

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- kin
  }

  return(out)
}

## BEN: ========================================================================

# Checks

# No dead parent at birth: deadcausei=0 when age_focal==0
# ff # fertility starts at age 13
# kin |> filter(kin == "m", age_focal ==0, age_kin >= 12)
#
# # pi when age_focal==0 and age_kin when fx>0:
# kin |> filter(kin == "m", age_kin >= 13, age_focal ==0)
# pif[14:101]
#
# # mother dying from cause i at age x when focal is age==1 comes from nber of
# # living mother age x when focal is age==1 multiplied by (1-pf[x])*(1/3)
# kin |> filter(kin == "m", age_kin == 14, age_focal ==1)
# 0.000246 * ((1-pf[15])*(1/3)) # mother
# 0.0000486 * ((1-pm[15])*(1/3)) # father
#
# # Store to compare with kin_time_invariant_2sex.R
# saveRDS(
#   kin,
#   here(
#     "checks",
#     "output_time_invariant_2sex.rds"
#   )
# )


## =============================================================================

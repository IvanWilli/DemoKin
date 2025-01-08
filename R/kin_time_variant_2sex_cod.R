#' Estimate kin counts in a time variant framework (dynamic rates) in a two-sex framework (Caswell, 2022)

#' @description Two-sex matrix framework for kin count estimates with varying rates.
#' This produces kin counts grouped by kin, age and sex of each relatives at each Focal´s age.
#' For example, male cousins from aunts and uncles from different sibling's parents are grouped in one male count of cousins. This also produces kin deaths grouped by kin, age, sex of
#' each relatives at each Focal´s age, and cause of death.
#' @details See Caswell (2022) for details on formulas.
#' @param pf numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param pm numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param ff numeric. Same as pf but for fertility rates.
#' @param fm numeric. Same as pm but for fertility rates.
#' @param Hf numeric. A list where each list element (being the name of each list element the year) contains a matrix with cause-specific hazards for females with rows as causes and columns as ages, being the name of each col the age.
#' @param Hm numeric. A list where each list element (being the name of each list element the year) contains a matrix with cause-specific hazards for males with rows as causes and columns as ages, being the name of each col the age.
#' @param sex_focal character. "f" for female or "m" for male.
#' @param pif numeric. For using some specific age distribution of childbearing for mothers (same length as ages). Default `NULL`.
#' @param pim numeric. For using some specific age distribution of childbearing for fathers (same length as ages). Default `NULL`.
#' @param nf numeric. Same as pf but for population distribution (counts or `%`). Optional.
#' @param nm numeric. Same as pm but for population distribution (counts or `%`). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth. This multiplies `f` argument. If `f` is already for female offspring, this needs to be set as 1.
#' @param list_output logical. Results as a list with years elements (as a result of `output_cohort` and `output_period` combination), with a second list of `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
#' @return A data.frame with year, cohort, Focal´s age, related ages, sex and type of kin (for example `d` is daughter, `oa` is older aunts, etc.), including living and dead kin at that age and sex.
#' @export

# BEN: Added hazard matrices as inputs.
#      Assume that input of cause-specific mortality will be in terms of
#      matrices of cause-specific hazards for the two sexes (causes * ages).
#      Alternative: a matrix (causes * ages) containing the ratio mxi/mx.
kin_time_variant_2sex_cod <- function(pf = NULL, pm = NULL,
                                  ff = NULL, fm = NULL,
                                  Hf = NULL, Hm = NULL,
                                  sex_focal = "f",
                                  birth_female = 1/2.04,
                                  pif = NULL, pim = NULL,
                                  nf = NULL, nm = NULL,
                                  output_cohort = NULL, output_period = NULL, output_kin = NULL,
                                  list_output = FALSE){

  # global vars
  .<-living<-dead<-age_kin<-age_focal<-cohort<-year<-total<-mean_age<-count_living<-sd_age<-count_dead<-mean_age_lost<-indicator<-value<-NULL

  # same input length

  # BEN: Now we should also check the dimensions of the cause-specific hazard
  #      matrices.
  if(!all(dim(pf) == dim(pm), dim(pf) == dim(ff), dim(pf) == dim(fm),
          nrow(Hf)==nrow(Hm), ncol(Hf)==ncol(Hm), ncol(Hf)==nrow(pf),
          length(Hf)==length(Hm), length(Hm)==ncol(pf))) stop("Dimension of P's, F's, and H's should match")

  # data should be from same interval years
  years_data <- as.integer(colnames(pf))
  if(stats::var(diff(years_data))!=0) stop("Data should be for same interval length years. Fill the gaps and run again")

  # utils
  age          <- 0:(nrow(pf)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  agess        <- ages*2
  om           <- max(age)

  # BEN: The zero matrix was deleted from line above and has
  #      to be made specific according to living/dead kin
  #      part of the block matrix Ut.
  causes <- nrow(Hf[[1]]) # number of causes of death
  zeros_l <- matrix(0, nrow = ages, ncol = (causes*ages)) # zero matrix for living kin part
  zeros_d = matrix(0, nrow = (causes*ages), ncol = (causes*ages)) # zero matrix for death kin part

  # age distribution at child born
  Pif <- pif; no_Pif <- FALSE
  Pim <- pim; no_Pim <- FALSE
  if(is.null(pif)){
    if(!is.null(nf)){
      Pif <- t(t(nf * ff)/colSums(nf * ff))
    }else{
      Pif <- matrix(0, nrow=ages, ncol=n_years_data)
      no_Pif <- TRUE
    }
  }
  if(is.null(pim)){
    if(!is.null(nm)){
      Pim <- t(t(nm * fm)/colSums(nm * fm))
    }else{
      Pim <- matrix(0, nrow=ages, ncol=n_years_data)
      no_Pim <- TRUE
    }
  }

  # get lists of matrix
  Ul = Fl = Fl_star = list()
  kin_all <- list()
  pb <- progress::progress_bar$new(
    format = "Running over input years [:bar] :percent",
    total = n_years_data + 1, clear = FALSE, width = 60)

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # BEN: First load function at the end of script
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  for(t in 1:n_years_data){
    # t = 1
    Uf = Um = Fft = Fmt = Mm = Mf = Gt = matrix(0, nrow=ages, ncol=ages)
    Uf[row(Uf)-1 == col(Uf)] <- pf[-ages,t]
    Uf[ages, ages] = pf[ages,t]
    Um[row(Um)-1 == col(Um)] <- pm[-ages,t]
    Um[ages, ages] = pm[ages,t]

    # BEN: Building of M, matrix of cause-specific prob. of dying.
    #      Hence, M = H D(h_tilde)^{-1} D(q)
    #      where h_tilde are the summed hazards for each age, and
    #      q = 1 - p
    sum_hf <- t(rep(1, causes)) %*% Hf[[t]] # h_tilde female
    sum_hm <- t(rep(1, causes)) %*% Hm[[t]] # h_tilde male
    Mf <- Hf[[t]] %*% solve(diag(c(sum_hf))) %*% diag(1-pf[,t])
    Mm <- Hm[[t]] %*% solve(diag(c(sum_hm))) %*% diag(1-pm[,t])
    # Mm <- diag(1-pm[,t])
    # Mf <- diag(1-pf[,t])

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

    Ul[[as.character(years_data[t])]] <- Ut
    Fft[1,] = ff[,t]
    Fmt[1,] = fm[,t]

    # BEN: Accounting for causes of death leads to have different dimensions
    #      in Ft and Ft_star.
    Ft <- Ft_star <- matrix(0, (agess + agess*causes), (agess + agess*causes))

    Ft[1:agess,1:agess] <- rbind(cbind(birth_female * Fft, birth_female * Fmt),
                                 cbind((1-birth_female) * Fft, (1-birth_female) * Fmt))
    Ft_star[1:agess,1:ages] <- rbind(birth_female * Fft, (1-birth_female) * Fft)
    Fl[[as.character(years_data[t])]] <- Ft
    Fl_star[[as.character(years_data[t])]] <- Ft_star
    A = Matrix::bdiag(Uf, Um) + Ft_star[1:agess,1:agess]

    # stable assumption at start
    if (t==1){
      p1f <- pf[,1]; p1m <- pm[,1]
      f1f <- ff[,1]; f1m <- fm[,1]
      # time boundary for pi
      A_decomp = eigen(A)
      lambda = as.double(A_decomp$values[1])
      w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      wf = w[1:ages]
      wm = w[(ages+1):(2*ages)]
      pif1 = wf * ff[,t] / sum(wf * ff[,t])
      pim1 = wm * fm[,t] / sum(wm * fm[,t])

      # BEN: Add Hf and Hm
      H1f <- Hf[[1]]
      H1m <- Hm[[1]]

      # BEN: cod version !!!
      kin_all[[1]] <- kin_time_invariant_2sex_cod(pf = p1f, pm = p1m,
                                                  ff = f1f, fm = f1m,
                                                  pif = pif1, pim = pim1,
                                                  Hf = H1f, Hm = H1m,
                                                  birth_female = birth_female, list_output = TRUE)
    }

    # project pi
    if(no_Pim | no_Pif){
      w <- A %*% w
      wf <- w[1:ages]
      wm <- w[(ages+1):(2*ages)]
      Pif[,t] <- wf * ff[,t] / sum(wf * ff[,t])
      Pim[,t] <- wm * fm[,t] / sum(wm * fm[,t])
    }
    pit <- c(Pif[,t], Pim[,t])

    # kin for next year
    kin_all[[t+1]] <- timevarying_kin_2sex_cod(Ut=Ut, Ft=Ft, Ft_star=Ft_star, causes,
                                               pit=pit, sex_focal, ages, pkin=kin_all[[t]])
    pb$tick()
  }

  # filter years and kin that were selected
  names(kin_all) <- as.character(years_data)

  # combinations to return
  out_selected <- output_period_cohort_combination(output_cohort, output_period, age = age, years_data = years_data)

  possible_kin <- c("d","gd","ggd","m","gm","ggm","os","ys","nos","nys","oa","ya","coa","cya")
  if(is.null(output_kin)){
    selected_kin_position <- 1:length(possible_kin)
  }else{
    selected_kin_position <- which(possible_kin %in% output_kin)
  }

  # first filter
  kin_list <- kin_all %>%
    purrr::keep(names(.) %in% as.character(unique(out_selected$year))) %>%
    purrr::map(~ .[selected_kin_position])
  # long format
  message("Preparing output...")
  kin <- lapply(names(kin_list), FUN = function(Y){
    X <- kin_list[[Y]]
    X <- purrr::map2(X, names(X), function(x,y){
      # reassign deaths to Focal experienced age
      x[(agess+1):(agess + agess*causes),1:(ages-1)] <- x[(agess+1):(agess + agess*causes),2:ages]
      x[(agess+1):(agess + agess*causes),ages] <- 0
      x <- data.table::as.data.table(x)
      x$year <- Y
      x$kin <- y
      x$sex_kin <- c(rep(c("f", "m"),each=ages), rep(c("f", "m"),each=ages*causes))
      x$age_kin <- c(rep(age,2), rep(rep(age,each=causes),2))
      x$alive <- c(rep("living",2*ages), rep(paste0("deadcause",1:causes),2*ages))
      return(x)
    }) %>%
      data.table::rbindlist() %>%
      stats::setNames(c(as.character(age), "year","kin","sex_kin","age_kin","alive")) %>%
      data.table::melt(id.vars = c("year","kin","sex_kin","age_kin","alive"), variable.name = "age_focal", value.name = "count")
    X$age_focal = as.integer(as.character(X$age_focal))
    X$year = as.integer(X$year)
    X$cohort = X$year - X$age_focal
    X <- X[X$age_focal %in% out_selected$age[out_selected$year==as.integer(Y)],]
    X <- data.table::dcast(X, year + kin + sex_kin + age_kin + age_focal + cohort ~ alive, value.var = "count", fun.aggregate = sum)
  }) %>% data.table::rbindlist()

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- kin
  }
  return(out)
}

#' one time projection kin

#' @description one time projection kin. internal function.
#'
#' @param Ut numeric. A matrix of survival probabilities (or ratios).
#' @param Ft numeric. A matrix of age-specific fertility rates.
#' @param Ft_star numeric. Ft but for female fertility.
#' @param causes integer. Number of causes of death included.
#' @param pit numeric. A matrix with distribution of childbearing.
#' @param sex_focal character. "f" for female or "m" for male.
#' @param ages numeric.
#' @param pkin numeric. A list with kin count distribution in previous year.
#' @return A list of 14 types of kin matrices (kin age by Focal age, blocked for two sex) projected one time interval.
#' @export
timevarying_kin_2sex_cod<- function(Ut, Ft, Ft_star, causes, pit, sex_focal, ages, pkin){

  agess <- ages*2
  om <- ages-1
  pif <- pit[1:ages]
  pim <- pit[(ages+1):agess]

  # BEN : Add the number of CoD - IW: already as argument (Hf is not an argument)
  # causes <- nrow(Hf[[1]])

  # BEN: Changed dimensions of lower part (dead kin) to account for death from causes.
  phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0, (agess + agess*causes), ages)

  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                   nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # G matrix moves focal by age
  G <- matrix(0, nrow=ages, ncol=ages)
  G[row(G)-1 == col(G)] <- 1

  # BEN: Changed dimensions
  Gt <- matrix(0, (agess + agess*causes), (agess + agess*causes))

  Gt[1:(agess), 1:(agess)] <- as.matrix(Matrix::bdiag(G, G))

  # locate focal at age 0 depending sex
  sex_index <- ifelse(sex_focal == "f", 1, ages+1)
  phi[sex_index, 1] <- 1

  # BEN: NOT SURE ABOUT WHAT IS HAPPENING BELOW
  # Rows are multiplied by the sum of the pi?

  # initial distribution
  m[1:agess,1]   = pit
  gm[1:agess,1]  = pkin[["m"]][1:agess,] %*% (pif + pim)
  ggm[1:agess,1] = pkin[["gm"]][1:agess,] %*% (pif + pim)
  oa[1:agess,1]  = pkin[["os"]][1:agess,] %*% (pif + pim)
  ya[1:agess,1]  = pkin[["ys"]][1:agess,] %*% (pif + pim)
  coa[1:agess,1] = pkin[["nos"]][1:agess,] %*% (pif + pim)
  cya[1:agess,1] = pkin[["nys"]][1:agess,] %*% (pif + pim)
  # atribuible to focal sex
  pios <- if(sex_focal == "f") pif else pim
  os[1:agess,1]  = pkin[["d"]][1:agess,] %*% pios
  nos[1:agess,1] = pkin[["gd"]][1:ages,] %*% pios

  for (ix in 1:om){
    phi[,ix+1] = Gt %*% phi[, ix]
    d[,ix+1]   = Ut %*% pkin[["d"]][,ix]   + Ft %*% phi[,ix]
    gd[,ix+1]  = Ut %*% pkin[["gd"]][,ix]  + Ft %*% pkin[["d"]][,ix]
    ggd[,ix+1] = Ut %*% pkin[["ggd"]][,ix] + Ft %*% pkin[["gd"]][,ix]
    m[,ix+1]   = Ut %*% pkin[["m"]][,ix]
    gm[,ix+1]  = Ut %*% pkin[["gm"]][,ix]
    ggm[,ix+1] = Ut %*% pkin[["ggm"]][,ix]
    os[,ix+1]  = Ut %*% pkin[["os"]][,ix]
    ys[,ix+1]  = Ut %*% pkin[["ys"]][,ix]  + Ft_star %*% pkin[["m"]][,ix]
    nos[,ix+1] = Ut %*% pkin[["nos"]][,ix] + Ft %*% pkin[["os"]][,ix]
    nys[,ix+1] = Ut %*% pkin[["nys"]][,ix] + Ft %*% pkin[["ys"]][,ix]
    oa[,ix+1]  = Ut %*% pkin[["oa"]][,ix]
    ya[,ix+1]  = Ut %*% pkin[["ya"]][,ix]  + Ft_star %*% pkin[["gm"]][,ix]
    coa[,ix+1] = Ut %*% pkin[["coa"]][,ix] + Ft %*% pkin[["oa"]][,ix]
    cya[,ix+1] = Ut %*% pkin[["cya"]][,ix] + Ft %*% pkin[["ya"]][,ix]
  }

  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                   nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  return(kin_list)
}

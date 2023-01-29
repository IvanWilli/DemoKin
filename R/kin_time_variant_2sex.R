#' Estimate kin counts in a time variant framework (dynamic rates) in a two-sex framework (Caswell, 2022)

#' @description Two-sex matrix framework for kin count estimates with varying rates.
#' This produces kin counts grouped by kin, age and sex of each relatives at each Focal´s age.
#' For example, male cousins from aunts and uncles from different sibling's parents are grouped in one male count of cousins.
#' @details See Caswell (2022) for details on formulas.
#' @param pf numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param pm numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param ff numeric. Same as pf but for fertility rates.
#' @param fm numeric. Same as pm but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates. Default `TRUE`.
#' @param sex_focal character. "f" for female or "m" for male.
#' @param pif numeric. For using some specific age distribution of childbearing for mothers (same length as ages). Default `NULL`.
#' @param pim numeric. For using some specific age distribution of childbearing for fathers (same length as ages). Default `NULL`.
#' @param nf numeric. Same as pf but for population distribution (counts or `%`). Optional.
#' @param nm numeric. Same as pm but for population distribution (counts or `%`). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth. This multiplies `f` argument. If `f` is already for female offspring, this needs to be set as 1.
#' @param stable logic. Deprecated. Use `time_invariant`.
#' @return A data.frame with year, cohort, Focal´s age, related ages, sex and type of kin (for example `d` is daughter, `oa` is older aunts, etc.), including living and dead kin at that age and sex.
#' @export

kin_time_variant_2sex <- function(pf = NULL, pm = NULL,
                                   ff = NULL, fm = NULL,
                                   sex_focal = "f",
                                   birth_female = 1/2.04,
                                   pif = NULL, pim = NULL,
                                   nf = NULL, nm = NULL,
                                   output_cohort = NULL, output_period = NULL, output_kin = NULL,
                                   list_output = FALSE){

  # same input length
  if(!all(dim(pf) == dim(pm), dim(pf) == dim(ff), dim(pf) == dim(fm))) stop("Dimension of P's and F's should be the same")

  # data should be from same interval years
  years_data <- as.integer(colnames(pf))
  if(var(diff(years_data))!=0) stop("Data should be for same interval length years. Fill the gaps and run again")

  # utils
  age          <- 0:(nrow(pf)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  agess        <- ages*2
  om           <- max(age)
  zeros        <- matrix(0, nrow=ages, ncol=ages)

  # age distribution at childborn
  Pif <- pif
  Pim <- pim
  if(is.null(pif)){
    if(!is.null(nf)){
      Pif <- rbind(t(t(nf * ff)/colSums(nf * ff)), matrix(0,ages,length(years_data)))
    }else{
      Pif <- matrix(0, nrow=ages, ncol=n_years_data)
      no_Pif <- TRUE
    }
  }
  if(is.null(pim)){
    if(!is.null(nm)){
      Pim <- rbind(t(t(nm * fm)/colSums(nm * fm)), matrix(0,ages,length(years_data)))
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
  for(t in 1:n_years_data){
    # t = 1
    Uf = Um = Fft = Fmt = Mm = Mf = Gt = zeros = matrix(0, nrow=ages, ncol=ages)
    Uf[row(Uf)-1 == col(Uf)] <- pf[-ages,t]
    Uf[ages, ages] = pf[ages,t]
    Um[row(Um)-1 == col(Um)] <- pm[-ages,t]
    Um[ages, ages] = pm[ages,t]
    Mm <- diag(1-pm[,t])
    Mf <- diag(1-pf[,t])
    Ut <- as.matrix(rbind(
      cbind(Matrix::bdiag(Uf, Um), Matrix::bdiag(zeros, zeros)),
      cbind(Matrix::bdiag(Mf, Mm), Matrix::bdiag(zeros, zeros))))
    Ul[[as.character(years_data[t])]] <- Ut
    Fft[1,] = ff[,t]
    Fmt[1,] = fm[,t]
    Ft <- Ft_star <- matrix(0, agess*2, agess*2)
    Ft[1:agess,1:agess] <- rbind(cbind(birth_female * Fft, birth_female * Fmt),
                                 cbind((1-birth_female) * Fft, (1-birth_female) * Fmt))
    Ft_star[1:agess,1:ages] <- rbind(birth_female * Fft, (1-birth_female) * Fft)
    Fl[[as.character(years_data[t])]] <- Ft
    Fl_star[[as.character(years_data[t])]] <- Ft_star
    if(no_Pif){
      A <- Uf + Fft
      A_decomp = eigen(A)
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      Pif[,t] <- w*A[1,]/sum(w*A[1,])
    }
    if(no_Pim){
      A <- Um + Fmt
      A_decomp = eigen(A)
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      Pim[,t] <- w*A[1,]/sum(w*A[1,])
    }
    # project
    Ut <- as.matrix(Ul[[t]])
    Ft <- as.matrix(Fl[[t]])
    Ft_star <- as.matrix(Fl_star[[t]])
    pitf <- Pif[,t]
    pitm <- Pim[,t]
    pit <- c(pitf, pitm)
    if (t==1){
      p1f <- pf[,1]
      p1m <- pm[,1]
      f1f <- ff[,1]
      f1m <- fm[,1]
      pif1 <- Pif[,1]
      pim1 <- Pim[,1]
      kin_all[[1]] <- kin_time_invariant_2sex(pf = p1f, pm = p1m,
                                              ff = f1f, fm = f1m,
                                              pif = pif1, pim = pim1,
                                              birth_female = birth_female, list_output = TRUE)
    }
    kin_all[[t+1]] <- timevarying_kin_2sex(Ut=Ut, Ft=Ft, Ft_star=Ft_star, pit=pit, sex_focal, ages, pkin=kin_all[[t]])
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
  kin <- lapply(names(kin_list), FUN = function(Y){
    X <- kin_list[[Y]]
    X <- purrr::map2(X, names(X), function(x,y){
      x <- as.data.frame(x)
      x$year <- Y
      x$kin <- y
      x$sex_kin <- rep(c(rep("f",ages), rep("m",ages)),2)
      x$age_kin <- rep(age,2)
      x$alive <- c(rep("living",ages), rep("dead",ages))
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
  pb$tick()

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
#' @param ft numeric. A matrix of age-specific fertility rates.
#' @param pit numeric. A matrix with distribution of childbearing.
#' @param ages numeric.
#' @param pkin numeric. A list with kin count distribution in previous year.
#
timevarying_kin_2sex<- function(Ut, Ft, Ft_star, pit, sex_focal, ages, pkin){

  agess <- ages*2
  om <- ages-1
  pif <- pit[1:ages]
  pim <- pit[(ages+1):agess]
  phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,agess*2,ages)
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # G matrix moves focal by age
  G <- matrix(0, nrow=ages, ncol=ages)
  G[row(G)-1 == col(G)] <- 1
  Gt <- matrix(0, agess*2, agess*2)
  Gt[1:(agess), 1:(agess)] <- as.matrix(Matrix::bdiag(G, G))

  # locate focal at age 0 depending sex
  sex_index <- ifelse(sex_focal == "f", 1, ages+1)
  phi[sex_index, 1] <- 1

  # initial distribution
  m[1:agess,1]   = pit
  gm[1:agess,1]  = pkin[["m"]][1:agess,] %*% (pif + pim)
  ggm[1:agess,1] = pkin[["gm"]][1:agess,] %*% (pif + pim)
  os[1:agess,1]  = pkin[["d"]][1:agess,] %*% pif
  nos[1:agess,1] = pkin[["gd"]][1:ages,] %*% pif
  oa[1:agess,1]  = pkin[["os"]][1:agess,] %*% (pif + pim)
  ya[1:agess,1]  = pkin[["ys"]][1:agess,] %*% (pif + pim)
  coa[1:agess,1] = pkin[["nos"]][1:agess,] %*% (pif + pim)
  cya[1:agess,1] = pkin[["nys"]][1:agess,] %*% (pif + pim)

  for (ix in 1:om){
    # ix = 1
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

#' defince apc combination to return

#' @description defince apc to return.
#'
output_period_cohort_combination <- function(output_cohort = NULL, output_period = NULL, age = NULL, years_data = NULL){

  # no specific
  if(is.null(output_period) & is.null(output_cohort)){
    message("No specific output was set. Return all period data.")
    output_period <- years_data
  }

  # cohort combination
  if(!is.null(output_cohort)){
    selected_cohorts_year_age <- data.frame(age = rep(age,length(output_cohort)),
                                            year = purrr::map(output_cohort,.f = ~.x+age) %>%
                                              unlist(use.names = F))
  }else{selected_cohorts_year_age <- c()}

  # period year combination
  if(!is.null(output_period)){selected_years_age <- expand.grid(age, output_period) %>% dplyr::rename(age=1,year=2)
  }else{selected_years_age <- c()}

  # end
  return(dplyr::bind_rows(selected_years_age,selected_cohorts_year_age) %>% dplyr::distinct())
}

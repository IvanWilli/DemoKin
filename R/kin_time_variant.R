#' Estimate kin counts in a time variant framework

#' @description Implementation of time variant Goodman-Keyfitz-Pullum equations based on Caswell (2021).
#'
#' @param U numeric. A matrix of survival ratios with rows as ages and columns as years. Column names must be equal interval.
#' @param f numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. Coincident with `U`.
#' @param N numeric. A matrix of population with rows as ages and columns as years. Coincident with `U`.
#' @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. Coincident with `U`.
#' @param output_cohort integer. Year of birth of focal to return as output. Could be a vector. Should be within input data years range.
#' @param output_period integer. Year for which to return kinship structure. Could be a vector. Should be within input data years range.
#' @param output_kin character. kin to return as output: "m" for mother, "d" for daughter,... See `vignette` for exahustive kin.
#' @param birth_female numeric. Female portion at birth.
#' @param list_output logical. Results as a list with years elements (as a result of `output_cohort` and `output_period` combination), with a second list of `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`

#' @return A data frame of population kinship structure, with focal's cohort, focal´s age, period year, type of relatives
#' (for example `d` is daughter, `oa` is older aunts, etc.), living and death kin counts, and age of (living or time deceased) relatives. If `list_output = TRUE` then this is a list.
#' @export

kin_time_variant <- function(U = NULL, f = NULL, N = NULL, pi = NULL,
                            output_cohort = NULL, output_period = NULL, output_kin = NULL,
                            birth_female = 1/2.04, list_output = FALSE){

  # check input
  if(is.null(U) | is.null(f)) stop("You need values on U and/or f.")

  # diff years
  if(!any(as.integer(colnames(U)) == as.integer(colnames(f)))) stop("Data should be from same years.")

  # data should be from same interval years
  years_data <- as.integer(colnames(U))
  if(var(diff(years_data))!=0) stop("Data should be for same interval length years. Fill the gaps and run again")

  # utils
  age          <- 0:(nrow(U)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  om           <- max(age)
  zeros        <- matrix(0, nrow=ages, ncol=ages)

  # age distribution at childborn
  if(is.null(pi)){
    if(is.null(N)){
      # create pi and fill it during the loop
      message("Stable assumption was made for calculating pi on each year because no input data.")
      pi_N_null_flag <- TRUE
      pi <- matrix(0, nrow=ages, ncol=n_years_data)
    }else{
      pi_N_null_flag <-FALSE
      pi <- rbind(t(t(N * f)/colSums(N * f)), matrix(0,ages,length(years_data)))
    }
  }

  # get lists of matrix
  Ul = fl = list()
  for(t in 1:n_years_data){
    Ut = Mt = Dcum = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- U[-ages,t]
    Ut[ages, ages]=U[ages,t]
    diag(Mt) = 1 - U[,t]
    Ul[[as.character(years_data[t])]] <- rbind(cbind(Ut,zeros),cbind(Mt,Dcum))
    ft = matrix(0, nrow=ages*2, ncol=ages*2)
    ft[1,1:ages] = f[,t] * birth_female
    fl[[as.character(years_data[t])]] <- ft
  }
  U <- Ul
  f <- fl

  # loop over years (more performance here)
  kin_all <- list()
  pb <- progress::progress_bar$new(
    format = "Running over input years [:bar] :percent",
    total = n_years_data, clear = FALSE, width = 60)
  for (iyear in 1:n_years_data){
    # print(iyear)
    Ut <- as.matrix(U[[iyear]])
    ft <- as.matrix(f[[iyear]])
    if(pi_N_null_flag){
      A <- Ut[1:ages,1:ages] + ft[1:ages,1:ages]
      A_decomp = eigen(A)
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      pit <- pi[,iyear] <- w*A[1,]/sum(w*A[1,])
    }else{
      pit <- pi[,iyear]
    }
    if (iyear==1){
      U1 <- c(diag(Ut[-1,])[1:om],Ut[om,om])
      f1 <- ft[1,][1:ages]
      pi1 <- pit[1:ages]
      kin_all[[1]] <- kin_time_invariant(U = U1, f = f1/birth_female, pi = pi1, birth_female = birth_female,
                                         list_output = TRUE)
    }
    kin_all[[iyear+1]] <- timevarying_kin(Ut=Ut,ft=ft,pit=pit,ages,pkin=kin_all[[iyear]])
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
  kin <- lapply(names(kin_list), function(Y){
    X <- kin_list[[Y]]
    X <- purrr::map2(X, names(X), function(x,y) as.data.frame(x) %>%
                dplyr::mutate(year = Y,
                        kin=y,
                        age_kin = rep(age,2),
                        alive = c(rep("living",ages), rep("dead",ages)),
                        .before=everything())) %>%
      dplyr::bind_rows() %>%
      stats::setNames(c("year","kin","age_kin","alive",as.character(age))) %>%
      tidyr::gather(age_focal, count,-age_kin, -kin, -year, -alive) %>%
      dplyr::mutate(age_focal = as.integer(age_focal),
                     year = as.integer(year),
                     cohort = year - age_focal) %>%
      dplyr::filter(age_focal %in% out_selected$age[out_selected$year==as.integer(Y)]) %>%
      tidyr::pivot_wider(names_from = alive, values_from = count)}) %>%
    dplyr::bind_rows()

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
timevarying_kin<- function(Ut, ft, pit, ages, pkin){

  # frequently used zero vector for initial condition
  zvec=rep(0,ages*2);
  I = matrix(0, ages * 2, ages * 2)
  diag(I[1:ages,1:ages]) = 1
  om=ages-1;
  d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,ages*2,ages)
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # initial distribution
  d[,1]=gd[,1]=ggd[,1]=ys[,1]=nys[,1]=zvec
  m[1:ages,1]  = pit[1:ages]
  gm[1:ages,1] = pkin[["m"]][1:ages,] %*% pit[1:ages]
  ggm[1:ages,1]= pkin[["gm"]][1:ages,] %*% pit[1:ages]
  os[1:ages,1] = pkin[["d"]][1:ages,] %*% pit[1:ages]
  nos[1:ages,1] = pkin[["gd"]][1:ages,] %*% pit[1:ages]
  oa[1:ages,1] = pkin[["os"]][1:ages,] %*% pit[1:ages]
  ya[1:ages,1] = pkin[["ys"]][1:ages,] %*% pit[1:ages]
  coa[1:ages,1]= pkin[["nos"]][1:ages,] %*% pit[1:ages]
  cya[1:ages,1]= pkin[["nys"]][1:ages,] %*% pit[1:ages]

  for (ix in 1:om){
    d[,ix+1]   = Ut %*% pkin[["d"]][,ix]   + ft %*% I[,ix]
    gd[,ix+1]  = Ut %*% pkin[["gd"]][,ix]  + ft %*% pkin[["d"]][,ix]
    ggd[,ix+1] = Ut %*% pkin[["ggd"]][,ix] + ft %*% pkin[["gd"]][,ix]
    m[,ix+1]   = Ut %*% pkin[["m"]][,ix]
    gm[,ix+1]  = Ut %*% pkin[["gm"]][,ix]
    ggm[,ix+1] = Ut %*% pkin[["ggm"]][,ix]
    os[,ix+1]  = Ut %*% pkin[["os"]][,ix]
    ys[,ix+1]  = Ut %*% pkin[["ys"]][,ix]  + ft %*% pkin[["m"]][,ix]
    nos[,ix+1] = Ut %*% pkin[["nos"]][,ix] + ft %*% pkin[["os"]][,ix]
    nys[,ix+1] = Ut %*% pkin[["nys"]][,ix] + ft %*% pkin[["ys"]][,ix]
    oa[,ix+1]  = Ut %*% pkin[["oa"]][,ix]
    ya[,ix+1]  = Ut %*% pkin[["ya"]][,ix]  + ft %*% pkin[["gm"]][,ix]
    coa[,ix+1] = Ut %*% pkin[["coa"]][,ix] + ft %*% pkin[["oa"]][,ix]
    cya[,ix+1] = Ut %*% pkin[["cya"]][,ix] + ft %*% pkin[["ya"]][,ix]
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

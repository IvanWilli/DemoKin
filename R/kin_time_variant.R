#' Estimate kin counts in a time variant framework (dynamic rates) for one-sex model (matrilineal/patrilineal)

#' @description Matrix implementation of time variant Goodman-Keyfitz-Pullum equations in a matrix framework.
#' @details See Caswell (2021) for details on formulas.
#' @param p numeric. A matrix of survival ratios with rows as ages and columns as years. Column names must be equal interval.
#' @param f numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. Coincident with `U`.
#' @param n numeric. A matrix of population with rows as ages and columns as years. Coincident with `U`.
#' @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. Coincident with `U`.
#' @param output_cohort integer. Year of birth of focal to return as output. Could be a vector. Should be within input data years range.
#' @param output_period integer. Year for which to return kinship structure. Could be a vector. Should be within input data years range.
#' @param output_kin character. kin to return as output: "m" for mother, "d" for daughter,... See `vignette` for exahustive kin.
#' @param birth_female numeric. Female portion at birth.
#' @param list_output logical. Results as a list with years elements (as a result of `output_cohort` and `output_period` combination), with a second list of `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`

#' @return A data frame of population kinship structure, with Focal's cohort, focal´s age, period year, type of relatives
#' (for example `d` is daughter, `oa` is older aunts, etc.), living and death kin counts, and age of (living or time deceased) relatives. If `list_output = TRUE` then this is a list.
#' @export

kin_time_variant <- function(p = NULL, f = NULL, pi = NULL, n = NULL,
                            output_cohort = NULL, output_period = NULL, output_kin = NULL,
                            birth_female = 1/2.04, list_output = FALSE){

  # global vars
  .<-living<-dead<-age_kin<-age_focal<-cohort<-year<-total<-mean_age<-count_living<-sd_age<-count_dead<-mean_age_lost<-indicator<-value<-NULL

  # check input
  if(is.null(p) | is.null(f)) stop("You need values on p and f.")

  # diff years
  if(!any(as.integer(colnames(p)) == as.integer(colnames(f)))) stop("Make sure that p and f are matrices and have the same column names.")

  # data should be from same interval years
  years_data <- as.integer(colnames(p))
  if(stats::var(diff(years_data))!=0) stop("The years given as column names in the p and f matrices must be equally spaced.")

  # utils
  age          <- 0:(nrow(p)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  om           <- max(age)
  zeros        <- matrix(0, nrow=ages, ncol=ages)

  # consider input data for age distribution at child born, or flag it
  no_Pi <- FALSE
  if(is.null(pi)){
    if(is.null(n)){
      # create pi and fill it during the loop
      no_Pi <- TRUE
      pi <- matrix(0, nrow=ages, ncol=n_years_data)
    }else{
      no_Pi <- FALSE
      pi <- rbind(t(t(n * f)/colSums(n * f)), matrix(0,ages,length(years_data)))
    }
  }

  # loop over years
  kin_all <- list()
  pb <- progress::progress_bar$new(
    format = "Running over input years [:bar] :percent",
    total = n_years_data + 1, clear = FALSE, width = 50)
  for (t in 1:n_years_data){
    # build set of matrix
    Ut <- Mt <- matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- p[-ages,t]
    Ut[ages, ages] <- p[ages,t]
    diag(Mt) <- 1 - p[,t]
    Ut <- rbind(cbind(Ut,zeros),cbind(Mt,zeros))
    ft <- matrix(0, nrow=ages*2, ncol=ages*2)
    ft[1,1:ages] <- f[,t] * birth_female
    A <- Ut[1:ages,1:ages] + ft[1:ages,1:ages]
    # stable assumption at start
    if (t==1){
      p1 <- c(diag(Ut[-1,])[1:om],Ut[om,om])
      f1 <- ft[1,][1:ages]/birth_female
      A_decomp <- eigen(A)
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      pit <- w*A[1,]/sum(w*A[1,])
      pi1 <- pit[1:ages]
      kin_all[[1]] <- kin_time_invariant(p = p1, f = f1, pi = pi1, birth_female = birth_female,
                                         list_output = TRUE)
    }
    # project pi
    if(no_Pi){
      w <- A %*% w
      pi[,t] <- w*A[1,]/sum(w*A[1,])
    }
    # kin for next year
    kin_all[[t+1]] <- timevarying_kin(Ut = Ut, ft = ft, pit = pi[,t], ages, pkin = kin_all[[t]])
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
      x[(ages+1):(ages*2),1:(ages-1)] <- x[(ages+1):(ages*2),2:ages]
      x[(ages+1):(ages*2),ages] <- 0
      x <- as.data.frame(x)
      x$year <- Y
      x$kin <- y
      x$age_kin <- rep(age,2)
      x$alive <- c(rep("living",ages), rep("dead",ages))
      return(x)
      }) %>%
      data.table::rbindlist() %>%
      stats::setNames(c(as.character(age), "year","kin","age_kin","alive")) %>%
      data.table::melt(id.vars = c("year","kin","age_kin","alive"), variable.name = "age_focal", value.name = "count")
    X$age_focal = as.integer(as.character(X$age_focal))
    X$year = as.integer(X$year)
    X$cohort = X$year - X$age_focal
    X[X$age_focal %in% out_selected$age[out_selected$year==as.integer(Y)],] %>%
      data.table::dcast(year + kin + age_kin + age_focal + cohort ~ alive, value.var = "count")
    }) %>% data.table::rbindlist()

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- kin
  }

  # end
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
#' @return A list of 14 types of kin matrices (kin age by Focal age) projected one time interval.
#' @export
timevarying_kin<- function(Ut, ft, pit, ages, pkin){

  # frequently used zero vector for initial condition
  zvec <- rep(0,ages*2);
  I <- matrix(0, ages * 2, ages * 2)
  diag(I[1:ages,1:ages]) <- 1
  om <- ages-1;
  d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya <- matrix(0,ages*2,ages)
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # initial distribution
  d[,1] = gd[,1] = ggd[,1] = ys[,1] = nys[,1] = zvec
  m[1:ages,1]   = pit[1:ages]
  gm[1:ages,1]  = pkin[["m"]][1:ages,] %*% pit[1:ages]
  ggm[1:ages,1] = pkin[["gm"]][1:ages,] %*% pit[1:ages]
  os[1:ages,1]  = pkin[["d"]][1:ages,] %*% pit[1:ages]
  nos[1:ages,1] = pkin[["gd"]][1:ages,] %*% pit[1:ages]
  oa[1:ages,1]  = pkin[["os"]][1:ages,] %*% pit[1:ages]
  ya[1:ages,1]  = pkin[["ys"]][1:ages,] %*% pit[1:ages]
  coa[1:ages,1] = pkin[["nos"]][1:ages,] %*% pit[1:ages]
  cya[1:ages,1] = pkin[["nys"]][1:ages,] %*% pit[1:ages]

  # focal´s trip
  for(ix in 1:om){
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

#' APC combination to return

#' @description define APC combination to return in `kin` and `kin2sex`.
#' @details Because returning all period and cohort data from a huge time-series would be hard memory consuming,
#' this function is an auxiliary one to deal with selection from inputs `output_cohort` and `output_period`.
#' @param output_cohort integer. A vector with selected calendar years.
#' @param output_period integer. A vector with selected cohort years.
#' @param age integer. A vector with ages from the kinship network to be filtered.
#' @param years_data integer. A vector with years from the time-varying kinship network to be filtered.
#' @return data.frame with years and ages to filter in `kin` and `kin_2sex` functions.
#' @export
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

  # period combination
  if(!is.null(output_period)){selected_years_age <- expand.grid(age, output_period) %>% dplyr::rename(age=1,year=2)
  }else{selected_years_age <- c()}

  # end
  return(dplyr::bind_rows(selected_years_age,selected_cohorts_year_age) %>% dplyr::distinct())
}


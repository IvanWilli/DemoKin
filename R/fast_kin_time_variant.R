# Fast (C++) implementation of kin_time_variant(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_time_variant.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).

# Estimate kin counts in a time variant framework (dynamic rates) for one-sex model (matrilineal/patrilineal)

# @description Matrix implementation of time variant Goodman-Keyfitz-Pullum equations in a matrix framework.
# @details See Caswell (2021) for details on formulas.
# @param p numeric. A matrix of survival ratios with rows as ages and columns as years. Column names must be equal interval.
# @param f numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. Coincident with `U`.
# @param n numeric. A matrix of population with rows as ages and columns as years. Coincident with `U`.
# @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. Coincident with `U`.
# @param output_cohort integer. Year of birth of focal to return as output. Could be a vector. Should be within input data years range.
# @param output_period integer. Year for which to return kinship structure. Could be a vector. Should be within input data years range.
# @param output_kin character. kin to return as output: "m" for mother, "d" for daughter,... See `vignette` for exahustive kin.
# @param birth_female numeric. Female portion at birth. This can be a vector of length equal to the number of years in the input data, or a single value that will be repeated for all years.
# @param list_output logical. Results as a list with years elements (as a result of `output_cohort` and `output_period` combination), with a second list of `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`

# @return A data frame of population kinship structure, with Focal's cohort, focal´s age, period year, type of relatives
# (for example `d` is daughter, `oa` is older aunts, etc.), living and death kin counts, and age of (living or time deceased) relatives. If `list_output = TRUE` then this is a list.

kin_time_variant_fast <- function(p = NULL, f = NULL, pi = NULL, n = NULL,
                            output_cohort = NULL, output_period = NULL, output_kin = NULL,
                            birth_female = 1/2.04, list_output = FALSE){

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
  year_step    <- diff(years_data)[1]

  # birth_female needs to be dynamic. Complete in case length is lower than data
  if(length(birth_female) < n_years_data){
    last_birth_female <- utils::tail(birth_female, n=1)
    n_birth_female <- length(birth_female)
    birth_female <- c(birth_female, rep(last_birth_female, n_years_data - n_birth_female))
    message("Length of birth_female lower than years of risk. Completed with last value...")
  }

  # consider input data for age distribution at child born, or flag it
  no_Pi <- FALSE
  if(is.null(pi)){
    if(is.null(n)){
      # create pi and fill it during the loop
      no_Pi <- TRUE
      pi <- matrix(0, nrow=ages, ncol=n_years_data)
    }else{
      no_Pi <- FALSE
      pi <- rbind(t(t(n * f)/colSums(n * f)), matrix(0,ages,n_years_data))
    }
  }

  # years and ages to return: decided up front so that only those are kept in memory
  out_selected <- suppressMessages(
    output_period_cohort_combination(output_cohort, output_period, age = age, years_data = years_data))
  ages_by_year <- selected_ages_by_year(out_selected)
  possible_kin <- kin_names_14
  if(is.null(output_kin)){
    selected_kin_position <- 1:length(possible_kin)
  }else{
    selected_kin_position <- which(possible_kin %in% output_kin)
  }
  years_all <- as.character(c(years_data, dplyr::last(years_data) + year_step))
  store <- list()
  keep_year <- function(pkin, Y){
    if(!Y %in% names(ages_by_year)) return(invisible(NULL))
    pkin <- pkin[selected_kin_position]
    if(list_output){
      store[[Y]] <<- pkin
    }else{
      A <- ages_by_year[[Y]]
      store[[Y]] <<- list(ages = A, kin = lapply(pkin, function(x) x[, A + 1, drop = FALSE]))
    }
  }

  # loop over years
  message("Running kinship model...")
  pkin <- NULL
  for (t in 1:n_years_data){
    # build set of matrix
    Ut <- Mt <- matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- p[-ages,t]
    Ut[ages, ages] <- p[ages,t]
    diag(Mt) <- 1 - p[,t]
    Ut <- rbind(cbind(Ut,zeros),cbind(Mt,zeros))
    ft <- matrix(0, nrow=ages*2, ncol=ages*2)
    ft[1,1:ages] <- f[,t] * birth_female[t]
    A <- Ut[1:ages,1:ages] + ft[1:ages,1:ages]
    # stable assumption at start
    if (t==1){
      p1 <- c(diag(Ut[-1,])[1:om],Ut[ages,ages])
      f1 <- ft[1,][1:ages]/birth_female[1]
      A_decomp <- eigen(A)
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      pit <- w*A[1,]/sum(w*A[1,])
      pi1 <- pit[1:ages]
      pkin <- kin_time_invariant_fast(p = p1, f = f1, pi = pi1, birth_female = birth_female[1],
                                 list_output = TRUE)
      keep_year(pkin, years_all[1])
    }
    # project pi
    if(no_Pi){
      w <- A %*% w
      pi[,t] <- w*A[1,]/sum(w*A[1,])
    }
    # kin for next year
    pkin <- project_kin_1sex(Ut, ft, pi[,t], ages, prev = pkin)
    keep_year(pkin, years_all[t+1])
  }

  message("Preparing output...")
  if(is.null(output_period) & is.null(output_cohort)) message("No specific output was set. Return all period data.")

  # results as list?
  if(list_output) {
    out <- store
  }else{
    # build only the output years (identical to building all and subsetting)
    keep <- if (!is.null(output_period) && any(as.integer(names(store)) %in% output_period)) output_period
    kin <- build_kin_df_tv(store, year_step, layout_1sex(age), keep_years = keep)
    if (!is.null(output_period) && is.null(keep)) {
      kin <- kin[kin$year %in% output_period, ]
    }
    out <- kin
  }

  # end
  return(out)
}

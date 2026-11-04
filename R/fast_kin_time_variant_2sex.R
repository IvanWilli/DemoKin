# Fast (C++) implementation of kin_time_variant_2sex(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_time_variant_2sex.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).

# Estimate kin counts in a time variant framework (dynamic rates) in a two-sex framework (Caswell, 2022)

# @description Two-sex matrix framework for kin count estimates with varying rates.
# This produces kin counts grouped by kin, age and sex of each relatives at each Focal´s age.
# For example, male cousins from aunts and uncles from different sibling's parents are grouped in one male count of cousins.
# @details See Caswell (2022) for details on formulas.
# @param pf numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
# @param pm numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
# @param ff numeric. Same as pf but for fertility rates.
# @param fm numeric. Same as pm but for fertility rates.
# @param sex_focal character. "f" for female or "m" for male.
# @param pif numeric. For using some specific age distribution of childbearing for mothers (same length as ages). Default `NULL`.
# @param pim numeric. For using some specific age distribution of childbearing for fathers (same length as ages). Default `NULL`.
# @param nf numeric. Same as pf but for population distribution (counts or `%`). Optional.
# @param nm numeric. Same as pm but for population distribution (counts or `%`). Optional.
# @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
# @param output_period integer. Vector of period years for returning results. Should be within input data years range.
# @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
# @param birth_female numeric. Female portion at birth. This can be a vector of length equal to the number of years in the input data, or a single value that will be repeated for all years.
# @param list_output logical. Results as a list with years elements (as a result of `output_cohort` and `output_period` combination), with a second list of `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
# @return A data.frame with year, cohort, Focal´s age, related ages, sex and type of kin (for example `d` is daughter, `oa` is older aunts, etc.), including living and dead kin at that age and sex.

kin_time_variant_2sex_fast <- function(pf = NULL, pm = NULL,
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
  if(stats::var(diff(years_data))!=0) stop("Data should be for same interval length years. Fill the gaps and run again")

  # utils
  age          <- 0:(nrow(pf)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  agess        <- ages*2
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

  # consider input data for age distribution at child born, or flag it to fill it
  Pif <- pif; no_Pif <- FALSE
  if(is.null(pif)){
    if(!is.null(nf)){
      Pif <- t(t(nf * ff)/colSums(nf * ff))
    }else{
      Pif <- matrix(0, nrow=ages, ncol=n_years_data)
      no_Pif <- TRUE
    }
  }
  Pim <- pim; no_Pim <- FALSE
  if(is.null(pim)){
    if(!is.null(nm)){
      Pim <- t(t(nm * fm)/colSums(nm * fm))
    }else{
      Pim <- matrix(0, nrow=ages, ncol=n_years_data)
      no_Pim <- TRUE
    }
  }

  # years, ages and kin to return: decided up front so that only those are kept
  out_selected <- suppressMessages(
    output_period_cohort_combination(output_cohort, output_period, age = age, years_data = years_data))
  ages_by_year <- selected_ages_by_year(out_selected)
  possible_kin <- kin_names_14
  if(is.null(output_kin)){
    selected_kin_position <- 1:length(possible_kin)
  }else{
    selected_kin_position <- which(possible_kin %in% output_kin)
  }
  years_all <- as.character(c(years_data, dplyr::last(years_data) + dplyr::last(diff(years_data))))
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

  # subdiagonal (ageing) positions of an ages x ages block, in column order
  sub_i <- 2:ages; sub_j <- 1:(ages-1); age1 <- seq_len(ages)

  # loop over years
  message("Running kinship model...")
  pkin <- NULL
  for(t in 1:n_years_data){
    # The projection matrices go to the C++ kernel as triplet lists of their
    # non-zero entries, built from the rates with the same arithmetic as the
    # dense matrices they replace (Ut: survival on the subdiagonal and in the
    # last age, deaths 1 - p below; Ft, Ft_star: birth_female * f in the rows
    # of newborn females and males). The kernel keeps exactly the non-zero
    # entries either way; the 404 x 404 matrices are no longer allocated.
    bf <- birth_female[t]
    Ut <- list(i = c(sub_i, ages, ages + sub_i, agess, agess + age1, agess + ages + age1),
               j = c(sub_j, ages, ages + sub_j, agess, age1, ages + age1),
               x = c(pf[-ages,t], pf[ages,t], pm[-ages,t], pm[ages,t], 1-pf[,t], 1-pm[,t]),
               nrow = 2L*agess, ncol = 2L*agess)
    Ft <- list(i = rep(c(1L, ages + 1L), each = agess), j = rep(seq_len(agess), 2),
               x = c(bf * ff[,t], bf * fm[,t], (1-bf) * ff[,t], (1-bf) * fm[,t]),
               nrow = 2L*agess, ncol = 2L*agess)
    Fs_f <- bf * ff[,t]; Fs_m <- (1-bf) * ff[,t]
    Ft_star <- list(i = rep(c(1L, ages + 1L), each = ages), j = rep(age1, 2),
                    x = c(Fs_f, Fs_m), nrow = 2L*agess, ncol = 2L*agess)
    # A = bdiag(Uf, Um) + Ft_star[1:agess, 1:agess], dense, as before
    A <- Fs <- matrix(0, agess, agess)
    A[cbind(sub_i, sub_j)] <- pf[-ages,t]
    A[ages, ages] <- pf[ages,t]
    A[cbind(ages + sub_i, ages + sub_j)] <- pm[-ages,t]
    A[agess, agess] <- pm[ages,t]
    Fs[1, age1] <- Fs_f
    Fs[ages + 1, age1] <- Fs_m
    A <- A + Fs

    # stable assumption at start
    if (t==1){
      p1f <- pf[,1]; p1m <- pm[,1]
      f1f <- ff[,1]; f1m <- fm[,1]
      # time boundary for pi
      A_decomp <- eigen(A)
      lambda <- as.double(A_decomp$values[1])
      w <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
      wf <- w[1:ages]
      wm <- w[(ages+1):(2*ages)]
      pif1 <- wf * ff[,t] / sum(wf * ff[,t])
      pim1 <- wm * fm[,t] / sum(wm * fm[,t])
      pkin <- kin_time_invariant_2sex_fast(pf = p1f, pm = p1m,
                                      ff = f1f, fm = f1m,
                                      sex_focal = sex_focal,
                                      pif = pif1, pim = pim1,
                                      birth_female = birth_female[1], list_output = TRUE)
      keep_year(pkin, years_all[1])
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
    pkin <- project_kin_2sex(Ut, Ft, Ft_star, pit, sex_focal, ages, prev = pkin)
    keep_year(pkin, years_all[t+1])
  }

  message(" Preparing output...")
  if(is.null(output_period) & is.null(output_cohort)) message("No specific output was set. Return all period data.")

  # results as list?
  if(list_output) {
    out <- store
  }else{
    # build only the output years (identical to building all and subsetting)
    keep <- if (!is.null(output_period) && any(as.integer(names(store)) %in% output_period)) output_period
    kin <- build_kin_df_tv(store, year_step, layout_2sex_tv(age), keep_years = keep)
    if (!is.null(output_period) && is.null(keep)) {
      kin <- kin[kin$year %in% output_period, ]
    }
    out <- kin
  }
  return(out)
}

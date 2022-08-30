#' Estimate kin counts

#' @description Implementation of Goodman-Keyfitz-Pullum equations in a matrix framework.
#' @details See Caswell (2019) and Caswell (2021) for details on formulas.
#' @param U numeric. A vector or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param f numeric. Same as U but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates. Default `TRUE`.
#' @param N numeric. Same as U but for population distribution (counts or `%`). Optional.
#' @param pi numeric. Same as U but for childbearing distribution (sum to 1). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth.
#' @return A list with:
#'  * `kin_full`: a data frame with focal´s age, related ages and type of kin (for example `d` is daughter, `oa` is older aunts, etc.), with living kin and death on that age.
#'  * `kin_summary`: a data frame with focal´s age, related ages and type of kin, with indicators obtained processing `kin_full`:
#'  - `count`: count of living kin at actual age of focal
#'  - `mean_age`: mean age of each type of living kin.
#'  - `sd_age`: standard deviation age of each type of living kin .
#'  - `count_death`: count of death kin at specific age of focal.
#'  - `count_cum_death`: cumulated count of death kin at specific age of focal.
#'  - `mean_age_lost`: mean age where focal lost her relative.
#' @export
#'
# get kin ----------------------------------------------------------------
kin <- function(U = NULL, f = NULL,
                 time_invariant = TRUE,
                 N = NULL, pi = NULL,
                 output_cohort = NULL, output_period = NULL, output_kin=NULL,
                 birth_female = 1/2.04)
  {

  age <- as.integer(rownames(U))
  years_data <- as.integer(colnames(U))

  # kin to return
  all_possible_kin <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
  if(is.null(output_kin)){
    output_kin <- all_possible_kin
  }else{
    output_kin <- match.arg(tolower(output_kin), all_possible_kin, several.ok = TRUE)
  }

  # if time dependent or not
  if(time_invariant){
      if(!is.vector(U)) {
        output_period <- min(years_data)
        U <- U[,as.character(output_period)]
        f <- f[,as.character(output_period)]
      }
      kin_full <- kin_time_invariant(U = U, f = f,
                                     output_kin = output_kin, birth_female = birth_female) %>%
                              dplyr::mutate(cohort = NA, year = NA)
  }else{
      if(!is.null(output_cohort) & !is.null(output_period)) stop("sorry, you can not select cohort and period. Choose one please")
      kin_full <- kin_time_variant(U = U, f = f, N = N, pi = pi,
                              output_cohort = output_cohort, output_period = output_period,
                              output_kin = output_kin,
                              birth_female = birth_female)
      message(paste0("Assuming stable population before ", min(years_data), "."))
  }

  # reorder
  kin_full <- kin_full %>% dplyr::select(year, cohort, age_focal, kin, age_kin, living, death)

  # summary
  kin_summary <- dplyr::bind_rows(
    kin_full %>%
      dplyr::rename(total=living) %>%
      dplyr::group_by(cohort, year, age_focal, kin) %>%
      dplyr::summarise(count = sum(total),
                mean_age = sum(total*age_kin)/sum(total),
                sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
      tidyr::pivot_longer(count:sd_age, names_to = "indicator", "value"),
    kin_full %>%
      dplyr::rename(total=death) %>%
      dplyr::group_by(cohort, year, age_focal, kin) %>%
      dplyr::summarise(count_death = sum(total)) %>%
      dplyr::group_by(cohort, year, kin) %>%
      dplyr::mutate(count_cum_death = cumsum(count_death),
                    mean_age_lost = cumsum(count_death * age_focal)/cumsum(count_death)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(count_death:mean_age_lost, names_to = "indicator", "value")) %>%
      dplyr::ungroup() %>%
      dplyr::select(cohort, age_focal, kin, indicator, value) %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value)

    # return
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)
  return(kin_out)
}

#' Estimate kin counts in a two-sex framework

#' @description Implementation of Goodman-Keyfitz-Pullum equations in a matrix framework.
#' @details See Caswell (2022) for details on formulas.
#' @param U numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param f numeric. Same as U but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates. Default `TRUE`.
#' @param N numeric. Same as U but for population distribution (counts or `%`). Optional.
#' @param pi numeric. Same as U but for childbearing distribution (sum to 1). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth. This multiplies `f` argument. If `f` is already for female offspring, this needs to be set as 1.
#' @param stable logic. Deprecated. Use `time_invariant`.
#' @return A list with:
#' \itemize{
#'  \item{kin_full}{ a data frame with year, cohort, Focal´s age, related ages and type of kin (for example `d` is daughter, `oa` is older aunts, etc.), including living and dead kin at that age.}
#'  \item{kin_summary}{ a data frame with Focal´s age, related ages and type of kin, with indicators obtained processing `kin_full`, grouping by cohort or period (depending on the given arguments):}
#'  {\itemize{
#'  \item{`count_living`}{: count of living kin at actual age of Focal}
#'  \item{`mean_age`}{: mean age of each type of living kin.}
#'  \item{`sd_age`}{: standard deviation of age of each type of living kin.}
#'  \item{`count_death`}{: count of dead kin at specific age of Focal.}
#'  \item{`count_cum_death`}{: cumulated count of dead kin until specific age of Focal.}
#'  \item{`mean_age_lost`}{: mean age where Focal lost her relative.}
#'  }
#'  }
#' }

#' @export
#'
# get kin ----------------------------------------------------------------
kin2sex <- function(pf = NULL, pm = NULL, ff = NULL, fm = NULL,
                 time_invariant = TRUE,
                 sex_focal = "f",
                 birth_female = 1/2.04,
                 pif = NULL, pim = NULL,
                 output_cohort = NULL, output_period = NULL, output_kin=NULL)
  {

  age <- as.integer(rownames(pf))
  years_data <- as.integer(colnames(pf))

  # kin to return
  all_possible_kin <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
  if(is.null(output_kin)){
    output_kin <- all_possible_kin
  }else{
    output_kin <- match.arg(tolower(output_kin), all_possible_kin, several.ok = TRUE)
  }

  # if time dependent or not
  if(time_invariant){
      if(!is.vector(pf)) {
        output_period <- min(years_data)
        pf <- pf[,as.character(output_period)]
        pm <- pm[,as.character(output_period)]
        ff <- ff[,as.character(output_period)]
        fm <- fm[,as.character(output_period)]
      }
      kin_full <- kin_time_invariant_2sex(pf, pm, ff, fm,
                                     sex_focal = sex_focal,
                                     birth_female = birth_female,
                                     pif = pif, pim = pim,
                                     output_kin = output_kin) %>%
                              dplyr::mutate(cohort = NA, year = NA)
  }else{
      if(!is.null(output_cohort) & !is.null(output_period)) stop("sorry, you can not select cohort and period. Choose one please")
      kin_full <- kin_time_variant_2sex(Pf = pf, Pm = pm,
                                   Ff = ff, Fm = fm,
                                   sex_focal = sex_focal,
                                   birth_female = birth_female,
                                   Pif = pif, Pim = pim,
                                   output_cohort = output_cohort, output_period = output_period,
                                   output_kin = output_kin)
      message(paste0("Assuming stable population before ", min(years_data), "."))
  }

  # reorder
  kin_full <- kin_full %>%
    dplyr::select(year, cohort, age_focal, sex, kin, age_kin, living, dead) %>%
    dplyr::mutate(kin_group =  dplyr::case_when(kin %in% c("ys", "os") ~ "s",
                                                kin %in% c("ya", "oa") ~ "a",
                                                kin %in% c("coa", "cya") ~ "c",
                                                kin %in% c("nys", "nos") ~ "n",
                                                T ~ kin))

  # summary
  # select period/cohort
  if(!is.null(output_cohort)){
    agrupar <- "cohort"
  } else if(!is.null(output_period)){
    agrupar <- "year"
  } else{
    agrupar <- c("year", "cohort")
  }
  agrupar_no_age_focal <- c("kin", "sex", agrupar)
  agrupar <- c("age_focal", "kin", "sex", agrupar)

  kin_summary <- dplyr::bind_rows(
    kin_full %>%
      dplyr::rename(total=living) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
      dplyr::summarise(count_living = sum(total),
                mean_age = sum(total*age_kin)/sum(total),
                sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
      tidyr::pivot_longer(count_living:sd_age, names_to = "indicator", "value"),
    kin_full %>%
      dplyr::rename(total=dead) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
      dplyr::summarise(count_dead = sum(total)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agrupar_no_age_focal))) %>%
      dplyr::mutate(count_cum_dead = cumsum(count_dead),
                    mean_age_lost = cumsum(count_dead * age_focal)/cumsum(count_dead)) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_longer(count_dead:mean_age_lost, names_to = "indicator", "value")) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value)

    # return
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)
  return(kin_out)
}

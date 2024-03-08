#' Estimate kin counts in a one-sex framework.

#' @description Implementation of Goodman-Keyfitz-Pullum equations in a matrix framework. This produce a matrilineal (or patrilineal)
#' kin count distribution by kin and age.
#' @details See Caswell (2019) and Caswell (2021) for details on formulas. One sex only (female by default).
#' @param p numeric. A vector (atomic) or  matrix with probabilities (or survival ratios, or transition between age class
#' in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param f numeric. Same as `p` but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates. Default `TRUE`.
#' @param n numeric. Only for `time_invariant = FALSE`. Same as `p` but for population distribution (counts or `%`). Optional.
#' @param pi numeric. Same as `U` but for childbearing distribution (sum to 1). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param output_age_focal integer. Vector of ages to select (and make faster the run).
#' @param birth_female numeric. Female portion at birth. This multiplies `f` argument. If `f` is already for female offspring,
#' @param summary_kin logical. Whether or not include `kin_summary` table (see output details). Default `TRUE`.
#' this needs to be set as 1.
#' @return A list with:
#' \itemize{
#'  \item{kin_full}{ a data frame with year, cohort, Focal´s age, related ages and type of kin (for example `d` is daughter,
#'  `oa` is older aunts, etc.), including living and dead kin at that age.}
#'  \item{kin_summary}{ a data frame with Focal´s age, related ages and type of kin, with indicators obtained processing `kin_full`,
#'  grouping by cohort or period (depending on the given arguments):}
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
#' @examples
#' # Kin expected matrilineal count for a Swedish female based on 2015 rates.
#' swe_surv_2015 <- swe_px[,"2015"]
#' swe_asfr_2015 <- swe_asfr[,"2015"]
#' # Run kinship models
#' swe_2015 <- kin(p = swe_surv_2015, f = swe_asfr_2015)
#' head(swe_2015$kin_summary)

kin <- function(p = NULL, f = NULL,
                time_invariant = TRUE,
                pi = NULL, n = NULL,
                output_cohort = NULL, output_period = NULL, output_kin=NULL, output_age_focal = NULL,
                birth_female = 1/2.04,
                summary_kin = TRUE)
  {

  # global vars
  living<-dead<-age_kin<-age_focal<-cohort<-year<-total<-mean_age<-count_living<-sd_age<-count_dead<-mean_age_lost<-indicator<-value<-NULL

  # kin to return
  all_possible_kin <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
  output_kin_asked <- output_kin
  if(is.null(output_kin)){
    output_kin <- all_possible_kin
  }else{
    if("s" %in% output_kin) output_kin <- c(output_kin, "os", "ys")
    if("c" %in% output_kin) output_kin <- c(output_kin, "coa", "cya")
    if("a" %in% output_kin) output_kin <- c(output_kin, "oa", "ya")
    if("n" %in% output_kin) output_kin <- c(output_kin, "nos", "nys")
    output_kin <- output_kin[!output_kin %in% c("s", "c", "a", "n")]
    output_kin <- match.arg(tolower(output_kin), all_possible_kin, several.ok = TRUE)
  }

  # if is time dependent or not
  age <- as.integer(rownames(p))
  years_data <- as.integer(colnames(p))
  if(time_invariant){
      if(!is.vector(p)) {
        output_period <- min(years_data)
        p <- p[,as.character(output_period)]
        f <- f[,as.character(output_period)]
      }
      kin_full <- kin_time_invariant(p = p, f = f, pi = pi,
                                     output_kin = output_kin, birth_female = birth_female) %>%
                              dplyr::mutate(cohort = NA, year = NA)
  }else{
      if(!is.null(output_cohort) & !is.null(output_period)) stop("sorry, you can not select cohort and period. Choose one please")
      kin_full <- kin_time_variant(p = p, f = f, pi = pi, n = n,
                              output_cohort = output_cohort, output_period = output_period,
                              output_kin = output_kin,
                              birth_female = birth_female)
      message(paste0("Assuming stable population before ", min(years_data), "."))
  }

  # re-group if grouped type is asked
  if(!is.null(output_kin_asked) & length(output_kin_asked)!=length(output_kin)){
    if("s" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("os", "ys")]   <- "s"
    if("c" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("coa", "cya")] <- "c"
    if("a" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("oa", "ya")]   <- "a"
    if("n" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("nos", "nys")] <- "n"
    kin_full <- kin_full %>%
      dplyr::summarise(living = sum(living), dead = sum(dead),
                       .by = c(kin, age_kin, age_focal, cohort, year))
  }

  # select period/cohort/age
  if(!is.null(output_age_focal) & all(output_age_focal %in% 1:120)){
    kin_full <- kin_full %>% dplyr::filter(age_focal %in% output_age_focal)
  }
  if(!is.null(output_cohort)){
    agrupar <- "cohort"
  } else if(!is.null(output_period)){
    agrupar <- "year"
  } else{
    agrupar <- c("year", "cohort")
  }
  agrupar_no_age_focal <- c("kin", agrupar)
  agrupar <- c("age_focal", "kin", agrupar)

  # get summary indicators based on group variables. If it is asked
  if(summary_kin){
    kin_summary <- dplyr::bind_rows(
      kin_full %>%
        dplyr::rename(total=living) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
        dplyr::summarise(count_living = sum(total),
                         mean_age = sum(total*age_kin)/sum(total),
                         sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
        tidyr::pivot_longer(count_living:sd_age, names_to = "indicator", values_to = "value"),
      kin_full %>%
        dplyr::rename(total=dead) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
        dplyr::summarise(count_dead = sum(total)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(agrupar_no_age_focal))) %>%
        dplyr::mutate(count_cum_dead = cumsum(count_dead),
                      mean_age_lost = cumsum(count_dead * age_focal)/cumsum(count_dead)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_longer(count_dead:mean_age_lost, names_to = "indicator", values_to = "value")) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value)

    # return
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)
  }else{
    # return
    kin_out <- kin_full
  }

  return(kin_out)
}

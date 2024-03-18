#' Estimate kin counts in a two-sex framework

#' @description Implementation of two-sex matrix kinship model. This produces kin counts grouped by kin, age and sex of
#' each relatives at each Focal´s age. For example, male cousins from aunts and uncles from different sibling's parents
#' are grouped in one male count of cousins. Note that the output labels relative following female notation: the label `m`
#' refers to either mothers or fathers, and column `sex_kin` determine the sex of the relative.
#' @details See Caswell (2022) for details on formulas.
#' @param pf numeric. A vector (atomic) or  matrix with female probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param pm numeric. A vector (atomic) or  matrix with male probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param ff numeric. Same as `pf` but for fertility rates.
#' @param fm numeric. Same as `pm` but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates. Default `TRUE`.
#' @param sex_focal character. "f" for female or "m" for male.
#' @param pif numeric. For using some specific age distribution of childbearing for mothers (same length as ages). Default `NULL`.
#' @param pim numeric. For using some specific age distribution of childbearing for fathers (same length as ages). Default `NULL`.
#' @param Hf numeric. A list where each list element (being the name of each list element the year) contains a matrix with cause-specific hazards for females with rows as causes and columns as ages, being the name of each col the age.
#' @param Hm numeric. A list where each list element (being the name of each list element the year) contains a matrix with cause-specific hazards for males with rows as causes and columns as ages, being the name of each col the age.
#' @param nf numeric. Only for `time_invariant = FALSE`. Same as `pf` but for population distribution (counts or `%`). Optional.
#' @param nm numeric. Only for `time_invariant = FALSE`. Same as `pm` but for population distribution (counts or `%`). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param output_age_focal integer. Vector of ages to select (and make faster the run).
#' @param birth_female numeric. Female portion at birth. This multiplies `f` argument. If `f` is already for female offspring, this needs to be set as 1.
#' @param summary_kin logical. Whether or not include `kin_summary` table (see output details). Default `TRUE`.
#' @return A list with:
#' \itemize{
#'  \item{kin_full}{ a data frame with year, cohort, Focal´s age, related ages and type of kin (for example `d` could be daughter or son depending `sex_kin`,
#'  `oa` is older aunts or uncles also depending `sex_kin` value, etc.), including living and dead kin at that age.}
#'  \item{kin_summary}{ a data frame with Focal´s age, related ages, sex and type of kin, with indicators obtained processing `kin_full`, grouping by cohort or period (depending on the given arguments):}
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
#' # Kin expected count by relative sex for a French female based on 2012 rates.
#' fra_fert_f <- fra_asfr_sex[,"ff"]
#' fra_fert_m <- fra_asfr_sex[,"fm"]
#' fra_surv_f <- fra_surv_sex[,"pf"]
#' fra_surv_m <- fra_surv_sex[,"pm"]
#' fra_2012 <- kin2sex(fra_surv_f, fra_surv_m, fra_fert_f, fra_fert_m)
#' head(fra_2012$kin_summary)
#'
# get kin ----------------------------------------------------------------
kin2sex <- function(pf = NULL, pm = NULL, ff = NULL, fm = NULL,
                 time_invariant = TRUE,
                 sex_focal = "f",
                 birth_female = 1/2.04,
                 pif = NULL, pim = NULL,
                 nf = NULL, nm = NULL,
                 Hf = NULL, Hm = NULL,
                 output_cohort = NULL, output_period = NULL, output_kin=NULL,output_age_focal = NULL,
                 summary_kin = TRUE)
  {

  # global vars
  living<-age_focal<-cohort<-year<-total<-mean_age<-count_living<-sd_age<-count_dead<-mean_age_lost<-indicator<-value<-sex_kin<-age_kin<-dead<-NULL
  age <- as.integer(rownames(pf))
  years_data <- as.integer(colnames(pf))

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

  # is cause of death specific or not
  is_cod <- !is.null(Hf) & !is.null(Hm)

  # if time dependent or not
  if(time_invariant){
      if(!is.vector(pf)) {
        output_period <- min(years_data)
        pf <- pf[,as.character(output_period)]
        pm <- pm[,as.character(output_period)]
        ff <- ff[,as.character(output_period)]
        fm <- fm[,as.character(output_period)]
      }
      if(is_cod){
        kin_full <- kin_time_invariant_2sex_cod(pf, pm, ff, fm,
                                            sex_focal = sex_focal,
                                            birth_female = birth_female,
                                            pif = pif, pim = pim,
                                            Hf = Hf, Hm = Hm,
                                            output_kin = output_kin) %>%
          dplyr::mutate(cohort = NA, year = NA)
      }else{
        kin_full <- kin_time_invariant_2sex(pf, pm, ff, fm,
                                            sex_focal = sex_focal,
                                            birth_female = birth_female,
                                            pif = pif, pim = pim,
                                            output_kin = output_kin) %>%
          dplyr::mutate(cohort = NA, year = NA)
      }

  }else{
      if(!is.null(output_cohort) & !is.null(output_period)) stop("sorry, you can not select cohort and period. Choose one please")
    if(is_cod){
      kin_full <- kin_time_variant_2sex_cod(pf = pf, pm = pm,
                                        ff = ff, fm = fm,
                                        sex_focal = sex_focal,
                                        birth_female = birth_female,
                                        pif = pif, pim = pim,
                                        nf = nf, nm = nm,
                                        Hf = Hf, Hm = Hm,
                                        output_cohort = output_cohort, output_period = output_period,
                                        output_kin = output_kin)
    }else{
      kin_full <- kin_time_variant_2sex(pf = pf, pm = pm,
                                        ff = ff, fm = fm,
                                        sex_focal = sex_focal,
                                        birth_female = birth_female,
                                        pif = pif, pim = pim,
                                        nf = nf, nm = nm,
                                        output_cohort = output_cohort, output_period = output_period,
                                        output_kin = output_kin)
    }
      message(paste0("Assuming stable population before ", min(years_data), "."))
  }

  # reorder
  kin_full <- kin_full %>% dplyr::select(year, cohort, age_focal, sex_kin, kin, age_kin, living, starts_with("dea"))

  # re-group if grouped type is asked
  if(!is.null(output_kin_asked) & length(output_kin_asked)!=length(output_kin)){
    if("s" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("os", "ys")]   <- "s"
    if("c" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("coa", "cya")] <- "c"
    if("a" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("oa", "ya")]   <- "a"
    if("n" %in% output_kin_asked) kin_full$kin[kin_full$kin %in% c("nos", "nys")] <- "n"
    kin_full <- kin_full %>%
      dplyr::group_by(kin, age_kin, age_focal, sex_kin, cohort, year) %>%
      dplyr::summarise_at(vars(c("living", dplyr::starts_with("dea"))), funs(sum)) %>%
      dplyr::ungroup()
  }

  # summary
  # select period/cohort/ge
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
  agrupar_no_age_focal <- c("kin", "sex_kin", agrupar)
  agrupar <- c("age_focal", "kin", "sex_kin", agrupar)

  # only return summary if is asked and is not cod
  if(summary_kin & !is_cod){
    kin_summary <- dplyr::bind_rows(
      as.data.frame(kin_full) %>%
        dplyr::rename(total=living) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(agrupar))) %>%
        dplyr::summarise(count_living = sum(total),
                         mean_age = sum(total*age_kin)/sum(total),
                         sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
        tidyr::pivot_longer(count_living:sd_age, names_to = "indicator", values_to = "value"),
      as.data.frame(kin_full) %>%
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
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)
  }else{
    kin_out <- kin_full
  }

  return(kin_out)
}

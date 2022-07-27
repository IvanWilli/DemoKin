#' Estimate kin counts

#' @description Implementation of Goodman-Keyfitz-Pullum equations in a stable or not framework.
#' @details See Caswell (2019) and Caswell (2021) for details on formulas.
#' @param U numeric. A vector or  matrix with probabilities (or survival ratios, or transition between age class in a more general perspective) with rows as ages (and columns as years in case of matrix, being the name of each col the year).
#' @param f numeric. Same as U but for fertility rates.
#' @param time_invariant logical. Constant assumption for a given `year` rates.
#' @param N numeric. Same as U but for population distribution (counts or `%`). Optional.
#' @param pi numeric. Same as U but for childbearing distribution (sum to 1). Optional.
#' @param output_cohort integer. Vector of year cohorts for returning results. Should be within input data years range.
#' @param output_period integer. Vector of period years for returning results. Should be within input data years range.
#' @param output_kin character. kin types to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth.
#' @param living logical. Only living kin counts `TRUE`, or including death kin `FALSE`.
#' @return A list with:
#'  * `kin`: a data frame with focal´s age, related ages and type of kin (for example `d` is daughter, `oa` is older aunts, etc.). The content (living or deaths) dependes on `alive` param.
#'  * A data frame with kin at each age of focal´s life. The content (living or deaths) dependes on `alive` param.
#'  * A data frame with kin at actual age of focal. The content (living or deaths) dependes on `alive` param.
#'  * Mean age of each type of kin at actual age of focal.  The content (living or deaths) dependes on `alive` param.
#'  * Total of each type of kin at actual age of focal.  The content (living or deaths) dependes on `alive` param.
#' @export
#' @examples
#' # If focal is 30 years old and lives in 1950. How much live kin would have if
#' # her relatives would experience mortality and fertility in that calendar year?
#' \dontrun{
#' swe30_1950_stable <- kin(U = swe_surv, f = swe_asfr, focal_year = 1950, stable = TRUE,selected_kin = c("m","gm"))
#' # How much live mothers and grandmothers would have if her relatives would experience
#' # mortality and fertility at each observed year?
#' swe30_1950_nonstable <- kin(U = swe_surv, f = swe_asfr,N = swe_pop,
#'                              stable = FALSE, focal_year = 1950, selected_kin = c("m","gm"))
#' # Difference in total by kin:
#' swe30_1950_stable$kin_by_age_focal %>% filter(age_focal==30) %>% select(kin, total)
#' swe30_1950_nonstable$kin_by_age_focal %>% filter(age_focal==30) %>% select(kin, total)
#' }
#'
# get kin ----------------------------------------------------------------
kin <- function(U = NULL, f = NULL,
                 time_invariant = TRUE,
                 N = NULL, pi = NULL,
                 output_cohort = NULL, output_period = NULL, output_kin=NULL,
                 birth_female = 1/2.04,
                 living = TRUE)
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
        focal_year <- min(years_data)
        U <- U[,as.character(focal_year)]
        f <- f[,as.character(focal_year)]
      }
      kin_full <- kin_time_invariant(U = U, f = f,
                                     output_kin = output_kin, birth_female = birth_female) %>%
                              mutate(cohort = NA, year = NA)
  }else{
      kin_full <- kin_time_variant(U = U, f = f, N = N, pi = pi,
                              output_cohort = output_cohort, output_period = output_period,
                              output_kin = output_kin,
                              birth_female = birth_female)
      message(paste0("Assuming stable population before ", min(years_data), "."))
  }

  # reorder
  kin_full <- kin_full %>% select(year, cohort, age_focal, kin, age_kin, alive, count)

  # return
  kin_summary <- kin_full %>%
    filter(alive=="yes") %>%
    rename(total=count) %>%
    group_by(year, cohort, age_focal, kin) %>%
    summarise(count    = sum(total),
              mean_age = sum(total*age_kin)/sum(total),
              sd_age  = (sum(total*age_kin^2)/sum(total)-mean_age^2)^.5) %>%
    ungroup()
  if(living){
    kin_out <- list(kin_full=kin_full %>% filter(alive=="yes"), kin_summary=kin_summary)
  }else{
    kin_summary <- kin_full %>%
      rename(total=count) %>%
      group_by(cohort, age_focal, kin, alive) %>%
      summarise(count = sum(total)) %>%
      ungroup() %>%
      group_by(cohort, kin, alive) %>%
      mutate(count_cum = cumsum(count),
             mean_age_lost = cumsum(count*age_focal)/cumsum(count)) %>%
      ungroup()
    kin_out <- list(kin_full=kin_full, kin_summary=kin_summary)
  }
  return(kin_out)
}

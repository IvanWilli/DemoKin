#' Estimate kin counts

#' @description Implementation of Goodman-Keyfitz-Pullum equations in a stable or not framework.
#' @details See Caswell (2019) and Caswell (2021) for details on formulas.
#' @param U numeric. A matrix of survival ratios with rows as ages and columns as years. The name of each col must be the year.
#' @param f numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. The name of each col must be the year.
#' @param N numeric. A matrix of population with rows as ages and columns as years. The name of each col must be the year.
#' @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. The name of each col must be the year.
#' @param ego_cohort integer. Year of birth of ego. Could be a vector. Should be within input data years range.
#' @param ego_year integer. Year of ego. Could be a vector. Should be within input data years range.
#' @param selected_kins character. Kins to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth.
#' @param Pb logic. Is given Pb as the first row in P?. If not, takes `P(0,1)` as `P(b,1)`. Useful for fertility matrix first row. Default `FALSE`.
#' @param stable logical. Stable assumption given `year` rates, taking ego_year in case is U and f are matrix.
#' @param living logial. Only living kin counts `TRUE`, or including death kins `FALSE`.
#' @return A list with:
#'  * `kins`: a data frame with ego´s age `x`, related ages `x_kin` and type of kin (for example `d` is daughter, `oa` is older aunts, etc.). The content (living or deaths) dependes on `alive` param.
#'  * A data frame with kins at each age of ego´s life. The content (living or deaths) dependes on `alive` param.
#'  * A data frame with kins at actual age of ego. The content (living or deaths) dependes on `alive` param.
#'  * Mean age of each type of kin at actual age of ego.  The content (living or deaths) dependes on `alive` param.
#'  * Total of each type of kin at actual age of ego.  The content (living or deaths) dependes on `alive` param.
#' @export
#' @examples
#' # If ego is 30 years old and lives in 1950. How much live kins would have if
#' # her relatives would experience mortality and fertility in that calendar year?
#' \dontrun{
#' swe30_1950_stable <- kins(U = swe_surv, f = swe_asfr, ego_year = 1950, stable = TRUE,selected_kins = c("m","gm"))
#' # How much live mothers and grandmothers would have if her relatives would experience
#' # mortality and fertility at each observed year?
#' swe30_1950_nonstable <- kins(U = swe_surv, f = swe_asfr,N = swe_pop,
#'                              stable = FALSE, ego_year = 1950, selected_kins = c("m","gm"))
#' # Difference in total by kin:
#' swe30_1950_stable$kins_by_age_ego %>% filter(age_ego==30) %>% select(kin, total)
#' swe30_1950_nonstable$kins_by_age_ego %>% filter(age_ego==30) %>% select(kin, total)
#' }
#'
# get kins ----------------------------------------------------------------
kins <- function(U = NULL, f = NULL, N = NULL, pi = NULL,
                 stable = TRUE,
                 ego_cohort = NULL, ego_year = NULL, selected_kins=NULL,
                 birth_female = 1/2.04,
                 Pb = FALSE,
                 living = TRUE)
  {

  age = as.integer(rownames(U))
  # if stable or not
  if(stable){
      if(!is.vector(U)) {
        U = U[,as.character(ego_year)]
        f = f[,as.character(ego_year)]
      }
      kins <- kins_stable(U = U, f = f, selected_kins=selected_kins, birth_female = birth_female) %>%
        mutate(cohort = 0, year = 0)
  }else{
      kins <- kins_non_stable(U = U, f = f, N = N, pi = pi,
                              ego_cohort = ego_cohort, ego_year = ego_year,
                              selected_kins=selected_kins,
                              birth_female = birth_female,
                              Pb = Pb)
  }

    # return
    kins_by_age_ego <- kins %>%
      filter(alive=="yes") %>%
      group_by(year, cohort, age_ego, kin) %>%
      summarise(total    = sum(count),
                mean_age = sum(count*age_kin)/sum(count),
                sd_age  = (sum(count*age_kin^2)/sum(count)-mean_age^2)^.5) %>%
      ungroup()
    if(living){
      kins_out <- list(kins=kins %>% filter(alive=="yes"), kins_by_age_ego=kins_by_age_ego)
    }else{
      kins_alive_death <- kins %>%
        group_by(cohort, age_ego, kin, alive) %>%
        summarise(total = sum(count)) %>%
        ungroup() %>%
        group_by(cohort, kin, alive) %>%
        mutate(total_cum = cumsum(total),
               mean_age_lost = cumsum(total*age_ego)/cumsum(total)) %>%
        ungroup()
      # add lost mean age
      kins_out <- list(kins=kins, kins_by_age_ego=kins_by_age_ego, kins_alive_death=kins_alive_death)
    }
  return(kins_out)
}

#' Estimate kin counts

#' @description Implementation of Goodman-Keyfitz-Pullum equations in an stable or not framework.
#' @details See Caswell (2019) for details on stable option.
#' @param ego_age integer. Age where ego is.
#' @param year integer. Year where ego is with `ego_age` age.
#' @param P numeric. A matrix of survival ratios with rows as ages and columns as years.
#' The name of each col must be the year.
#' @param asfr numeric. A matrix of age-specific fertility rates with rows as ages and
#' columns as years. The name of each col must be the year.
#' @param N numeric. A matrix of population with rows as ages and columns as years.
#' The name of each col must be the year.
#' @param age integer. Ages, assuming last one as an open age group.
#' @param birth_female numeric. Female portion at birth.
#' @param stable logical. Stable assumption given `year` rates.
#'
#' @return A list with:
#'  *
#'  * A data frame with ego´s age `x`, related ages `x_kin` and kind of kin (for example `d` is daughter, `oa` is older aunts, etc.).
#'  * A data frame with available kins at each age of ego.
#'  * A data frame with available kins at actual age of ego.
#'  * Mean age of each type of kin at actual age of ego.
#'  * Total of each type of kin at actual age of ego.
#' @export
#' @examples
#' # If ego is 30 years old and lives in 1950. How much live kins would have if
#' # her relatives would experience mortality and fertility in that calendar year?
#' \dontrun{
#' swe30_1950_stable <- kins(ego_age = 30, year = 1950,
#'                           P = swe_surv, asfr = swe_asfr,
#'                           stable = TRUE)
#' # How much live kins would have if her relatives would experience
#' # mortality and fertility at each observed year?
#' swe30_1950_nonstable <- kins(ego_age = 30, year = 1950,
#'                           P = swe_surv, asfr = swe_asfr,
#'                           stable = FALSE)
#' # Difference in total by kin:
#' swe30_1950_stable[["kins_total"]] - swe30_1950_nonstable[["kins_total"]]
#' }
#'
# get kins ----------------------------------------------------------------
kins <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                     P = NULL, asfr = NULL, N = NULL,
                     stable = FALSE,
                     age = 0:100,
                     birth_female = 1/2.04)
  {

  # if stable or not
  if(stable){
      kins <- kins_stable(p = P[,as.character(year)],
                          f = asfr[,as.character(year)],
                          cum_deaths = FALSE) %>%
              filter(x <= ego_age)
  }else{
      kins <- kins_non_stable(ego_age = ego_age, year = year,
                              P = P, asfr = asfr, N = N)
  }

  # summarise results
  kins_by_age_ego <- kins %>% group_by(x) %>% select(-x_kin) %>% summarise_all(sum)
  kins_by_age_kin <- kins[kins$x == ego_age,]
  kins_mean_age   <- colSums(kins_by_age_kin[,3:ncol(kins)]*0:100)/colSums(kins_by_age_kin[,3:ncol(kins)])
  kins_sd_age     <- colSums(kins_by_age_kin[,3:ncol(kins)]*(0:100)^2)/colSums(kins_by_age_kin[,3:ncol(kins)]) - kins_mean_age^2
  kins_total      <- colSums(kins_by_age_kin[,c(3:ncol(kins))])
  return(list(kins = kins,
              kins_by_age_ego = kins_by_age_ego,
              kins_by_age_kin = kins_by_age_kin,
              kins_mean_age = kins_mean_age,
              kins_total = kins_total))
}

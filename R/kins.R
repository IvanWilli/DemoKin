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
#' @param alive character. Only living kin counts `yes`, death kins `no`, or both other char.
#' @return A list with:
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
                     birth_female = 1/2.04,
                     alive = "yes")
  {

  # if stable or not
  if(stable){
      kins <- kins_stable(P = P[,as.character(year)],
                          asfr = asfr[,as.character(year)],
                          birth_female = birth_female) %>%
              filter(x <= ego_age)
  }else{
      kins <- kins_non_stable(ego_age = ego_age, year = year,
                              P = P, asfr = asfr, N = N,
                              birth_female = birth_female)
  }

  # living results
  alive_yes <- kins %>% filter(alive=="yes")
    alive_yes$alive <- NULL
    kins_by_age_ego <- alive_yes %>% group_by(x) %>% select(-x_kin) %>% summarise_all(sum)
    kins_by_age_kin <- alive_yes[alive_yes$x == ego_age,]
    kins_mean_age   <- colSums(kins_by_age_kin[,3:ncol(alive_yes)]*0:100)/colSums(kins_by_age_kin[,3:ncol(alive_yes)])
    kins_var_age    <- colSums(kins_by_age_kin[,3:ncol(alive_yes)]*(0:100)^2)/colSums(kins_by_age_kin[,3:ncol(alive_yes)]) - kins_mean_age^2
    kins_total      <- colSums(kins_by_age_kin[,c(3:ncol(alive_yes))])
    out_yes <- list(kins = alive_yes,
                    kins_by_age_ego = kins_by_age_ego,
                    kins_by_age_kin = kins_by_age_kin,
                    kins_mean_age = kins_mean_age,
                    kins_total = kins_total)

  # death results
  alive_no <- kins %>% filter(alive=="no")
    alive_no$alive <- NULL
    freq_d_by_age_ego <- alive_no %>% group_by(x) %>% select(-x_kin) %>% summarise_all(sum)
    cum_d_by_age_ego  <- freq_d_by_age_ego %>% ungroup() %>%
                                summarise_at(.vars = 2:ncol(freq_d_by_age_ego),cumsum) %>%
                                mutate(x=freq_d_by_age_ego$x)
    cum_d_total       <- cum_d_by_age_ego %>% filter(x == ego_age) %>% select(-x)
    lost_mean_age     <- colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)]*freq_d_by_age_ego$x)/
                         colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)])
    lost_var_age      <- colSums(freq_d_by_age_ego[,2:ncol(freq_d_by_age_ego)]*freq_d_by_age_ego$x^2)/lost_mean_age^2
    out_no <- list(kins = alive_no,
                   kins_death_by_age_ego = freq_d_by_age_ego,
                   kins_cum_death_by_age_ego = cum_d_by_age_ego,
                   kins_cum_death_total = cum_d_total,
                   lost_mean_age = lost_mean_age,
                   lost_var_age  = lost_var_age)

  # not return all
  if(alive=="yes"){
    kins=out_yes
  }else if(alive=="no"){
    kins=out_no
  } else{
    kins=list(kins_living=out_yes,kins_death=out_no)
  }
  return(kins)
}

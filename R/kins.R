# estimate kin counts by age for female population, given an age and year for ego.
# U, f and W are lists of matrix (age x age) probabilities with year as name

# get kins ----------------------------------------------------------------
kins <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                     P = NULL, asfr = NULL, N = NULL,
                     stable = FALSE,
                     age = 0:100,
                     birth_female = 1/2.04,
                     verbose=T)
  {

  # if stable or not
  if(stable){
      kins <- kins_stable(p = P[,as.character(year)],
                          f = asfr[,as.character(year)]) %>%
              filter(x <= ego_age) %>%
              spread(kin,count)
  }else{
    kins <- kins_non_stable(ego_age = ego_age, year = year,
                            P = P, asfr = asfr, N = N,
                            stable = stable)
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

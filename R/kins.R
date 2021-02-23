# estimate kin counts by age for female population, given an age and year for ego.
# U, f and W are lists of matrix (age x age) probabilities with year as name

# get kins ----------------------------------------------------------------
kins <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                     U = NULL, f = NULL, W = NULL,
                     age = 0:100,
                     birth_female = 1/2.04,
                     verbose=T)
  {

  # half year and half age
  years_data  <- as.integer(names(U))
  ages        <- length(age)
  w           <- last(age)
  ego_cohort  <- year - ego_age

  # complete data on the right(left) with last (first) available year
  # sonn forecats option here
  if((year-1) > max(years_data)){
    for(y in (max(years_data)+1):(year-1)){
      U[[as.character(y)]] = U[[as.character(max(years_data))]]
      f[[as.character(y)]] = f[[as.character(max(years_data))]]
      W[as.character(y)]   = W[,as.character(max(years_data))]
    }
  }
    for(y in 1700:(min(years_data)-1)){
      U[[as.character(y)]] = U[[as.character(min(years_data))]]
      f[[as.character(y)]] = f[[as.character(min(years_data))]]
      W[as.character(y)]   = W[,as.character(min(years_data))]
    }

  # female births only
  f = lapply(f, function(x) x * birth_female)

  # index martix
  e = matrix(0, ages, ages+1)
  diag(e) = 1

  # conditional matrix on mother´s age (cols) at ego´s birth
  osM = nosM = oaM = coaM = yaM = cyaM = osM = gmMy = gmM = ggmMy = ggmM = matrix(0, ages, ages)

  # conditional matrix on grandmother´s age (cols) at mothers´s birth
  oaMy  = coaMy  = yaMy  = cyaMy  = osMy = nosMy = matrix(0, ages, ages)

  # mother´s age at ego´s birth
  for(m_age in age){
    m          = W[,as.character(ego_cohort)]
    m_cohort   = ego_cohort - m_age - 1
    gm         = W[,as.character(m_cohort)]
    ya = cya = os = nos = rep(0,ages)

    for(y in m_cohort:ego_cohort){
      Ut = U[[as.character(y)]]
      ft = f[[as.character(y)]]
      gm = Ut %*% gm
      ya = Ut %*% ya + ft %*% gm
      cya = Ut %*% cya + ft %*% ya
      os = Ut %*% os + ft %*% e[, y - m_cohort + 1]
      nos = Ut %*% nos + ft %*% os
    }
    # conditionated to mother´s age
    gmM[, m_age+1] = gm
    yaM[, m_age+1] = ya
    cyaM[,m_age+1] = cya
    osM[, m_age+1] = os
    nosM[,m_age+1] = nos

    # grandmother´s age at mother´s birth
    for(gm_age in age){
      gm_cohort  = m_cohort - gm_age - 1
      ggm        = W[,as.character(gm_cohort)]
      oa = coa = rep(0,ages)

      # before mother born
      for(y in gm_cohort:m_cohort){
        Ut = U[[as.character(y)]]
        ft = f[[as.character(y)]]
        ggm = Ut %*% ggm
        oa = Ut %*% oa + ft %*% e[,y - gm_cohort + 1]
        coa = Ut %*% coa + ft %*% oa
      }

      # after mother born
      for(y in m_cohort:ego_cohort){
        Ut = U[[as.character(y)]]
        ft = f[[as.character(y)]]
        ggm = Ut %*% ggm
        oa = Ut %*% oa
        coa = Ut %*% coa + ft %*% oa
      }
      # conditionated to granmother´s age
      ggmMy[,gm_age+1] = ggm
      oaMy[,gm_age+1] = oa
      coaMy[,gm_age+1] = coa
    }
    # expected count for possible granmother´s age
    ggmM[,m_age+1] = ggmMy %*% gm
    oaM[,m_age+1]  = oaMy  %*% gm
    coaM[,m_age+1] = coaMy %*% gm
  }
  # expected count for possible mother´s age
  gm  = gmM  %*% m
  ggm = ggmM %*% m
  oa  = oaM  %*% m
  coa = coaM %*% m
  ya  = yaM  %*% m
  cya = cyaM %*% m
  os  = osM  %*% m
  nos = nosM %*% m

  # ego´s trip
  e = matrix(0, ages, ages)

  # initial descendants
  d = gd = ys = nys = matrix(rep(0, ages), byrow = F)

  # collect kins
  kins = data.frame(x=0, x_kin = age,
                    ggm=ggm,
                    gm=gm,
                    oa=oa,
                    m=m,
                    ya=ya,
                    coa=coa, cya=cya,
                    os=os, ys=ys,
                    nos=nos, nys=nys,
                    d=d, gd=gd)

  for(x in 1:ego_age){
    Ut = U[[as.character(ego_cohort + x - 1)]]
    ft = f[[as.character(ego_cohort + x - 1)]]
    ggm = Ut %*% ggm
    gm  = Ut %*% gm
    oa  = Ut %*% oa
    m   = Ut %*% m
    ya  = Ut %*% ya  + ft %*% gm
    coa = Ut %*% coa + ft %*% oa
    cya = Ut %*% cya + ft %*% ya
    os  = Ut %*% os
    ys  = Ut %*% ys  + ft %*% m
    nos = Ut %*% nos + ft %*% os
    nys = Ut %*% nys + ft %*% ys
    d   = Ut %*% d   + ft %*% e[,x]
    gd  = Ut %*% gd  + ft %*% d

    # bind results
    kins <- rbind(kins, data.frame(x=x, x_kin = 0:100,
                                  ggm=ggm,
                                  gm=gm,
                                  oa=oa,
                                  m=m,
                                  ya=ya,
                                  coa=coa, cya=cya,
                                  os=os, ys=ys,
                                  nos=nos, nys=nys,
                                  d=d, gd=gd))

  }

  # summarise results
  kins_by_age_ego <- kins %>% group_by(x) %>% summarise_at(vars(ggm:gd),sum)
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

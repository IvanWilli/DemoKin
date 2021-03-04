

kins_non_stable <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                            P = NULL, asfr = NULL, N = NULL,
                            stable = FALSE,
                            age = 0:100,
                            birth_female = 1/2.04){

  # half year and half age
  years_data  <- as.integer(colnames(P))
  ages        <- length(age)
  w           <- last(age)
  ego_cohort  <- year - ego_age

  # get lists
  # arrange lists of matrixs
  U = f = list()
  for(t in 1:length(years_data)){
    Ut = ft = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- P[-101,t]
    Ut[ages, ages]=P[101,t]
    ft[1,] = asfr[,t]
    U[[as.character(years_data[t])]] <- Ut
    f[[as.character(years_data[t])]] <- ft
  }

  # age distribution at childborn
  W <- t(t(N * asfr)/colSums(N * asfr))

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
    W   = cbind(W, W[,as.character(min(years_data))])
    colnames(W)[ncol(W)] = as.character(y)
  }

  # female births only
  f = lapply(f, function(x) x * birth_female)

  # conditional matrix on mother´s age (cols) at ego´s birth
  osM = nosM = oaM = coaM = yaM = cyaM = osM = gmMy = gmM = ggmMy = ggmM = matrix(0, ages, ages)

  # conditional matrix on grandmother´s age (cols) at mothers´s birth
  oaMy  = coaMy  = yaMy  = cyaMy  = osMy = nosMy = matrix(0, ages, ages)

  # index martix
  e = matrix(0, ages, ages+1)
  diag(e) = 1

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
  diag(e[1:ages,1:ages]) = 1

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
  return(kins)
}

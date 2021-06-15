#' Estimate kin counts in a non stable framework

#' @description Implementation of Goodman-Keyfitz-Pullum equations but as a
#' weigthed average of possible ages at childbear of mothers and grandmothers of ego.

#' @param ego_age integer. Age where ego is.
#' @param year integer. Year where ego is with `ego_age` age.
#' @param P numeric. A matrix of survival ratios with rows as ages and columns as years. The name of each col must be the year.
#' @param asfr numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. The name of each col must be the year.
#' @param N numeric. A matrix of population with rows as ages and columns as years. The name of each col must be the year.
#' @param age integer. Ages, assuming last one as an open age group.
#' @param birth_female numeric. Female portion at birth.
#' @param Pb logic. Is given Pb as the first row in P?. If not take `P(0,1)` as `P(b,1)`. Default `FALSE`.
#'
#' @return A data frame with ego´s age `x`, related ages `x_kin` and type of kin
#' (for example `d` is daughter, `oa` is older aunts, etc.), alive and death.
#' @export

kins_non_stable <- function(ego_age = NULL, year = NULL, # where ego is, it will be at half year
                            P = NULL, asfr = NULL, N = NULL,
                            age = 0:100,
                            birth_female = 1/2.04,
                            Pb = FALSE){


  # check input
  stopifnot(!is.null(P)&!is.null(asfr)&!is.null(N))

  # diff years
  if(!any(as.integer(colnames(P))==as.integer(colnames(asfr))))stop("Data should be from same years.")

  # half year and half age
  years_data  <- as.integer(colnames(P))
  ages        <- length(age)
  w           <- last(age)
  ego_cohort  <- year - ego_age
  zeros = matrix(0, nrow=ages, ncol=ages)
  if(Pb){
    stopifnot(length(years_data)==ncol(Pb))
  }else{
    Pb = P[1,]
  }

  # get lists
  # arrange lists of matrixs
  U = f = list()
  for(t in 1:length(years_data)){
    Ut = Mt = Dcum = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- P[-101,t]
    Ut[ages, ages]=P[101,t]
    diag(Mt) = 1 - P[,t]
    # diag(Dcum) = 1
    U[[as.character(years_data[t])]] <- rbind(cbind(Ut,zeros),cbind(Mt,Dcum))
    ft = matrix(0, nrow=ages*2, ncol=ages*2)
    ft[1,1:ages] = asfr[,t] * birth_female * (1+P[,t])/2 * Pb[1,t]
    f[[as.character(years_data[t])]] <- ft
  }

  # age distribution at childborn
  W <- rbind(t(t(N * asfr)/colSums(N * asfr)),
             matrix(0,ages,length(years_data)))

  # complete data on the right(left) with last (first) available year
  cat(paste0("Rates before ",min(years_data)," assumed as constant."))
  for(y in 1500:(min(years_data)-1)){
    U[[as.character(y)]] = U[[as.character(min(years_data))]]
    f[[as.character(y)]] = f[[as.character(min(years_data))]]
    W = cbind(W, W[,as.character(min(years_data))])
    colnames(W)[ncol(W)] = as.character(y)
  }
  if((year-1) > max(years_data)){
    cat(paste0("Rates after ",max(years_data)," assumed as constant."))
    for(y in (max(years_data)+1):(year-1)){
      U[[as.character(y)]] = U[[as.character(max(years_data))]]
      f[[as.character(y)]] = f[[as.character(max(years_data))]]
      W   = cbind(W, W[,as.character(max(years_data))])
      colnames(W)[ncol(W)] = as.character(y)
    }
  }

  # conditional matrix on mother´s age (cols) at ego´s birth
  osM = nosM = oaM = coaM = yaM = cyaM = osM = gmMy = gmM = ggmMy = ggmM = matrix(0, ages * 2, ages)

  # conditional matrix on grandmother´s age (cols) at mothers´s birth
  oaMy  = coaMy  = yaMy  = cyaMy  = osMy = nosMy = matrix(0, ages * 2, ages)

  # index martix
  e = matrix(0, ages * 2, ages * 2)
  diag(e[1:ages,1:ages]) = 1

  # mother´s age at ego´s birth
  for(m_age in age){
    m          = W[,as.character(ego_cohort)]
    m_cohort   = ego_cohort - m_age - 1
    gm         = W[,as.character(m_cohort)]
    ya = cya = os = nos = rep(0,ages*2)

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
      oa = coa = rep(0, ages * 2)

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
    ggmM[,m_age+1] = ggmMy %*% W[1:ages,as.character(m_cohort)]
    oaM[,m_age+1]  = oaMy  %*% W[1:ages,as.character(m_cohort)]
    coaM[,m_age+1] = coaMy %*% W[1:ages,as.character(m_cohort)]
  }
  # expected count for possible mother´s age
  gm  = gmM  %*% m[1:ages]
  ggm = ggmM %*% m[1:ages]
  oa  = oaM  %*% m[1:ages]
  coa = coaM %*% m[1:ages]
  ya  = yaM  %*% m[1:ages]
  cya = cyaM %*% m[1:ages]
  os  = osM  %*% m[1:ages]
  nos = nosM %*% m[1:ages]

  # ego´s trip

  # initial descendants
  d = gd = ys = nys = matrix(0,ages*2,1)

  # no matters deaths before ego borns
  ggm[(ages+1):(2*ages)] = 0
  gm[(ages+1) :(2*ages)] = 0
  m[(ages+1)  :(2*ages)] = 0
  oa[(ages+1) :(2*ages)] = 0
  ya[(ages+1) :(2*ages)] = 0
  coa[(ages+1):(2*ages)] = 0
  cya[(ages+1):(2*ages)] = 0
  os[(ages+1) :(2*ages)] = 0
  nos[(ages+1):(2*ages)] = 0

  # collect kins
  kins = data.frame(x=0, x_kin = rep(age,2),
                    alive = c(rep("yes",ages),rep("no",ages)),
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
    kins <- rbind(kins, data.frame(x=x, x_kin = rep(age,2),
                                   alive = c(rep("yes",ages),rep("no",ages)),
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

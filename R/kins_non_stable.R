#' Estimate kin counts in a non stable framework

#' @description Implementation of non-stable Goodman-Keyfitz-Pullum equations based on Caswell (2021).
#'
#' @param U numeric. A matrix of survival ratios with rows as ages and columns as years. The name of each col must be the year.
#' @param f numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. The name of each col must be the year.
#' @param N numeric. A matrix of population with rows as ages and columns as years. The name of each col must be the year.
#' @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. The name of each col must be the year.
#' @param ego_cohort integer. Year of birth of ego. Could be a vector. Should be within input data years range.
#' @param ego_year integer. Year of ego. Could be a vector. Should be within input data years range.
#' @param selected_kins character. Kins to return: "m" for mother, "d" for daughter,...
#' @param birth_female numeric. Female portion at birth.
#' @param Pb logic. Is given Pb as the first row in P?. If not, takes `P(0,1)` as `P(b,1)`. Useful for fertility matrix first row. Default `FALSE`.
#'
#' @return A data frame with ego´s age, related ages and type of kin
#' (for example `d` is daughter, `oa` is older aunts, etc.), alive and death.
#' @export

kins_non_stable <- function(U = NULL, f = NULL, N = NULL, pi = NULL,
                            ego_cohort = NULL, ego_year = NULL, selected_kins = NULL,
                            birth_female = 1/2.04,
                            Pb = FALSE){

  # check input
  stopifnot(!is.null(U)&!is.null(f)&any(!is.null(N)|!is.null(pi)))

  # diff years
  if(!any(as.integer(colnames(U))==as.integer(colnames(f))))stop("Data should be from same years.")

  # data should be from consequtive years
  years_data = as.integer(colnames(U))
  if(any(diff(years_data)!=1))stop("Data should be for consecutive years.")

  # half year and half age
  age          <- 0:(nrow(U)-1)
  n_years_data <- length(years_data)
  ages         <- length(age)
  om           <- max(age)
  zeros        <- matrix(0, nrow=ages, ncol=ages)
  if(Pb){
    stopifnot(length(years_data)==ncol(Pb))
  }else{
    Pb = U[1,]
  }

  # age distribution at childborn
  if(is.null(pi)){
    if(is.null(N)){
      stop("You need data on pi or N.")
    }else{
      pi <- rbind(t(t(N * f)/colSums(N * f)),
                  matrix(0,ages,length(years_data)))
    }
  }

  # get lists of matrix
  Ul = fl = list()
  for(t in 1:n_years_data){
    Ut = Mt = Dcum = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- U[-ages,t]
    Ut[ages, ages]=U[ages,t]
    diag(Mt) = 1 - U[,t]
    Ul[[as.character(years_data[t])]] <- rbind(cbind(Ut,zeros),cbind(Mt,Dcum))
    ft = matrix(0, nrow=ages*2, ncol=ages*2)
    ft[1,1:ages] = f[,t] * birth_female * (1+U[,t])/2 * as.numeric(Pb[t])
    fl[[as.character(years_data[t])]] <- ft
  }
  U <- Ul
  f <- fl

  # loop over years (more performance here)
  kins_all <- list()
  for (iyear in 1:n_years_data){
    # print(iyear)
    Ut <- as.matrix(U[[iyear]]);
    ft <- as.matrix(f[[iyear]]);
    pit <- pi[,iyear];
    if (iyear==1){
      U1 <- c(diag(Ut[-1,])[1:om],Ut[om,om])
      f1 <- ft[1,][1:ages]
      pi1 <- pit[1:ages]
      kins_all[[1]] <- kins_stable(U = U1, f = f1, pi = pi1, birth_female = birth_female, list_output = TRUE)
    }
    kins_all[[iyear+1]] <- timevarying_kin(Ut=Ut,ft=ft,pit=pit,ages,pkin=kins_all[[iyear]])
  }

  # filter years and kins that were selected
  names(kins_all) <- as.character(years_data)
  if(!is.null(ego_cohort)){
    selected_cohorts_year_age <- data.frame(age = rep(age,length(ego_cohort)),
                                            year = map(ego_cohort,.f = ~.x+age) %>%
                                              unlist(use.names = F))
  }else{selected_cohorts_year_age <- c()}
  if(!is.null(ego_year)){selected_years_age <- expand.grid(age, ego_year) %>% rename(age=1,year=2)
  }else{selected_years_age <- c()}
  if(is.null(ego_year) & is.null(ego_cohort)){
    ego_year = years_data
  }
  out_selected <- bind_rows(selected_years_age,selected_cohorts_year_age) %>% distinct()
  possible_kins <- c("d","gd","ggd","m","gm","ggm","os","ys","nos","nys","oa","ya","coa","cya")
  if(is.null(selected_kins)){
    selected_kins_position <- 1:length(possible_kins)
  }else{
    selected_kins_position <- which(possible_kins %in% selected_kins)
  }

  # first filter
  kins_all <- kins_all %>%
    keep(names(.) %in% as.character(unique(out_selected$year))) %>%
    map(~ .[selected_kins_position])

  # long format
  kins <- lapply(names(kins_all), function(Y){
    X <- kins_all[[Y]]
    X <- map2(X, names(X), function(x,y) as.data.frame(x) %>%
                mutate(year = Y,
                        kin=y,
                        age_kin = rep(age,2),
                        alive = c(rep("yes",ages), rep("no",ages)),
                        .before=everything())) %>%
      bind_rows() %>%
      setNames(c("year","kin","age_kin","alive",as.character(age))) %>%
      gather(age_ego, count,-age_kin, -kin, -year, -alive) %>%
      mutate(age_ego = as.integer(age_ego),
             year = as.integer(year),
             cohort = year - age_ego) %>%
      filter(age_ego %in% out_selected$age[out_selected$year==as.integer(Y)])}) %>%
    bind_rows()
  return(kins)
}

# main fun for each year: based on caswell matlab code

#' @description one time projection kin. internal function.
#'
#' @param Ut numeric. A matrix of survival ratios with rows as ages and columns as years. The name of each col must be the year.
#' @param ft numeric. A matrix of age-specific fertility rates with rows as ages and columns as years. The name of each col must be the year.
#' @param pi numeric. A matrix with distribution of childbearing with rows as ages and columns as years. The name of each col must be the year.
#' @param ages numeric.
#' @param pkin numeric. A list with kin count distribution in previous year.
#
timevarying_kin<- function(Ut,ft,pit,ages, pkin){

  # frequently used zero vector for initial condition
  zvec=rep(0,ages*2);
  I = matrix(0, ages * 2, ages * 2)
  diag(I[1:ages,1:ages]) = 1
  om=ages-1;
  d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,ages*2,ages)
  kins_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # initial distribution
  d[,1]=gd[,1]=ggd[,1]=ys[,1]=nys[,1]=zvec
  m[1:ages,1]  = pit[1:ages]
  gm[1:ages,1] = pkin[["m"]][1:ages,] %*% pit[1:ages]
  ggm[1:ages,1]= pkin[["gm"]][1:ages,] %*% pit[1:ages]
  os[1:ages,1] = pkin[["d"]][1:ages,] %*% pit[1:ages]
  ys[1:ages,1] = pkin[["gd"]][1:ages,] %*% pit[1:ages]
  oa[1:ages,1] = pkin[["os"]][1:ages,] %*% pit[1:ages]
  ya[1:ages,1] = pkin[["ys"]][1:ages,] %*% pit[1:ages]
  coa[1:ages,1]= pkin[["nos"]][1:ages,] %*% pit[1:ages]
  cya[1:ages,1]= pkin[["nys"]][1:ages,] %*% pit[1:ages]

  for (ix in 1:om){
    d[,ix+1]   = Ut %*% pkin[["d"]][,ix]   + ft %*% I[,ix]
    gd[,ix+1]  = Ut %*% pkin[["gd"]][,ix]  + ft %*% pkin[["d"]][,ix]
    ggd[,ix+1] = Ut %*% pkin[["ggd"]][,ix] + ft %*% pkin[["gd"]][,ix]
    m[,ix+1]   = Ut %*% pkin[["m"]][,ix]
    gm[,ix+1]  = Ut %*% pkin[["gm"]][,ix]
    ggm[,ix+1] = Ut %*% pkin[["ggm"]][,ix]
    os[,ix+1]  = Ut %*% pkin[["os"]][,ix]
    ys[,ix+1]  = Ut %*% pkin[["ys"]][,ix]  + ft %*% pkin[["m"]][,ix]
    nos[,ix+1] = Ut %*% pkin[["nos"]][,ix] + ft %*% pkin[["os"]][,ix]
    nys[,ix+1] = Ut %*% pkin[["nys"]][,ix] + ft %*% pkin[["ys"]][,ix]
    oa[,ix+1]  = Ut %*% pkin[["oa"]][,ix]
    ya[,ix+1]  = Ut %*% pkin[["ya"]][,ix]  + ft %*% pkin[["gm"]][,ix]
    coa[,ix+1] = Ut %*% pkin[["coa"]][,ix] + ft %*% pkin[["oa"]][,ix]
    cya[,ix+1] = Ut %*% pkin[["cya"]][,ix] + ft %*% pkin[["ya"]][,ix]
  }

  kins_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  return(kins_list)
}

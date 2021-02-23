get_HMDHFD <- function(country = "SWE",
                         min_year = 1900,
                         max_year = 2018,
                         user = NULL,
                         pass = NULL){
  # source HMD HFD - use SWE now-----------------------------------------------------------------
  library(HMDHFDplus)

  pop <- readHMDweb(CNTRY = country, "Population", user, pass, fixup = TRUE) %>%
          select(Year, Age, N = Female1)%>%
          filter(Year >= min_year, Year <= max_year)
  lt <- readHMDweb(country, "fltper_1x1", user, pass, fixup = TRUE) %>%
          filter(Year >= min_year, Year <= max_year)
  asfr <- readHFDweb(country, "asfrRR", user, "52962", fixup = TRUE)%>%
          filter(Year >= min_year, Year <= max_year)

  # list of yearly Leslie matrix ---------------------------------------------------

  age = 0:100
  ages = length(age)
  w = last(age)
  last_year = max(lt$Year)
  years = min_year:last_year

  # survival probabilities
  L <- lt %>%
    filter(Age<=w) %>%
    mutate(Lx = ifelse(Age==w, Tx, Lx)) %>%
    select(Year, Age, Lx) %>%
    spread(Year, Lx) %>%
    select(-Age)

  P <- rbind(L[c(-1,-101),]/L[-c(100:101),],
             L[101,]/(L[100,]+L[101,]),
             L[101,]/(L[100,]+L[101,]))

  # fertility
  asfr <- asfr %>%
    filter(Year >= min_year) %>%
    select(-OpenInterval) %>%
    rbind(
      expand.grid(Year = years,
                  Age = c(0:11,56:100),
                  ASFR = 0)) %>%
    arrange(Year, Age) %>%
    spread(Year, ASFR) %>%
    select(-Age)

  # population
  N <- pop %>%
    mutate(Age = ifelse(Age>100, 100, Age)) %>%
    group_by(Year, Age) %>% summarise(N = sum(N)) %>%
    filter(Age<=100, Year >= min_year) %>%
    arrange(Year, Age) %>%
    spread(Year, N) %>%
    select(-Age)

  # mean age at childborn
  W <- t(t(N * asfr)/colSums(N * asfr))

  # arrange lists of matrixs
  U = f = list()
  for(t in 1:length(years)){
    # t = 1
    Ut = ft = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- P[-101,t]
    Ut[ages, ages]=P[101,t]
    ft[1,] = asfr[,t]
    U[[as.character(years[t])]] <- Ut
    f[[as.character(years[t])]] <- ft

  }

  return(list(U=U, f=f, W=as.data.frame(W)))
}

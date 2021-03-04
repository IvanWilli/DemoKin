#' Get time serie matrix data from HMD/HFD

#' @description Wrapper function to get data.frames of survival, fertlity and population
#' of selected countriy on selected period.

#' @param country numeric. Country code from rom HMD/HFD.
#' @param min_year numeric. Latest year to get data.
#' @param min_year integer. Older year to get data.
#' @param user character. From HMD/HFD.
#' @param pass character. From HMD/HFD.
#'
#' @return A list wiith survival, fertility and poopulation age specific matrixes, with calendar year as colnames.
#' @export

get_HMDHFD <- function(country = "SWE",
                         min_year = 1900,
                         max_year = 2018,
                         user = NULL,
                         pass = NULL){

  # source HMD HFD - use SWE now-----------------------------------------------------------------
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

  return(list(P=P, asfr=asfr, N=N))
}

#' Get time serie matrix data from HMD/HFD

#' @description Wrapper function to get data of female survival, fertlity and population
#' of selected country on selected period.

#' @param country numeric. Country code from rom HMD/HFD.
#' @param max_year numeric. Latest year to get data.
#' @param min_year integer. Older year to get data.
#' @param user_HMD character. From HMD.
#' @param user_HFD character. From HFD.
#' @param pass_HMD character. From HMD.
#' @param pass_HFD character. From HFD.
#' @param OAG numeric. Open age group to standarize output.
#' @return A list wiith female survival probability, survival function, fertility and poopulation age specific matrixes, with calendar year as colnames.
#' @export

get_HMDHFD <- function(country = "SWE",
                         min_year = 1900,
                         max_year = 2018,
                         user_HMD = NULL,
                         pass_HMD = NULL,
                         user_HFD = NULL,
                         pass_HFD = NULL,
                         OAG = 100){

  if(any(c(is.null(user_HMD), is.null(user_HFD), is.null(pass_HMD), is.null(pass_HFD)))){
    stop("The function needs HMD and HMF access.")
  }

  # source HMD HFD -----------------------------------------------------------------
  pop <- HMDHFDplus::readHMDweb(CNTRY = country, "Population", user_HMD, pass_HMD, fixup = TRUE) %>%
          dplyr::select(Year, Age, N = Female1)%>%
          dplyr::filter(Year >= min_year, Year <= max_year)
  lt <- HMDHFDplus::readHMDweb(country, "fltper_1x1", user_HMD, pass_HMD, fixup = TRUE) %>%
          dplyr::filter(Year >= min_year, Year <= max_year)
  asfr <- HMDHFDplus::readHFDweb(country, "asfrRR", user_HFD, pass_HFD, fixup = TRUE)%>%
          dplyr::filter(Year >= min_year, Year <= max_year)

  # list of yearly Leslie matrix ---------------------------------------------------

  age = 0:OAG
  ages = length(age)
  w = last(age)
  last_year = max(lt$Year)
  years = min_year:last_year

  # survival probability
  px <- lt %>%
    dplyr::filter(Age<=OAG) %>%
    dplyr::mutate(px = 1 - qx,
           px = ifelse(Age==OAG, 0, px)) %>%
    dplyr::select(Year, Age, px) %>%
    tidyr::pivot_wider(names_from = "Year", values_from = "px") %>%
    dplyr::select(-Age) %>%
    as.matrix()
  rownames(px) = 0:OAG

  # survival function
  Lx <- lt %>%
    dplyr::filter(Age<=OAG) %>%
    dplyr::mutate(Lx = ifelse(Age==OAG, Tx, Lx)) %>%
    dplyr::select(Year, Age, Lx) %>%
    tidyr::pivot_wider(names_from = "Year", values_from = "Lx") %>%
    dplyr::select(-Age) %>%
    as.matrix()

  Sx <- rbind(Lx[c(-1,-ages),]/Lx[-c(w:ages),],
             Lx[ages,]/(Lx[w,]+Lx[ages,]),
             Lx[ages,]/(Lx[w,]+Lx[ages,]))
  rownames(Sx) = 0:w

  # fertility
  fx <- asfr %>%
    dplyr::filter(Year >= min_year) %>%
    dplyr::select(-OpenInterval) %>%
    rbind(
      expand.grid(Year = years,
                  Age = c(0:(min(asfr$Age)-1),(max(asfr$Age)+1):OAG),
                  ASFR = 0)) %>%
    dplyr::arrange(Year, Age) %>%
    tidyr::spread(Year, ASFR) %>%
    dplyr::select(-Age) %>%
    as.matrix()
  rownames(fx) = 0:OAG

  # population
  Nx <- pop %>%
    dplyr::mutate(Age = ifelse(Age>OAG, OAG, Age)) %>%
    dplyr::group_by(Year, Age) %>% summarise(N = sum(N)) %>%
    dplyr::filter(Age<=OAG, Year >= min_year) %>%
    dplyr::arrange(Year, Age) %>%
    tidyr::spread(Year, N) %>%
    dplyr::select(-Age) %>%
    as.matrix()
  rownames(Nx) = 0:OAG

  # only return data with values
  if(any(is.na(colSums(Sx)))){
    warning("Asked for data out of HMDHFD range")
    Sx <- Sx[,!is.na(colSums(Sx))]
  }
  if(any(is.na(colSums(fx)))){
    warning("Asked for data out of HMDHFD range")
    fx <- fx[,!is.na(colSums(fx))]
  }
  if(any(is.na(colSums(Nx)))){
    warning("Asked for data out of HMDHFD range")
    Nx <- Nx[,!is.na(colSums(Nx))]
  }

  return(list(px=px,
              Sx=Sx,
              fx=fx,
              Nx=Nx))
}

# save data
  # swe_px <- swe_data$px
  # swe_Sx <- swe_data$Sx
  # swe_asfr <-swe_data$fx
  # swe_pop <- swe_data$Nx
  # save(swe_px, file = "data/swe_px.rda")
  # save(swe_Sx, file = "data/swe_Sx.rda")
  # save(swe_asfr, file = "data/swe_asfr.rda")
  # save(swe_pop, file = "data/swe_pop.rda")

#' UK female fertility from 1965 to 2022
#' @docType data
#' @format
#'   list of age by stage matrices, entries give female fert. List starting 1965 ending 2022.
#' @source
#'   HFD and ONS
"Female_parity_fert_list_UK"

#' UK female parity transitions from 1965 to 2022
#' @docType data
#' @format
#'   list of age by stage matrices, entries give female parity transitions. List starting 1965 ending 2022.
#' @source
#'   HFD and ONS
"Parity_transfers_by_age_list_UK"

#' UK female parity mortality from 1965 to 2022
#' @docType data
#' @format
#'   list of age by stage matrices, entries give female parity mortality List starting 1965 ending 2022.
#' @source
#'   HFD and ONS
"Female_parity_mortality_list_UK"

#' UK male parity mortality from 1965 to 2022
#' @docType data
#' @format
#'   list of age by stage matrices, entries give male parity mortality List starting 1965 ending 2022.
#' @source
#'   HFD and ONS
"Male_parity_mortality_list_UK"

#' UK parity assign parity at birth
#' @docType data
#' @format
#'   list of matrices which redistributes newborns to age-class 1 and parity 0. No time-variation.
#' @source
#'   None
"Redistribution_by_parity_list_UK"

#' Female swedish survival ratios from 1900 to 2015
#' @docType data
#' @format
#'   A matrix with years as cols and ages (0 to 100 as OAG) as rows.
#'
#' @source
#'   HMD/HFD
"swe_Sx"

#' Female swedish survival probabilities from 1900 to 2015
#' @docType data
#' @format
#'   A matrix with years as cols and ages (0 to 100 as OAG) as rows.
#'
#' @source
#'   HMD/HFD
"swe_px"

#' Swedish age-specific fertility rates from 1900 to 2015
#'
#' Swedish age-specific fertility rates from 1900 to 2015
#' @docType data
#' @format
#'   A matrix with years as cols and ages (0 to 100 as OAG) as rows.
#'
#' @source
#'   HMD/HFD
"swe_asfr"

#' Female swedish population from 1900 to 2015
#'
#' Female swedish population from 1900 to 2015
#' @docType data
#' @format
#'   A matrix with years as cols and ages (0 to 100 as OAG) as rows.
#'
#' @source
#'   HMD/HFD
"swe_pop"

#' Historic and projected survival ratios from Sweden used in Caswell (2021)
#'
#' Historic and projected survival ratios from Sweden used in Caswell (2021)
#' @docType data
#' @format
#'   A matrix U with years as cols and ages as rows.
#'
#' @source
#'   Caswell (2019)
"U_caswell_2021"

#' Historic and projected fertility ratios from Sweden used in Caswell (2021)
#'
#' Historic and projected fertility ratios from Sweden used in Caswell (2021)
#' @docType data
#' @format
#'   A matrix f with years as cols and ages as rows.
#'
#' @source
#'   Caswell (2019)
"f_caswell_2021"

#' Historic and projected mother´s age distribution of childbearing from Sweden used in Caswell (2021)
#'
#' Historic and projected mother´s age distribution of childbearing from Sweden used in Caswell (2021)
#' @docType data
#' @format
#'   A matrix pi with years as cols and ages as rows.
#'
#' @source
#'   Caswell (2019)
"pi_caswell_2021"

#' Female Slovakian survival probabilities by parity stage in 1990 (Caswell, 2021)
#'
#' Female Slovakian survival probabilities by parity stage in 1990 (Caswell, 2021)
#' @docType data
#' @format
#'   A matrix of px with stages as cols and ages as rows.
#'
#' @source
#'   Caswell (2021)
"svk_pxs"

#' Female Slovakian fertility rates by parity stage in 1990 (Caswell, 2021)
#'
#' Female Slovakian fertility rates by parity stage in 1990 (Caswell, 2021)
#' @docType data
#' @format
#'   A matrix of fx with stages as cols and ages as rows.
#'
#' @source
#'   Caswell (2021)
"svk_fxs"

#' Age where assign offspring of individuals in each partity stage (Caswell, 2021). All to zero age in this case.
#'
#' Age where assign offspring of individuals in each partity stage (Caswell, 2021). All to zero age in this case.
#' @docType data
#' @format
#'   A matrix of ones in ages where assign offspring individuals, with stages as cols and ages as rows.
#'
#' @source
#'   Caswell (2021)
"svk_Hxs"

#' Probability of transition among parity stage for Slovakia in 1990, for each age, conditional on survival (Caswell, 2021).
#'
#' Probability of transition among parity stage for Slovakia in 1990, for each age, conditional on survival (Caswell, 2021).
#' @docType data
#' @format
#' A list of column-stochastic matrix with probabilities of transition among parity stage, for each age, conditional on survival.
#'
#' @source
#'   Caswell (2021)
"svk_Uxs"

#' Output for Slovakia 1990 in Caswell (2020).
#'
#' Output for Slovakia 1990 in Caswell (2020).
#' @docType data
#' @format
#' A list with specific kin types age-stage matrix
#'
#' @source
#'   Caswell (2021)
"kin_svk1990_caswell2020"

#' Fertility for France (2012) by sex in Caswell (2022).
#'
#' Fertility for France (2012) by sex in Caswell (2022).
#' @docType data
#' @format
#' A data.frame with age specific fertility rates by age and sex.
#'
#' @source
#'   Caswell (2022)
"fra_asfr_sex"

#' Survival probability for France (2012) by sex in Caswell (2022).
#'
#' Survival probability for France (2012) by sex in Caswell (2022).
#' @docType data
#' @format
#' A data.frame with survival probabilities by age and sex.
#'
#' @source
#'   Caswell (2022)
"fra_surv_sex"

#' DemoKin codes, Caswell (2020) codes, and useful labels.
#'
#' DemoKin codes, Caswell (2020) codes, and useful labels.
#' @docType data
#' @format
#' A data.frame with codes and labels for distinction between kin types.

"demokin_codes"

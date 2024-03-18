
rm(list = ls())

## Compare output

library(DemoKin)

source("./R/kin_time_invariant_2sex_cod.R")
source("./R/kin_time_variant_2sex_cod.R")

source("./R/kin_time_invariant_2sex.R")
source("./R/kin_time_variant_2sex.R")

## Example from Vignette 2 sex but few years for speed

years_all <- ncol(swe_px)
years_test <- 1:20
ages <- nrow(swe_px)
pf <- swe_px[,years_test]
pm <- swe_px[,years_test] ^ 1.5 # artificial perturbation for this example
ff <- swe_asfr[, years_test]
fm <- rbind(matrix(0, 5, max(years_test)),
                           swe_asfr[-((ages-4):ages), years_test]) * 1.05 # artificial perturbation for this example

# Create a fictitious hazard matrix with three causes of death, where each
# year is a list item (Hazard matrix needs to be (causes * ages) for the
# matrix algebra to work well with existing code).
H <- matrix(c(0.5, 1, 2), nrow = 3, ncol = nrow(pf))
Hf <- Hm <- sapply(colnames(pf), function(x) {
  return(H)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
  )

## COMPARISON

start.time <- Sys.time()
no_cod <- kin_time_variant_2sex(pf = pf, pm = pm,
                      ff = ff, fm = fm,
                      sex_focal = "f",
                      birth_female = 1/2.04,
                      pif = NULL, pim = NULL,
                      nf = NULL, nm = NULL,
                      output_cohort = NULL, output_period = NULL, output_kin = NULL,
                      list_output = FALSE)
end.time <- Sys.time()
time.taken.no.cod <- end.time - start.time


start.time <- Sys.time()
cod <- kin_time_variant_2sex_cod(pf = pf, pm = pm,
                                 ff = ff, fm = fm,
                                 Hf = Hf, Hm = Hm,
                                 sex_focal = "f",
                                 birth_female = 1/2.04,
                                 pif = NULL, pim = NULL,
                                 nf = NULL, nm = NULL,
                                 output_cohort = NULL, output_period = NULL, output_kin = NULL,
                                 list_output = FALSE)
end.time <- Sys.time()
time.taken.cod <- end.time - start.time

no_cod
cod |>
  mutate(
    dead = deadcause1 + deadcause2 + deadcause3
    )


no_cod |>
  filter(
    kin == "gm",
    age_focal %in% 30:35,
    age_kin > 60
  )
cod |>
  mutate(
    dead = deadcause1 + deadcause2 + deadcause3
    ) |>
  filter(
    kin == "gm",
    age_focal %in% 30:35,
    age_kin > 60
  )

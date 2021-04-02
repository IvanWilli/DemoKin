
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DemoKin

This package uses matrix demographic methods to estimate kin counts and
the age distribution of relatives in stable and non-stable populations.
The package is an implementation of Caswell (2019) and draws on previous
theoretical development by Goodman, Keyfitz and Pullum (1974). `DemoKin`
is giving its first steps, so please contact us or submit a pull request
if you have any suggestion.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("IvanWilli/DemoKin")
```

## Example

Consider an average Swedish woman aged 50 in year 2015. How many
relatives can this woman expect to have? The function `kins` can be used
to obtain the expected number of Ego’s relatives for the following types
of relative:

| Code | Relative                                   |
| :--- | :----------------------------------------- |
| coa  | Cousins (through aunt older than mother)   |
| cya  | Cousins (through aunt younger than mother) |
| d    | Daughter                                   |
| gd   | Grand-daughter                             |
| ggm  | Great-grandmother                          |
| gm   | Grandmother                                |
| m    | Mother                                     |
| nos  | Nieces through older sister                |
| nys  | Nieces through younger sister              |
| oa   | Aunt older than mother                     |
| ya   | Aunt younger than mother                   |
| os   | Older sister                               |
| ys   | Younger sister                             |

For this example, we assume demographic stability (i.e., we assume that
the womans and her relatives experienced the mortality and fertility
rates from 2015 at each age throughout their life):

``` r
library(DemoKin)
swe50_2015_stable <- kins(ego_age = 50, year = 2015,
                             P = swe_surv, asfr = swe_asfr,
                             stable = TRUE)

plot_diagram(swe50_2015_stable[["kins_total"]],ego_age = 50)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

For more details, including an extension to non-stable populations, see
`vignette("Use")`.

## Citation

Williams, Iván and Diego Alburez-Gutierrez. (2021) DemoKin: An R package
to implement matrix kinship models in stable and non-stable populations.
URL: <https://github.com/IvanWilli/DemoKin>.

## References

Caswell, H. (2019). The formal demography of kinship: A matrix
formulation. Demographic Research 41:679–712.
<doi:10.4054/DemRes.2019.41.24>.

Caswell, H. (2020). The formal demography of kinship II: Multistate
models, parity, and sibship. Demographic Research 42:1097–1146.
<doi:10.4054/DemRes.2020.42.38>.

Goodman, L.A., Keyfitz, N., and Pullum, T.W. (1974). Family formation
and the frequency of various kinship relationships. Theoretical
Population Biology 5(1):1–27. <doi:10.1016/0040-5809(74)90049-5>.

## Next steps:

  - Improve performance of `kins_non_stable` function.
  - Give an option to forecast mortality and fertility to locate Ego and
    relatives in the future.
  - Add more functionalities to the diagram, like colors by kin degree
    and box size weighted by kin amount.
  - Add cumulative deaths by kin in the non-stable case, and also stage
    properties as in Caswell (2020).

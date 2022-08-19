
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DemoKin

`DemoKin` uses matrix demographic methods to compute expected (average)
kin counts from demographic rates under a range of scenarios and
assumptions. The package is an R-language implementation of Caswell
(2019) and Caswell and Song (2021). It draws on previous theoretical
development by Goodman, Keyfitz and Pullum (1974).

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("IvanWilli/DemoKin")
```

## Example

Consider an average Swedish woman in year 2015 (i.e., ‘Focal’).
Relatives for the `selected_kin` argument are identified by a unique
code. Note that the relationship codes used in `DemoKin` differ from
those in Caswell (2019). The equivalence between the two set of codes is
given in the following table:

| DemoKin | Caswell | Label                      |
|:--------|:--------|:---------------------------|
| coa     | t       | Cousins from older aunt    |
| cya     | v       | Cousins from younger aunt  |
| d       | a       | Daughter                   |
| gd      | b       | Grand-daughter             |
| ggd     | c       | Great-grand-daughter       |
| ggm     | h       | Great-grandmother          |
| gm      | g       | Grandmother                |
| m       | d       | Mother                     |
| nos     | p       | Nieces from older sister   |
| nys     | q       | Nieces from younger sister |
| oa      | r       | Aunt older than mother     |
| ya      | s       | Aunt younger than mother   |
| os      | m       | Older sister               |
| ys      | n       | Younger sister             |

Equivalence between relative codes between DemoKin and Caswell (2019).

Let’s show a quick example. We assume a female closed population in
which everyone experiences the Swedish 2015 mortality and fertility
rates at each age throughout their life. We then ask:

> How can we characterize the kinship network of an average member of
> the population (‘Focal’)?

For this exercise, we’ll use the Swedish data pre-loaded with `DemoKin`.

``` r
library(DemoKin)
swe_surv_2015 <- swe_px[,"2015"]
swe_asfr_2015 <- swe_asfr[,"2015"]
swe_2015 <- kin(U = swe_surv_2015, f = swe_asfr_2015, time_invariant = TRUE)
```

*px* is the survival probability by age from a life table and *f* are
the age specific fertility ratios by age (see `?kin` for details).

For more details, including an extension to time varying populations
rates and relative´s death distribution, see `vignette("Use")`. Note
that if the vignette does not load, you may need to install the package
as `devtools::install_github("IvanWilli/DemoKin", build_vignettes = T)`.

## Citation

Williams, Iván; Alburez-Gutierrez, Diego; Song, Xi; and Hal Caswell.
(2021) DemoKin: An R package to estimate kinship networks in stable and
non-stable populations. URL: <https://github.com/IvanWilli/DemoKin>.

## References

Caswell, H. (2019). The formal demography of kinship: A matrix
formulation. Demographic Research 41:679–712.
<doi:10.4054/DemRes.2019.41.24>.

Caswell, Hal and Xi Song. 2021. “The Formal Demography of Kinship. III.
Kinship Dynamics with Time-Varying Demographic Rates.” Demographic
Research 45: 517–46. <doi:10.4054/DemRes.2021.45.16>.

Goodman, L.A., Keyfitz, N., and Pullum, T.W. (1974). Family formation
and the frequency of various kinship relationships. Theoretical
Population Biology 5(1):1–27. <doi:10.1016/0040-5809(74)90049-5>.

## Next steps:

1.  Implement multi-stage and two-sex matrix kinship models
2.  Create a hex logo for package
3.  Improve kinship diagram visualization
4.  Improve documentation and vignette of package

## Get involved!

`DemoKin` is giving its first steps. If you’re interested in
contributing, please get in touch, create an issue, or submit a pull
request. We look forward to hearing from you!

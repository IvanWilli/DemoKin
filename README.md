
# DemoKin

<div class="columns">

<div class="column" width="60%">

`DemoKin` uses matrix demographic methods to compute expected (average)
kin counts from demographic rates under a range of scenarios and
assumptions. The package is an R-language implementation of Caswell
(2019, 2020, 2022), and Caswell and Song (2021). It draws on previous
theoretical development by Goodman, Keyfitz and Pullum (1974).

</div>

<div class="column" width="40%">

<img src="man/figures/DemoKin-Logo.png" align="right" width="200" />

</div>

</div>

## Installation

Download the stable version [from
CRAN](https://cran.r-project.org/web/packages/DemoKin/):

``` r
install.packages("DemoKin")
```

Or you can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("IvanWilli/DemoKin")
```

## Usage

Consider an average Swedish woman called ‘Focal.’ For this exercise, we
assume a female closed population in which everyone experiences the
Swedish 2015 mortality and fertility rates at each age throughout their
life; i.e., the ‘time-invariant’ assumption in Caswell (2019).

We then ask:

> What is the expected number of relatives of Focal over her life
> course?

Let’s explore this using the Swedish data already included with
`DemoKin`.

``` r
library(DemoKin)
swe_surv_2015 <- swe_px[,"2015"]
swe_asfr_2015 <- swe_asfr[,"2015"]
swe_2015 <- kin(p = swe_surv_2015, f = swe_asfr_2015, time_invariant = TRUE)
```

*p* is the survival probability by age from a life table and *f* are the
age specific fertility ratios by age (see `?kin` for details).

Now, we can visualize the implied kin counts (i.e., the average number
of living kin) of Focal at age 35 using a network or ‘Keyfitz’ kinship
diagram with the function `plot_diagram`:

``` r
# We need to reformat the data a little bit
kin_total <- swe_2015$kin_summary
# Keep only data for Focal's age 35
kin_total <- kin_total[kin_total$age_focal == 35 , c("kin", "count_living")]
names(kin_total) <- c("kin", "count")
plot_diagram(kin_total, rounding = 2)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Relatives are identified by a unique code:

| DemoKin | Labels_female               | Labels_male                   | Labels_2sex                       |
|:--------|:----------------------------|:------------------------------|:----------------------------------|
| coa     | Cousins from older aunts    | Cousins from older uncles     | Cousins from older aunts/uncles   |
| cya     | Cousins from younger aunts  | Cousins from younger uncles   | Cousins from younger aunts/uncles |
| c       | Cousins                     | Cousins                       | Cousins                           |
| d       | Daughters                   | Brothers                      | Siblings                          |
| gd      | Grand-daughters             | Grand-sons                    | Grand-childrens                   |
| ggd     | Great-grand-daughters       | Great-grand-sons              | Great-grand-childrens             |
| ggm     | Great-grandmothers          | Great-grandfathers            | Great-grandfparents               |
| gm      | Grandmothers                | Grandfathers                  | Grandparents                      |
| m       | Mother                      | Father                        | Parents                           |
| nos     | Nieces from older sisters   | Nephews from older brothers   | Niblings from older siblings      |
| nys     | Nieces from younger sisters | Nephews from younger brothers | Niblings from younger siblings    |
| n       | Nieces                      | Nephews                       | Niblings                          |
| oa      | Aunts older than mother     | Uncles older than fathers     | Aunts/Uncles older than parents   |
| ya      | Aunts younger than mother   | Uncles younger than father    | Aunts/Uncles younger than parents |
| a       | Aunts                       | Uncles                        | Aunts/Uncles                      |
| os      | Older sisters               | Older brothers                | Older siblings                    |
| ys      | Younger sisters             | Younger brothers              | Younger siblings                  |
| s       | Sisters                     | Brothers                      | Siblings                          |

## Vignette

For more details, including an extension to time-variant rates, deceased
kin, and multi-state models in a one-sex framework, see the
[Reference_OneSex](https://cran.r-project.org/web/packages/DemoKin/vignettes/Reference_OneSex.html)
vignette; also accessible from DemoKin:
`vignette("Reference_OneSex", package = "DemoKin")`. For two-sex models,
see the
[Reference_TwoSex](https://cran.r-project.org/web/packages/DemoKin/vignettes/Reference_TwoSex.html)
vignette; also accessible from DemoKin:
`vignette("Reference_TwoSex", package = "DemoKin")`. If the vignette
does not load, you may need to install the package as
`devtools::install_github("IvanWilli/DemoKin", build_vignettes = T)`.

## Citation

Williams, Iván; Alburez-Gutierrez, Diego; Song, Xi; and Hal Caswell.
(2021) DemoKin: An R package to implement demographic matrix kinship
models. URL: <https://github.com/IvanWilli/DemoKin>.

## Acknowledgments

We thank Silvia Leek from the Max Planck Institute for Demographic
Research for designing the DemoKin logo. The logo includes elements that
have been taken or adapted [from this
file](https://commons.wikimedia.org/wiki/File:Escudo_de_la_Orden_de_San_Jer%C3%B3nimo.svg),
originally by Ansunando, [CC BY-SA
4.0](https://creativecommons.org/licenses/by-sa/4.0) via Wikimedia
Commons. Sha Jiang provided useful comments for improving the package.

## Get involved!

`DemoKin` is under constant development. If you’re interested in
contributing, please get in touch, create an issue, or submit a pull
request. We look forward to hearing from you!

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-caswell_formal_2019" class="csl-entry">

Caswell, Hal. 2019. “The Formal Demography of Kinship: A Matrix
Formulation.” *Demographic Research* 41 (September): 679–712.
<https://doi.org/10.4054/DemRes.2019.41.24>.

</div>

<div id="ref-caswell_formal_2020" class="csl-entry">

———. 2020. “The Formal Demography of Kinship II: Multistate Models,
Parity, and Sibship.” *Demographic Research* 42 (June): 1097–1146.
<https://doi.org/10.4054/DemRes.2020.42.38>.

</div>

<div id="ref-caswell_formal_2022" class="csl-entry">

———. 2022. “The Formal Demography of Kinship IV: Two-Sex Models and
Their Approximations.” *Demographic Research* 47 (September): 359–96.
<https://doi.org/10.4054/DemRes.2022.47.13>.

</div>

<div id="ref-caswell_formal_2021" class="csl-entry">

Caswell, Hal, and Xi Song. 2021. “The Formal Demography of Kinship III:
Kinship Dynamics with Time-Varying Demographic Rates.” *Demographic
Research* 45 (August): 517–46.
<https://doi.org/10.4054/DemRes.2021.45.16>.

</div>

<div id="ref-goodman_family_1974" class="csl-entry">

Goodman, Leo A, Nathan Keyfitz, and Thomas W. Pullum. 1974. “Family
Formation and the Frequency of Various Kinship Relationships.”
*Theoretical Population Biology*, 27.
<https://doi.org/10.1016/0040-5809(74)90049-5>.

</div>

</div>

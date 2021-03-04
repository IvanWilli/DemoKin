
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DemoKin

This is a package for estimating kin count in from a demographic point
of view (Goodman et. al, 1978).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("IvanWilli/DemoKin")
```

## Example

Thinking in a swedish female aged 50 in 2015, the expected kins are
this:

``` r
swe50_2015_stable <- kins(ego_age = 50, year = 2015,
                             P = swe_surv, asfr = swe_asfr,
                             stable = TRUE)
plot_diagramm(swe50_2015_stable[["kins_total"]],ego_age = 50)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Next steps are:

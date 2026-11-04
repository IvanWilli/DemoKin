# DemoKin 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Change stable/non-stable references to time varying/non-varying rates.
* Add multi-state process.

# DemoKin 1.0.1

* Submitted to CRAN
* Death counts are placed in the age where Focal experience the death.
* Aggregated kin types are allowed (`s` for older and younger sisters, for example).

# DemoKin 1.0.2

* Added the Two sex multi-state time-varying framework. 

# DemoKin 1.0.3

* Website
* Archive CRAN
* Typos

# DemoKin 1.1.0

* New argument `fast` (default `FALSE`) in `kin()`, `kin2sex()`, `kin_multi_stage()`,
  `kin_multi_stage_time_variant_2sex()`, the six `kin_time_*()` engines and the three
  `timevarying_kin*()` one-step projectors. With `fast = TRUE` the projections run in
  compiled C++ code and the output tables are built directly; the result is identical to
  `fast = FALSE` (same values, columns, row order and classes), typically 5 to 200 times
  faster. The R implementation run by default is unchanged. Equivalence is checked in the
  test suite and on a 46-scenario harness (`dev/equivalence/`).
* The package now contains compiled code (`Rcpp` in Imports and LinkingTo), so installing
  from source needs a C++ compiler (Rtools on Windows, Xcode command line tools on macOS).

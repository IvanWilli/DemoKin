% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kin_time_invariant.R
\name{kin_time_invariant}
\alias{kin_time_invariant}
\title{Estimate kin counts in a time invariant framework for one-sex model (matrilineal/patrilineal)}
\usage{
kin_time_invariant(
  p = NULL,
  f = NULL,
  birth_female = 1/2.04,
  pi = NULL,
  output_kin = NULL,
  list_output = FALSE
)
}
\arguments{
\item{p}{numeric. A vector of survival probabilities with same length as ages.}

\item{f}{numeric. A vector of age-specific fertility rates with same length as ages.}

\item{birth_female}{numeric. Female portion at birth.}

\item{pi}{numeric. For using some specific non-stable age distribution of childbearing (same length as ages). Default \code{NULL}.}

\item{output_kin}{character. kin to return. For example "m" for mother, "d" for daughter. See \code{vignette} for all kin types.}

\item{list_output}{logical. Results as a list with \code{output_kin} elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default \code{FALSE}}
}
\value{
A data frame with focal´s age, related ages and type of kin
(for example \code{d} is daughter, \code{oa} is older aunts, etc.), alive and death. If \code{list_output = TRUE} then this is a list.
}
\description{
Mtrix implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019).
}

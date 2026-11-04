# Fast (C++) implementation of kin_time_invariant_2sex(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_time_invariant_2sex.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).

# Estimate kin counts in a time invariant framework for two-sex model.

# @description Two-sex matrix framework for kin count estimates.This produces kin counts grouped by kin, age and sex of
# each relatives at each Focal´s age. For example, male cousins from aunts and uncles from different sibling's parents
# are grouped in one male count of cousins.
# @details See Caswell (2022) for details on formulas.
# @param pf numeric. A vector of survival probabilities for females with same length as ages.
# @param ff numeric. A vector of age-specific fertility rates for females with same length as ages.
# @param pm numeric. A vector of survival probabilities for males with same length as ages.
# @param fm numeric. A vector of age-specific fertility rates for males with same length as ages.
# @param sex_focal character. "f" for female or "m" for male.
# @param birth_female numeric. Female portion at birth.
# @param pif numeric. For using some specific non-stable age distribution of childbearing for mothers (same length as ages). Default `NULL`.
# @param pim numeric. For using some specific non-stable age distribution of childbearing for fathers (same length as ages). Default `NULL`.
# @param output_kin character. kin to return, considering matrilineal names. For example "m" for parents, "d" for children, etc. See the `vignette` for all kin types.
# @param list_output logical. Results as a list with `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
#
# @return A data frame with focal´s age, related ages and type of kin
# (for example `d` is children, `oa` is older aunts/uncles, etc.), sex, alive and death. If `list_output = TRUE` then this is a list.

kin_time_invariant_2sex_fast <- function(pf = NULL, pm = NULL,
                                    ff = NULL, fm = NULL,
                                    sex_focal = "f",
                                    birth_female = 1/2.04,
                                    pif = NULL, pim = NULL,
                                    output_kin = NULL,
                                    list_output = FALSE){

  # same input length
  if(!all(length(pf)==length(pm), length(pf)==length(ff), length(pf)==length(fm))) stop("Lengths of p's and f's should be the same")

  # make matrix transition from vectors. Include death counts with matrix M
  age = 0:(length(pf)-1)
  ages = length(age)
  agess = ages * 2
  Uf = Um = Ff = Fm = Gt = zeros = matrix(0, nrow=ages, ncol=ages)
  Uf[row(Uf)-1 == col(Uf)] <- pf[-ages]
  Uf[ages, ages] = pf[ages]
  Um[row(Um)-1 == col(Um)] <- pm[-ages]
  Um[ages, ages] = pm[ages]
  Mm <- diag(1-pm)
  Mf <- diag(1-pf)
  Ut <- rbind(  # dense blocks, as as.matrix() of the Matrix::bdiag() version
        cbind(bdiag_dense(Uf, Um), bdiag_dense(zeros, zeros)),
        cbind(bdiag_dense(Mf, Mm), bdiag_dense(zeros, zeros)))
  Ff[1,] = ff
  Fm[1,] = fm
  Ft <- Ft_star <- matrix(0, agess*2, agess*2)
  Ft[1:agess,1:agess] <- rbind(cbind(birth_female * Ff, birth_female * Fm),
                                     cbind((1-birth_female) * Ff, (1-birth_female) * Fm))

  # mother and father do not reproduce independently to produce focal´s siblings. Assign to mother
  Ft_star[1:agess,1:ages] <- rbind(birth_female * Ff, (1-birth_female) * Ff)

  # parents age distribution under stable assumption in case no input
  if(is.null(pim) | is.null(pif)){
    A = bdiag_dense(Uf, Um) + Ft_star[1:agess,1:agess]  # dense, as Matrix returned
    A_decomp = eigen(A)
    lambda = as.double(A_decomp$values[1])
    w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
    wf = w[1:ages]
    wm = w[(ages+1):(2*ages)]
    pif = wf * ff / sum(wf * ff)
    pim = wm * fm / sum(wm * fm)
  }

  # focal´s trip (C++ kernel)
  kin_list <- project_kin_2sex(Ut, Ft, Ft_star, c(pif, pim), sex_focal, ages)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list[names(kin_list) %in% output_kin]
  }

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- build_kin_df_ti(kin_list, age, layout_2sex_ti(age))
  }

  return(out)
}

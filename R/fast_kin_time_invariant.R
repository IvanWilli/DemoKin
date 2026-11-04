# Fast (C++) implementation of kin_time_invariant(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_time_invariant.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).

# Estimate kin counts in a time invariant framework for one-sex model (matrilineal/patrilineal)

# @description Mtrix implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019).

# @param p numeric. A vector of survival probabilities with same length as ages.
# @param f numeric. A vector of age-specific fertility rates with same length as ages.
# @param birth_female numeric. Female portion at birth.
# @param pi numeric. For using some specific non-stable age distribution of childbearing (same length as ages). Default `NULL`.
# @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See `vignette` for all kin types.
# @param list_output logical. Results as a list with `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
#
# @return A data frame with focal´s age, related ages and type of kin
# (for example `d` is daughter, `oa` is older aunts, etc.), alive and death. If `list_output = TRUE` then this is a list.

kin_time_invariant_fast <- function(p = NULL, f = NULL,
                        birth_female = 1/2.04,
                        pi = NULL,
                        output_kin = NULL,
                        list_output = FALSE){

  # make matrix transition from vectors
  age = 0:(length(p)-1)
  ages = length(age)
  Ut = Mt = zeros = matrix(0, nrow=ages, ncol=ages)
  Ut[row(Ut)-1 == col(Ut)] <- p[-ages]
  Ut[ages, ages] = p[ages]
  diag(Mt) = 1 - p
  Ut = rbind(cbind(Ut,zeros),
             cbind(Mt,zeros))
  ft = matrix(0, nrow=ages*2, ncol=ages*2)
  ft[1,1:ages] = f * birth_female

  # stable age distribution in case no pi is given
  if(is.null(pi)){
    A = Ut[1:ages,1:ages] + ft[1:ages,1:ages]
    A_decomp = eigen(A)
    lambda = as.double(A_decomp$values[1])
    w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
    pi = w*A[1,]/sum(w*A[1,])
  }

  # focal´s trip (C++ kernel, same recursions as Caswell 2019)
  kin_list <- project_kin_1sex(Ut, ft, pi, ages)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list[names(kin_list) %in% output_kin]
  }

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- build_kin_df_ti(kin_list, age, layout_1sex(age))
  }
  return(out)
}

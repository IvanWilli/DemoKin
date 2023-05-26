#' Estimate kin counts in a time invariant framework for one-sex model (matrilineal/patrilineal)

#' @description Mtrix implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019).

#' @param p numeric. A vector of survival probabilities with same length as ages.
#' @param f numeric. A vector of age-specific fertility rates with same length as ages.
#' @param birth_female numeric. Female portion at birth.
#' @param pi numeric. For using some specific non-stable age distribution of childbearing (same length as ages). Default `NULL`.
#' @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See `vignette` for all kin types.
#' @param list_output logical. Results as a list with `output_kin` elements, with focal´s age in columns and kin ages in rows (2 * ages, last chunk of ages for death experience). Default `FALSE`
#'
#' @return A data frame with focal´s age, related ages and type of kin
#' (for example `d` is daughter, `oa` is older aunts, etc.), alive and death. If `list_output = TRUE` then this is a list.
#' @export

kin_time_invariant <- function(p = NULL, f = NULL,
                        birth_female = 1/2.04,
                        pi = NULL,
                        output_kin = NULL,
                        list_output = FALSE){

  # global vars
  .<-alive<-age_kin<-alive<-age_focal<-count<-NULL

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

  # identity matrix
  e = matrix(0, ages * 2, ages * 2)
  diag(e[1:ages,1:ages]) = 1

  # initial vectors
  d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0, ages * 2, ages)

  # focal´s trip
  m[,1] = c(pi, rep(0,ages))
  for(i in 1:(ages-1)){
    d[,i+1]   = Ut %*% d[,i] + ft %*% e[,i]
    gd[,i+1]  = Ut %*% gd[,i] + ft %*% d[,i]
    ggd[,i+1]  = Ut %*% ggd[,i] + ft %*% gd[,i]
    m[,i+1]   = Ut %*% m[,i]
    ys[,i+1]  = Ut %*% ys[,i] + ft %*% m[,i]
    nys[,i+1] = Ut %*% nys[,i] + ft %*% ys[,i]
  }
  gm[1:ages,1] = m[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    gm[,i+1]  = Ut %*% gm[,i]
  }
  ggm[1:ages,1] = gm[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    ggm[,i+1]  = Ut %*% ggm[,i]
  }
  os[1:ages,1]  = d[1:ages,] %*% pi
  nos[1:ages,1] = gd[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    os[,i+1]  = Ut %*% os[,i]
    nos[,i+1] = Ut %*% nos[,i] + ft %*% os[,i]
  }
  oa[1:ages,1]  = os[1:ages,] %*% pi
  ya[1:ages,1]  = ys[1:ages,] %*% pi
  coa[1:ages,1] = nos[1:ages,] %*% pi
  cya[1:ages,1] = nys[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    oa[,i+1]  = Ut %*% oa[,i]
    ya[,i+1]  = Ut %*% ya[,i]  + ft %*% gm[,i]
    coa[,i+1] = Ut %*% coa[,i] + ft %*% oa[,i]
    cya[,i+1] = Ut %*% cya[,i] + ft %*% ya[,i]
  }

  # get results
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list %>% purrr::keep(names(.) %in% output_kin)
  }

  # reshape as data.frame
  kin <- purrr::map2(kin_list, names(kin_list),
               function(x,y){
                    # reassign deaths to Focal experienced age
                    x[(ages+1):(ages*2),1:(ages-1)] <- x[(ages+1):(ages*2),2:ages]
                    x[(ages+1):(ages*2),ages] <- 0
                    out <- as.data.frame(x)
                    colnames(out) <- age
                    out %>%
                      dplyr::mutate(kin = y,
                                   age_kin = rep(age,2),
                                   alive = c(rep("living",ages), rep("dead",ages))) %>%
                      tidyr::pivot_longer(c(-age_kin, -kin, -alive), names_to = "age_focal", values_to = "count") %>%
                      dplyr::mutate(age_focal = as.integer(age_focal)) %>%
                      tidyr::pivot_wider(names_from = alive, values_from = count)
                    }
               ) %>%
              purrr::reduce(rbind)

  # results as list?
  if(list_output) {
    out <- kin_list
    }else{
      out <- kin
    }
  return(out)
}

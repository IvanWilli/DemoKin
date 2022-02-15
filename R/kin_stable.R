#' Estimate kin counts in a stable framework

#' @description Implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019).

#' @param U numeric. A vector of survival ratios.
#' @param f numeric. A vector of age-specific fertility rates.
#' @param birth_female numeric. Female portion at birth.
#' @param pi numeric. For using the non-stable distribution of childbearing. Default `NULL`.
#' @param selected_kin character. kin to return: "m" for mother, "d" for daughter,...
#' @param pi_stable logical. Want mean age at childbearing as a result too. Default `FALSE`
#' @param list_output logical. Want kin results as a list. Default `FALSE`
#'
#' @return A data frame with ego´s age, related ages and type of kin
#' (for example `d` is daughter, `oa` is older aunts, etc.), alive and death.
#' @export

kin_stable <- function(U = NULL, f = NULL,
                        birth_female = 1/2.04,
                        pi = NULL,
                        selected_kin = NULL,
                        pi_stable = FALSE,
                        list_output = FALSE){

  # make matrix transition from vectors
  age = 0:(length(U)-1)
  ages = length(age)
  Ut = Mt = zeros = Dcum = matrix(0, nrow=ages, ncol=ages)
  Ut[row(Ut)-1 == col(Ut)] <- U[-ages]
  Ut[ages, ages] = U[ages]
  diag(Mt) = 1 - U
  Ut = rbind(cbind(Ut,zeros),
             cbind(Mt,Dcum))
  ft = matrix(0, nrow=ages*2, ncol=ages*2)

  # Caswell's assumption
  ft[1,1:ages] = f * U * birth_female

  # stable age distr
  if(is.null(pi)){
    A = Ut[1:ages,1:ages] + ft[1:ages,1:ages]
    A_decomp = eigen(A)
    lambda = as.double(A_decomp$values[1])
    w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
    pi = w*A[1,]/sum(w*A[1,])
  }

  # identity
  e = matrix(0, ages * 2, ages * 2)
  diag(e[1:ages,1:ages]) = 1

  # initial vectors
  d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0, ages * 2, ages)

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
  if(!is.null(selected_kin)){
    kin_list <- kin_list %>% keep(names(.) %in% selected_kin)
  }

  kin <- map2(kin_list, names(kin_list),
               function(x,y){
                    out = as.data.frame(x)
                    colnames(out) = age
                    out %>%
                      mutate(kin = y,
                             age_kin = rep(age,2),
                             alive = c(rep("yes",ages), rep("no",ages))) %>%
                      gather(age_ego,count,-age_kin, -kin, -alive) %>%
                      mutate(age_ego = as.integer(age_ego))
                    }
               ) %>%
              reduce(rbind)

  if(pi_stable){
    out <- list(kin=kin, pi_stable=pi)
  }else{
    out <- kin
  }

  if(list_output) out <- kin_list

  return(out)
}

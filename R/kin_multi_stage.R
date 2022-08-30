#' Estimate kin counts by age and stage in a time invariant framework

#' @description Implementation of Goodman-Keyfitz-Pullum equations adapted by Caswell (2019). Stages are implied in length of input lists.

#' @param U list. age elemnts with column-stochastic transition matrix with dimension for the state space, conditional on survival.
#' @param f matrix. state-specific fertility (age in rows and states in columns).
#' @param D matrix. survival probabilities by state (age in rows and states in columns)
#' @param H matrix. assigns the offspring of individuals in some stage to the appropriate age class with 1 (age in rows and states in columns).
#' @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See the `vignette` for all kin types.
#' @param list_output logical. Results as a list. Default `FALSE`.

#' @return A data frame with focal´s age, related ages and type of kin
#' (for example `d` is daughter, `oa` is older aunts, etc.), alive and death, and specific stage. If `list_output = TRUE` then this is a list.
#' @export
#'

kin_multi_stage <- function(U = NULL, f = NULL, D = NULL, H = NULL,
                               birth_female = 1,
                               output_kin = NULL,
                               list_output = FALSE){

  if(!is.list(U)) stop("U must be a list with age length of elements, and stage transitiotn matrix for each one.")

  # stages and ages
  s <- ncol(U[[1]])
  ages <- length(U)
  age <- (1:ages)-1

  # build matrix structure from data.frame input
  H <- purrr::map(colnames(D), function(Y){
    Ht = matrix(0, nrow=ages, ncol=ages)
    Ht[1,] <- 1
    Ht
  })
  D <- purrr::map(colnames(D), function(Y){
      X <- D[,Y]
      Dt = matrix(0, nrow=ages, ncol=ages)
      Dt[row(Dt)-1 == col(Dt)] <- X[-ages]
      Dt[ages, ages] = X[ages]
      Dt
    })
  f <- purrr::map(1:ages, function(Y){
    X <- f[Y,]
    ft = matrix(0, nrow=s, ncol=s)
    ft[1,] <- X
    ft
  })

  # build block matrix
  bbU <- Matrix::bdiag(U)
  bbF <- Matrix::bdiag(f) * birth_female
  bbD <- Matrix::bdiag(D)
  bbH <- Matrix::bdiag(H)

  # rearrange with conmutation matrix
  K <- matrixcalc::commutation.matrix(s, ages)
  Ut <- t(K) %*% bbD %*% K %*% bbU
  ft <-  t(K) %*% bbH %*% K %*% bbF
  Gt <-  Ut%*% MASS::ginv(diag(colSums(as.matrix(Ut))))

  # stable distribution mothers: age x stage
  At <- Ut + ft
  A_decomp <- eigen(At)
  lambda <- as.double(A_decomp$values[1])
  wt <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
  pi <- wt*At[1,]/sum(wt*At[1,])

  # marginal mothers age
  Iom <- diag(1,ages, ages);
  ones <- t(rep(1,s))
  onesom <- t(rep(1,s*ages))
  piage <- kronecker(Iom,ones) %*% pi

  # momarray is an array with pit in each column
  momarray <- pi %*% matrix(1,1,ages)
  Iom = diag(1, ages)
  Is = diag(1, s)
  Isom = diag(1, s*ages)
  zsom = matrix(0, s*ages, s*ages)
  Z=Is;
  Z[1,1]=0;
  for(i in 1:ages){
    # imom = 1
    E <- Iom[,i] %*% t(Iom[i,]); # al cuadrado?
    momarray[,i] <- kronecker(E,Z) %*% momarray[,i]
  }
  # re-scale
  momarray <- momarray %*% MASS::ginv(diag(colSums(momarray)))

  # considering deaths
  phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,ages*s*2,ages)
  phi[1,1] = 1
  Mtt <- diag(as.numeric(onesom - onesom %*% Ut))
  Utt <- rbind(cbind(Ut,zsom), cbind(Mtt,Isom)) %>% as.matrix()
  ftt <- rbind(cbind(ft,zsom), cbind(zsom,zsom)) %>% as.matrix()
  Gtt <- rbind(cbind(Gt,zsom), cbind(zsom,zsom)) %>% as.matrix()
  sages <- 1:(ages*s)

  # no considering deaths
    # phi = d = gd = ggd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0,ages*s,ages)
    # phi[1,1] = 1
    # Utt = Ut %>% as.matrix()
    # ftt = ft %>% as.matrix()
    # Gtt = Gt %>% as.matrix()

  # focal´s trip
  m[sages,1] = momarray %*% piage;
  for(i in 1:(ages-1)){
    phi[,i+1] = Gtt %*% phi[,i]
    d[,i+1]   = Utt %*% d[,i] + ftt %*% phi[,i]
    gd[,i+1]  = Utt %*% gd[,i] + ftt %*% d[,i]
    ggd[,i+1] = Utt %*% ggd[,i] + ftt %*% gd[,i]
    m[,i+1]   = Utt %*% m[,i]
    ys[,i+1]  = Utt %*% ys[,i] + ftt %*% m[,i]
    nys[,i+1] = Utt %*% nys[,i] + ftt %*% ys[,i]
  }

  gm[,1] = m %*% piage
  for(i in 1:(ages-1)){
    gm[,i+1]  = Utt %*% gm[,i]
  }

  ggm[,1] = gm %*% piage
  for(i in 1:(ages-1)){
    ggm[,i+1]  = Utt %*% ggm[,i]
  }

  os[,1]  = d %*% piage
  nos[,1] = gd[1:ages,] %*% piage
  for(i in 1:(ages-1)){
    os[,i+1]  = Utt %*% os[,i]
    nos[,i+1] = Utt %*% nos[,i] + ftt %*% os[,i]
  }

  oa[,1]  = os %*% piage
  ya[,1]  = ys %*% piage
  coa[,1] = nos %*% piage
  cya[,1] = nys %*% piage
  for(i in 1:(ages-1)){
    oa[,i+1]  = Utt %*% oa[,i]
    ya[,i+1]  = Utt %*% ya[,i]  + ftt %*% gm[,i]
    coa[,i+1] = Utt %*% coa[,i] + ftt %*% oa[,i]
    cya[,i+1] = Utt %*% cya[,i] + ftt %*% ya[,i]
  }

  # get results
  kin_list <- list(d=d,gd=gd,ggd=ggd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                   nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list %>% keep(names(.) %in% output_kin)
  }

  # as data.frame
  kin <- purrr::map2(kin_list, names(kin_list),
              function(x,y){
                out <- as.data.frame(x)
                colnames(out) <- age
                out %>%
                  dplyr::mutate(kin = y,
                         age_kin = rep(sort(rep(age,s)),2),
                         stage_kin = rep(rep(1:s,ages),2),
                         alive = c(rep("living",s*ages),rep("death",s*ages))
                         # age_kin = sort(rep(age,s)),
                         # stage_kin = rep(1:s,ages),
                         # alive = c(rep("yes",s*ages))
                         ) %>%
                  tidyr::pivot_longer(c(-age_kin, -stage_kin, -kin, -alive), names_to = "age_focal", values_to = "count") %>%
                  dplyr::mutate(age_focal = as.integer(age_focal)) %>%
                  tidyr::pivot_wider(names_from = alive, values_from = count)
              }) %>%
    purrr::reduce(rbind)

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    out <- kin
  }

  return(out)
}



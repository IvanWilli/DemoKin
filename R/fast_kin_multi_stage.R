# Fast (C++) implementation of kin_multi_stage(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_multi_stage.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).

# Estimate kin counts by age and stage in a time invariant framework

# @description Implementation of age-stage kin estimates (multi-state) by Caswell (2020). Stages are implied in length of input lists.

# @param U list. age elements with column-stochastic transition matrix with dimension for the state space, conditional on survival.
# @param f matrix. state-specific fertility (age in rows and states in columns). Is accepted also a list with for each age-class.
# @param D matrix. survival probabilities by state (age in rows and states in columns). Is accepted also a list for each state with survival matrices.
# @param H matrix. assigns the offspring of individuals in some stage to the appropriate age class (age in rows and states in columns). Is accepted also a list with a matrix for each state.
# @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See the `vignette` for all kin types.
# @param birth_female numeric. Female portion at birth.
# @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
# @param list_output logical. Results as a list. Default `FALSE`.

# @return A data frame with focal´s age, related ages and type of kin
# (for example `d` is daughter, `oa` is older aunts, etc.), living and death kin counts, and specific stage. If `list_output = TRUE` then this is a list with elements as kin types.

kin_multi_stage_fast <- function(U = NULL, f = NULL, D = NULL, H = NULL,
                            birth_female = 1/2.04,
                            output_kin = NULL,
                            parity = FALSE,
                            list_output = FALSE){

  # mandatory U as a list
  if(!is.list(U)) stop("U must be a list with age length of elements, and stage transitiotn matrix for each one.")

  # stages and age-classes
  s <- ncol(U[[1]])
  ages <- length(U)
  age <- (1:ages)-1

  # build H if it is not already a list
  if(!is.list(H)){
    H <- lapply(1:s, function(Y){
      Ht = matrix(0, nrow=ages, ncol=ages)
      Ht[1,] <- 1
      Ht
    })
  }

  # build D if it is not already a list
  if(!is.list(D)){
    D <- lapply(1:s, function(Y){
      X <- D[,Y]
      Dt = matrix(0, nrow=ages, ncol=ages)
      Dt[row(Dt)-1 == col(Dt)] <- X[-ages]
      Dt[ages, ages] = X[ages]
      Dt
    })
  }

  # build f if it is not already a list
  if(!is.list(f)){
    f <- lapply(1:ages, function(Y){
      X <- f[Y,]
      ft = matrix(0, nrow=s, ncol=s)
      ft[1,] <- X
      ft
    })
  }

  # build block matrices
  bbU <- Matrix::bdiag(U)
  bbF <- Matrix::bdiag(f) * birth_female
  bbD <- Matrix::bdiag(D)
  bbH <- Matrix::bdiag(H)

  # order transitions: first state within age, then age given state
  K <- commutation_matrix(s, ages)
  Ut <- t(K) %*% bbD %*% K %*% bbU
  ft <-  t(K) %*% bbH %*% K %*% bbF

  # focal transition but conditioned to survive
  Gt <-  Ut%*% ginv_of_diag(colSums(as.matrix(Ut)))  # = MASS::ginv(diag(...)), exactly

  # stable distribution of mothers
  At <- Ut + ft
  A_decomp <- eigen(At)
  lambda <- as.double(A_decomp$values[1])
  wt <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
  pi <- wt*At[1,]/sum(wt*At[1,])

  # useful vectors and matrices
  ones   <- t(rep(1,s))
  onesom <- t(rep(1,s*ages))
  onesa  <- t(rep(1,ages))
  Iom    <- diag(1, ages)
  Is     <- diag(1, s)
  Isom   <- diag(1, s*ages)
  zsom   <- matrix(0, s*ages, s*ages)

  # momarray is an array with pi in each column
  piage  <- kronecker(Iom,ones) %*% pi
  momarray <- pi %*% onesa

  # considering deaths (no cumulated): reacreate block struct matrices
  Mtt <- diag(as.numeric(onesom - onesom %*% Ut))
  # dense blocks built directly (Ut, ft and Gt are Matrix objects; binding them
  # with base matrices went through Matrix's rbind/cbind). The three only feed
  # the C++ kernel, which keeps the non-zero entries, and those are the same.
  Utt <- rbind(cbind(as.matrix(Ut),zsom), cbind(Mtt,Isom))
  ftt <- rbind(cbind(as.matrix(ft),zsom), cbind(zsom,zsom))
  Gtt <- rbind(cbind(as.matrix(Gt),zsom), cbind(zsom,zsom))
  sages <- 1:(ages*s)

  # if parity: restriction to no initial mothers with 0 parity
  if(parity){
    Z=Is
    Z[1,1]=0
    if (all(is.finite(momarray))) {
      # kronecker(E, Z) %*% momarray[, i] keeps rows (i-1)*s + 2..s of column i
      # and zeroes the rest. With finite values the BLAS product returns those
      # entries exactly: each is 1 * x plus 0 * x terms added to a sum started
      # at +0, which leaves x unchanged except -0, which becomes +0 (hence
      # "+ 0"). The 606 x 606 kronecker product per age is not built.
      keep <- matrix(FALSE, nrow(momarray), ncol(momarray))
      if (s > 1) for (i in 1:ages) keep[(i - 1) * s + 2:s, i] <- TRUE
      momarray[keep] <- momarray[keep] + 0
      momarray[!keep] <- 0
    } else {
      for(i in 1:ages){
        E <- Iom[,i] %*% t(Iom[i,])
        momarray[,i] <- kronecker(E,Z) %*% momarray[,i]
      }
    }
    # re-scale
    momarray <- momarray %*% ginv_of_diag(colSums(momarray))  # = MASS::ginv(diag(...)), exactly
    # no 0 parity mothers: (momarray %*% piage)[seq(1,600,6)]
    m0 <- as.vector(momarray %*% piage)
  }else{
    m0 <- pi
  }

  # focal´s trip (C++ kernel). Initial conditions use all rows (living and dead).
  n <- ages * s * 2
  piage <- as.vector(piage)
  k <- kin_project_cpp(Utt, ftt, ftt, Gtt,
                       phi0 = c(1, numeric(n - 1)), m0 = m0,
                       pi_all = piage, pi_os = piage,
                       ncol = ages, nlive = n, prev = NULL)
  kin_list <- k[c("Focal", kin_names_14)]
  names(kin_list)[1] <- "focal"

  # only selected kin
  if(!is.null(output_kin)){
    kin_list <- kin_list[names(kin_list) %in% output_kin]
  }

  # results as list?
  if(list_output) {
    out <- kin_list
  }else{
    lay <- list(ids = data.frame(age_kin = rep(age, each = s), stage_kin = rep(1:s, ages)),
                live = sages, dead = list(dead = ages * s + sages))
    out <- build_kin_df_ti(kin_list, age, lay)
  }

  # end
  return(out)
}

# Commutation (vec-permutation) matrix K such that K %*% vec(A) = vec(t(A))
# for an r x c matrix A. Identical to matrixcalc::commutation.matrix(r, c).
commutation_matrix <- function(r, c) {
  K <- matrix(0, r * c, r * c)
  row <- rep(1:r, times = c)
  col <- rep(1:c, each = r)
  K[cbind((row - 1) * c + col, (col - 1) * r + row)] <- 1
  K
}

# Fast (C++) implementation of kin_multi_stage_time_variant_2sex(), reached with `fast = TRUE`.
# The default `fast = FALSE` runs the original R code in R/kin_multi_stage_time_variant_2sex.R.
# Same inputs, same output object; the projection recursions run in
# src/kin_project.cpp and the tables are built directly (R/fast_internal*.R).


# Estimate kin counts by age, stage, and sex, in a time variant framework

# @description Implementation of combined formal demographic models: Caswell II,III,IV.

# @param U_list_females list with matrix entries: period-specific female survival probabilities. Age in rows and states in columns.
# @param U_list_males list with matrix entries: period-specific male survival probabilities. Age in rows and states in columns.
# @param F_list_females list with matrix with elements: period-specific female fertility (age in rows and states in columns).
# @param F_list_males list with matrix entries: period-specific male fertility (age in rows and states in columns).
# @param T_list_females list of lists with matrix entries: each outer list entry is period-specific, and composed of
#                     a list of stochastic matrices which describe age-specific female probabilities of transferring stage
# @param T_list_males list of lists with matrix entries: each outer list entry is period-specific, and composed of
#                     a list of stochastic matrices which describe age-specific male probabilities of transferring stage
# @param birth_female numeric. Female portion at birth. This can be a vector of length equal to the number of years in the input data, or a single value that will be repeated for all years.
# @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
# @param output_kin vector. A vector of particular kin one wishes to obtain results for, e.g., c("m","d","oa"). Default is all kin types.
# @param summary_kin logical. Results as a data frame of accumulated kin by age of Focal if TRUE, and kin by their age*stage distribution by age of Focal if FALSE.
# @param sex_Focal character. Female or Male as the user requests.
# @param initial_stage_Focal Numeric in Natural number set {1,2,...,}. The stage which Focal is born into (e.g., 1 for parity 0)
# @param output_years vector. The times at which we wish to count kin: start year = output_years[1], and end year = output_years[length.]
#
# @return A data frame with focal age, kin age, kin stage, kin sex, year, cohort, and expected number of kin given these restrictions.

#
kin_multi_stage_time_variant_2sex_fast <- function(U_list_females = NULL,
                                              U_list_males = NULL,
                                              F_list_females = NULL,
                                              F_list_males = NULL,
                                              T_list_females = NULL,
                                              T_list_males = NULL,
                                              birth_female = 0.49, ## Sex ratio -- note is 1 - alpha
                                              parity = FALSE,
                                              output_kin = NULL, # enter a vector of specific kin if we only want to analyse these (e.g., c("m","d"))
                                              summary_kin = TRUE, # Set to FALSE if we want only a full age*stage distribution of kin
                                              sex_Focal = "Female", # Female Focal is default
                                              initial_stage_Focal = NULL,
                                              output_years){

  # SD_fast() keeps its last eigen decomposition for the boundary year; drop it on exit
  on.exit({environment(SD_fast)$last_M <- NULL; environment(SD_fast)$last_eig <- NULL}, add = TRUE)
  no_years <- length(U_list_females)
  na <- nrow(U_list_females[[1]])
  ns <- ncol(U_list_females[[1]])

  # Ensure inputs are lists of matrices and that the timescale same length
  if(!is.list(U_list_females) | !is.list(U_list_males)){stop("U's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(F_list_females) | !is.list(F_list_males)){stop("F's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(T_list_females) | !is.list(T_list_males)){stop("T's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}

  # birth_female needs to be dynamic. Complete in case length is lower than data
  if(length(birth_female) < no_years){
    last_birth_female <- utils::tail(birth_female, n=1)
    n_birth_female <- length(birth_female)
    birth_female <- c(birth_female, rep(last_birth_female, no_years - n_birth_female))
    message("Length of birth_female lower than years of risk. Completed with last value...")
  }

  # periods (indices into the series of kin arrays) requested as output; only those are kept
  start_year <- output_years[1]
  needed <- as.numeric(output_years) - as.numeric(start_year) + 1
  kept <- list()

  message("Running kinship model...")
  prev <- NULL
  for(year in 1:no_years){
    ## U_tilde transfers across stage and advances age; F_tilde makes newborns from stage/age and puts them in age 0
    U_tilde_females <- build_U_tilde(U_list_females[[year]], T_list_females[[year]], na, ns)
    U_tilde_males   <- build_U_tilde(U_list_males[[year]],   T_list_males[[year]],   na, ns)
    F_tilde_females <- build_F_tilde(F_list_females[[year]], na, ns)
    F_tilde_males   <- build_F_tilde(F_list_males[[year]],   na, ns)

    ## if year == 1 we are at the boundary condition t=0 apply time-invariant kinship projections
    if(year == 1){
      prev <- all_kin_dy_fast(U_tilde_females, U_tilde_males, F_tilde_females, F_tilde_males,
                         1-birth_female[1], na, ns, parity, sex_Focal, initial_stage_Focal)
      if(1 %in% needed) kept[["1"]] <- prev
    }
    ## Output of the time-variant model
    kin_out <- all_kin_dy_TV_fast(U_tilde_females, U_tilde_males, F_tilde_females, F_tilde_males,
                             1-birth_female[year], na, ns, parity, sex_Focal, initial_stage_Focal, prev)
    if((year + 1) %in% needed) kept[[as.character(year + 1)]] <- kin_out
    prev <- kin_out
  }

  message("Preparing output...")

  relative_names <- c("Focal", "d", "gd", "ggd", "m", "gm", "ggm", "ys", "os", "ya", "oa", "nys", "nos", "cya", "coa")
  kin_full <- create_full_dists_df_fast(kept, relative_names, output_years, start_year, na, ns, output_kin)
  if(summary_kin){
    kin_summary <- create_cumsum_df_fast(kept, relative_names, output_years, start_year, na, ns, output_kin)
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)}
  else{
    kin_out <- kin_full
  }
  return(kin_out)
}

# U_tilde = t(K) %*% bdiag_s(U_s) %*% K %*% bdiag_a(T_a) written out: every entry is the
# single product U[a, s'] * T_a[s', s] (age advances a -> a+1, last age is absorbing),
# so building it directly gives exactly the same matrix as the permutation products.
build_U_tilde <- function(Umat, Tlist, na, ns){
  Umat <- as.matrix(Umat)
  a_from <- rep(1:na, each = ns * ns)
  a_to   <- pmin(a_from + 1, na)
  s_row  <- rep(rep(1:ns, times = ns), na)   # s' (stage of destination), fastest
  s_col  <- rep(rep(1:ns, each = ns), na)    # s  (stage of origin)
  Tv <- unlist(lapply(1:na, function(a) as.vector(as.matrix(Tlist[[a]]))), use.names = FALSE)
  x <- Umat[cbind(a_from, s_row)] * Tv
  keep <- x != 0
  Matrix::sparseMatrix(i = ((a_to - 1) * ns + s_row)[keep], j = ((a_from - 1) * ns + s_col)[keep],
                       x = x[keep], dims = c(na * ns, na * ns))
}

# F_tilde = t(K) %*% bdiag_s(H) %*% K %*% bdiag_a(F_a): newborns of every age*stage enter
# age 1, stage 1, at rate 1 * F[a, s].
build_F_tilde <- function(Fmat, na, ns){
  x <- as.vector(t(as.matrix(Fmat)))
  keep <- x != 0
  Matrix::sparseMatrix(i = rep(1L, sum(keep)), j = which(keep), x = x[keep], dims = c(na * ns, na * ns))
}

# triplet form of a sparse Matrix for the C++ kernel
csr_input <- function(x){
  x <- methods::as(methods::as(x, "CsparseMatrix"), "TsparseMatrix")
  list(i = as.integer(x@i + 1L), j = as.integer(x@j + 1L), x = as.numeric(x@x), nrow = nrow(x), ncol = ncol(x))
}

# Two-sex block projection matrices (sparse), as in the original all_kin_dy_fast()
two_sex_blocks <- function(Uf, Um, Ff, Fm, alpha){
  n <- nrow(Uf)
  Z <- Matrix::Matrix(0, n, n, sparse = TRUE)
  Uproj <- Matrix::bdiag(Uf, Um)
  Fproj <- rbind(cbind((1-alpha)*Ff, (1-alpha)*Fm), cbind(alpha*Ff, alpha*Fm))
  Fprojstar <- rbind(cbind((1-alpha)*Ff, Z), cbind(alpha*Ff, Z))
  list(Uproj = Uproj, Fproj = Fproj, Fprojstar = Fprojstar)
}

# Title time invariant two-sex multi-state kin projections
#
# @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
# @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
# @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
# @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
# @param alpha scalar. birth ratio (male:female)
# @param na scalar. number of ages.
# @param ns scalar. number of stages.
# @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
# @param sex_Focal logical. Female or Male
# @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
#
# @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
# nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
# yielding the age*stage distribution of kin for each age of Focal

all_kin_dy_fast <- function(Uf, Um, Ff, Fm, alpha, na, ns, Parity, sex_Focal, Initial_stage_Focal){

  n <- nrow(Uf) ## number of ages * stages for each sex
  B <- two_sex_blocks(Uf, Um, Ff, Fm, alpha)

  ## The stable population structure is an age*stage*sex vector
  population_age_stage_structure <- SD_fast(B$Uproj + B$Fprojstar)

  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1
    population_age_stage_of_parenting <- pi_mix_parity_fast(Uf, Um, Ff, Fm, alpha, na, ns)
  }
  else{
    population_age_stage_of_parenting <- pi_mix_fast(Uf, Um, Ff, Fm, alpha, na, ns)
  }
  parents_joint_age_stage <- as.vector(as.matrix(population_age_stage_of_parenting[[1]]))
  mothers_age_dist <- as.vector(as.matrix(population_age_stage_of_parenting[[4]]))
  fathers_age_dist <- as.vector(as.matrix(population_age_stage_of_parenting[[5]]))

  ### Focal requires its own dynamic: G_tilde tracks Focal's age*stage advancement
  G_tilde <- matrix(0, 2*n, 2*n)
  G_tilde[1:n, 1:n] <- as.matrix(get_G_fast(Uf, na, ns))
  G_tilde[(n+1):(2*n), (n+1):(2*n)] <- as.matrix(get_G_fast(Um, na, ns))
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }

  out <- kin_project_ms2sex_cpp(csr_input(B$Uproj), csr_input(B$Fproj), csr_input(B$Fprojstar), G_tilde,
                                IC_Focal, parents_joint_age_stage, mothers_age_dist, fathers_age_dist, na, NULL)
  out[["ps"]] <- population_age_stage_structure
  out
}

# Title time-variant two-sex multi-state kin projections
#
# @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
# @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
# @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
# @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
# @param alpha scalar. birth ratio (male:female)
# @param na scalar. number of ages.
# @param ns scalar. number of stages.
# @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
# @param sex_Focal logical. Female or Male
# @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
# @param prev list. last years kinship output (kin matrices and population structure `ps`).
#
# @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
# nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
# yielding the age*stage distribution of kin for each age of Focal
#
all_kin_dy_TV_fast <- function(Uf, Um, Ff, Fm, alpha, na, ns, Parity, sex_Focal, Initial_stage_Focal, prev){

  n <- nrow(Uf)
  B <- two_sex_blocks(Uf, Um, Ff, Fm, alpha)

  population_age_stage_structure <- prev[["ps"]]
  population_age_stage_structure <- population_age_stage_structure/sum(population_age_stage_structure)
  population_age_stage_structure_next <- (B$Uproj + B$Fprojstar)%*%population_age_stage_structure

  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1
    population_age_stage_of_parenting <- pi_mix_TV_parity_fast(Ff, Fm, alpha, na, ns, population_age_stage_structure)
  }
  else{
    population_age_stage_of_parenting <- pi_mix_TV_fast(Ff, Fm, alpha, na, ns, population_age_stage_structure)
  }
  parents_joint_age_stage <- as.vector(as.matrix(population_age_stage_of_parenting[[1]]))
  mothers_age_dist <- as.vector(as.matrix(population_age_stage_of_parenting[[4]]))
  fathers_age_dist <- as.vector(as.matrix(population_age_stage_of_parenting[[5]]))

  ### Focal requires its own dynamic: G_tilde tracks Focal's age*stage advancement
  G_tilde <- matrix(0, 2*n, 2*n)
  G_tilde[1:n, 1:n] <- as.matrix(get_G_fast(Uf, na, ns))
  G_tilde[(n+1):(2*n), (n+1):(2*n)] <- as.matrix(get_G_fast(Um, na, ns))
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }

  out <- kin_project_ms2sex_cpp(csr_input(B$Uproj), csr_input(B$Fproj), csr_input(B$Fprojstar), G_tilde,
                                IC_Focal, parents_joint_age_stage, mothers_age_dist, fathers_age_dist, na, prev)
  out[["ps"]] <- population_age_stage_structure_next
  out
}

################## Create data frame output

# Title Accumulated kin by each age of Focal, for each time period, and cohort of birth
#
# @param kin_matrix_lists list by period index (character names) of kin matrix lists.
# @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
# @param years vector. The timescale on which we implement the kinship model.
# @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
# @param na numeric. Number of ages.
# @param ns numeric. Number of stages.
# @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#
# @return A data frame which gives for each age of Focal at each year in the timescale, Focal's experienced number kin demarcated by stages (summed over all ages)
#
create_cumsum_df_fast <- function(kin_matrix_lists, kin_names, years, start_year, na, ns, specific_kin = NULL){
  # only the kin asked for are built (the same rows, in the same order, as
  # building all of them and filtering); if none matches, the original path
  # below keeps the empty table's columns
  if(!is.null(specific_kin) && any(unlist(kin_names) %in% specific_kin)){
    kin_names <- kin_names[unlist(kin_names) %in% specific_kin]
    specific_kin <- NULL
  }
  n <- na * ns
  nk <- length(kin_names)
  nblock <- na * ns * 2
  stage_f <- factor(rep(seq(1, ns), na), levels = seq(1, ns))
  # sums over age within (Focal age, stage, sex): rows ordered by Focal age, then stage, then f/m
  living <- unlist(lapply(years, function(j) {
    ii <- as.numeric(j) - start_year + 1
    unlist(lapply(kin_names, function(kin_member) {
      X <- kin_matrix_lists[[as.character(ii)]][[kin_member]]
      sf <- stage_sums_cpp(X[1:n, , drop = FALSE], ns)
      sm <- stage_sums_cpp(X[(n+1):(2*n), , drop = FALSE], ns)
      as.vector(rbind(as.vector(sf), as.vector(sm)))
    }), use.names = FALSE)
  }), use.names = FALSE)
  ny <- length(years)
  df <- dplyr::tibble(age_focal = rep(rep(as.numeric(seq(0, na - 1)), each = 2 * ns), nk * ny),
                      stage_kin = stage_f[rep(rep(rep(seq_len(ns), each = 2), na), nk * ny)],
                      living = living,
                      sex_kin = rep(c("f", "m"), na * ns * nk * ny),
                      year = rep(years, each = nblock * nk),
                      kin = rep(rep(unlist(kin_names), each = nblock), ny))
  df$cohort <- as.numeric(df$year) - as.numeric(df$age_focal)   # = dplyr::mutate()
  if(!is.null(specific_kin)){
    df <- df %>% dplyr::filter(kin %in% specific_kin)
  }
  return(df)
}

# Title joint age*stage distributions of kin by each age of Focal, for each time period, and cohort of birth
#
# @param kin_matrix_lists list by period index (character names) of kin matrix lists.
# @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
# @param years vector. The timescale on which we implement the kinship model.
# @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
# @param na numeric. Number of ages.
# @param ns numeric. Number of stages.
# @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#
# @return A data frame which gives for each age of Focal at each year in the timescale, the full age*stage dist of kin
#
create_full_dists_df_fast <- function(kin_matrix_lists, kin_names, years, start_year, na, ns, specific_kin = NULL){
  # only the kin asked for are built (the same rows, in the same order, as
  # building all of them and filtering); if none matches, the original path
  # below keeps the empty table's columns
  if(!is.null(specific_kin) && any(unlist(kin_names) %in% specific_kin)){
    kin_names <- kin_names[unlist(kin_names) %in% specific_kin]
    specific_kin <- NULL
  }
  n <- na * ns
  nk <- length(kin_names)
  ny <- length(years)
  nblock <- 2 * n * na
  age_kin <- rep(rep(seq(0, (na-1)), each = ns), 2)
  stage_f <- factor(rep(seq(1, ns), na), levels = seq(1, ns))
  stage_kin <- stage_f[rep(seq_len(n), 2)]
  sex_kin <- rep(c("f", "m"), each = n)
  living <- unlist(lapply(years, function(j) {
    ii <- as.numeric(j) - start_year + 1
    unlist(lapply(kin_names, function(kin_member) as.vector(kin_matrix_lists[[as.character(ii)]][[kin_member]])),
           use.names = FALSE)
  }), use.names = FALSE)
  df <- data.frame(age_focal = rep(rep(as.numeric(seq(0, na - 1)), each = 2 * n), nk * ny),
                   age_kin = rep(age_kin, na * nk * ny),
                   stage_kin = stage_kin[rep(seq_len(2 * n), na * nk * ny)],
                   living = living,
                   sex_kin = rep(sex_kin, na * nk * ny),
                   year = rep(years, each = nblock * nk),
                   kin = rep(rep(unlist(kin_names), each = nblock), ny),
                   stringsAsFactors = FALSE)
  df$cohort <- as.numeric(df$year) - as.numeric(df$age_focal)   # = dplyr::mutate()
  if(!is.null(specific_kin)){
    df <- df %>% dplyr::filter(kin %in% specific_kin)
  }
  return(df)
}

# Zero every entry of column i except age i, stages 2..ns (parity > 0)
#
# For a dense dgeMatrix (what pi_f %*% t(1) gives) the kept entries are
# copied and the rest set to +0 directly in the value slot, which is what the
# column loop below did through `[<-` (101 S4 sub-assignments, a large share
# of the parity model's run time). Anything else takes the loop.
keep_parity_block <- function(M, na, ns){
  if (identical(class(M), structure("dgeMatrix", package = "Matrix")) &&
      nrow(M) == na * ns && ncol(M) == na && length(M@factors) == 0L) {
    keep <- matrix(FALSE, na * ns, na)
    if (ns > 1) for (i in 1:na) keep[(i - 1) * ns + 2:ns, i] <- TRUE
    x <- M@x
    x[!as.vector(keep)] <- 0
    M@x <- x
    return(M)
  }
  for(i in 1:na){
    v <- M[,i]
    v2 <- numeric(length(v))
    if(ns > 1){
      keep <- (i-1)*ns + 2:ns
      v2[keep] <- v[keep]
    }
    M[,i] <- v2
  }
  M
}

## Construct a matrix which transfers Focal across stages, while ensuring Focal survives with probability 1
get_G_fast <- function(U, na, ns){
  sig <- Matrix::t(rep(1,na*ns)) %*% U
  # MASS::ginv(Matrix::diag(sig[1,])), exactly, without the SVD
  G <- U %*% ginv_of_diag(sig[1,], Matrix::diag)
  return(G)
}

# Mixing distributions for the time-invariant multi-state 2-sex model: Non-parity case
#
# @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
# @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
# @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
# @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
# @param alpha scalar. Birth ratio male:female
# @param na scalar. Number of age-classes
# @param ns scalar. Number of stages
#
# @return list (of vectors). list[[1]] = full age*stage*sex distribution
#                            list[[2]] = female age*stage distribution normalised
#                            list[[3]] = male age*stage distribution normalised
#                            list[[4]] = female marginal age distribution normalised
#                            list[[5]] = male marginal age distribution normalised
#
pi_mix_fast <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD_fast(A)
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(rbind(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

# Mixing distributions for the time-variant multi-state 2-sex model: Non-parity case
#
# @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
# @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
# @param alpha scalar. Birth ratio male:female
# @param na scalar. Number of age-classes
# @param ns scalar. Number of stages
# @param previous_age_stage_dist vector. Last years population structure (age*stage*sex full distribution)
#
# @return list (of vectors). list[[1]] = full age*stage*sex distribution
#                            list[[2]] = female age*stage distribution normalised
#                            list[[3]] = male age*stage distribution normalised
#                            list[[4]] = female marginal age distribution normalised
#                            list[[5]] = male marginal age distribution normalised
#
pi_mix_TV_fast <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(rbind(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

# Mixing distributions for the time-invariant multi-state 2-sex model: Parity-specific case
#
# @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
# @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
# @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
# @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
# @param alpha scalar. Birth ratio male:female
# @param na scalar. Number of age-classes
# @param ns scalar. Number of stages
#
# @return list (of vectors). list[[1]] = full age*stage*sex distribution
#                            list[[2]] = female age*stage distribution normalised
#                            list[[3]] = male age*stage distribution normalised
#                            list[[4]] = female marginal age distribution normalised
#                            list[[5]] = male marginal age distribution normalised
#
pi_mix_parity_fast <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD_fast(A)
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  # keep, in column i, only age i and stages 2..ns (exactly what the kronecker
  # projection did: entries are copied, everything else is zero)
  m_mat <- keep_parity_block(m_mat, na, ns)
  d_mat <- keep_parity_block(d_mat, na, ns)
  out_mum <- m_mat %*% ginv_of_diag(Matrix::colSums(m_mat), Matrix::diag)
  out_dad <- d_mat %*% ginv_of_diag(Matrix::colSums(d_mat), Matrix::diag)
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(rbind(pi_f,pi_m), pi_f, pi_m, pi_F, pi_M))
}

# Mixing distributions for the time-variant multi-state 2-sex model: Parity-specific case
#
# @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
# @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
# @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
# @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
# @param alpha scalar. Birth ratio male:female
# @param na scalar. Number of age-classes
# @param ns scalar. Number of stages
# @param previous_age_stage_dist vector. Last years population structure (age*stage*sex full distribution)
#
# @return list (of vectors). list[[1]] = full age*stage*sex distribution
#                            list[[2]] = female age*stage distribution normalised
#                            list[[3]] = male age*stage distribution normalised
#                            list[[4]] = female marginal age distribution normalised
#                            list[[5]] = male marginal age distribution normalised
#
pi_mix_TV_parity_fast <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F  <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  # keep, in column i, only age i and stages 2..ns (exactly what the kronecker
  # projection did: entries are copied, everything else is zero)
  m_mat <- keep_parity_block(m_mat, na, ns)
  d_mat <- keep_parity_block(d_mat, na, ns)
  out_mum <- m_mat %*% ginv_of_diag(Matrix::colSums(m_mat), Matrix::diag)
  out_dad <- d_mat %*% ginv_of_diag(Matrix::colSums(d_mat), Matrix::diag)
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(rbind(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

######################################################### Some useful utility functions required

# Find the column-eigenvector corresponding to the spectral radius (Stable population structure in Demographics)
#
# eigen() works on unname(as.matrix(PM)), so the result depends on that dense
# matrix only. The last one is kept with its result: the boundary year of
# kin_multi_stage_time_variant_2sex() asks for the same matrix twice (the
# population structure in all_kin_dy_fast() and the parents' distribution in
# pi_mix_fast()/pi_mix_parity_fast(), built differently but identical), and the second
# eigen() was a large share of the run. A different matrix is decomposed as
# before.
SD_fast <- local({
  last_M <- NULL
  last_eig <- NULL
  function(PM) {
    M <- unname(as.matrix(PM))
    if (!is.null(last_M) && identical(M, last_M)) {
      spectral_stuff <- last_eig
    } else {
      spectral_stuff <- eigen(M)
      last_M <<- M
      last_eig <<- spectral_stuff
    }
    spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
    # normalise...
    vec_lambda <- spectral_stuff/sum(spectral_stuff)
    return(vec_lambda)
  }
})
# block_diag_function() and E_matrix() are shared with the original code (R/kin_multi_stage_time_variant_2sex.R).

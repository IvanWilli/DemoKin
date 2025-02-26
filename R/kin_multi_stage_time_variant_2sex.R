

#' Estimate kin counts by age, stage, and sex, in a time variant framework

#' @description Implementation of combined formal demographic models: Caswell II,III,IV.

#' @param U_list_females list with matrix entries: period-specific female survival probabilities. Age in rows and states in columns.
#' @param U_list_males list with matrix entries: period-specific male survival probabilities. Age in rows and states in columns.
#' @param F_list_females list with matrix with elements: period-specific female fertility (age in rows and states in columns).
#' @param F_list_males list with matrix entries: period-specific male fertility (age in rows and states in columns).
#' @param T_list_females list of lists with matrix entries: each outer list entry is period-specific, and composed of
#'                     a list of stochastic matrices which describe age-specific female probabilities of transferring stage
#' @param T_list_males list of lists with matrix entries: each outer list entry is period-specific, and composed of
#'                     a list of stochastic matrices which describe age-specific male probabilities of transferring stage
#' @param H_list list with matrix entries: redistribution of newborns across each stage to a specific age-class
#' @param birth_female numeric. birth ratio of females to males in population
#' @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
#' @param output_kin vector. A vector of particular kin one wishes to obtain results for, e.g., c("m","d","oa"). Default is all kin types.
#' @param summary_kin logical. Results as a data frame of accumulated kin by age of Focal if TRUE, and kin by their age*stage distribution by age of Focal if FALSE.
#' @param sex_Focal character. Female or Male as the user requests.
#' @param initial_stage_Focal Numeric in Natural number set {1,2,...,}. The stage which Focal is born into (e.g., 1 for parity 0)
#' @param output_years vector. The times at which we wish to count kin: start year = output_years[1], and end year = output_years[length.]
#' @param model_years. The full timescale on which we run the matrix model. From these periods we extract the ``output_years''.
#'                     Note that if we use abridged life-tables: e.g., 1960,1965,1970 to run the model, by default age_increment = 5
#' @param age_year_consistent logical. Null sets age-bridge to be equal to year
#' @param age_increment. numeric. If age_year_consisent FALSE set own age-gap
#' @return A data frame with focal age, kin age, kin stage, kin sex, year, cohort, and expected number of kin given these restrictions.

#' @export
#'
kin_multi_stage_time_variant_2sex <- function(U_list_females = NULL,
                                              U_list_males = NULL,
                                              F_list_females = NULL,
                                              F_list_males = NULL,
                                              T_list_females = NULL,
                                              T_list_males = NULL,
                                              H_list = NULL,
                                              birth_female = 0.49, ## Sex ratio -- note is 1 - alpha
                                              parity = FALSE,
                                              output_kin = NULL, # enter a vector of specific kin if we only want to analyse these (e.g., c("m","d"))
                                              summary_kin = TRUE, # Set to FALSE if we want a full age*stage distribution of kin
                                              sex_Focal = "Female", # Female Focal is default
                                              initial_stage_Focal = NULL,
                                              output_years,
                                              model_years,
                                              age_year_consitent = TRUE,
                                              age_increment = NULL){

  no_years <- (-1+length(U_list_females))
  na <- nrow(U_list_females[[1]])
  ns <- ncol(U_list_females[[1]])
  if(age_year_consitent){age_increment <- as.numeric(model_years[2]-model_years[1])}
  # Ensure inputs are lists of matrices and that the timescale same length

  if(length(U_list_females) != (length(model_years))){stop("Proposed timescale longer than demographic timescale")} ## this is due to my struggles with counting! ( e.g., seq(10, 20, 1) != list(1 : 10) )
  if(output_years[length(output_years)] > model_years[length(model_years)]){stop("Output years longer than model run")}

  if(!is.list(U_list_females) | !is.list(U_list_males)){stop("U's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(F_list_females) | !is.list(F_list_males)){stop("F's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(T_list_females) | !is.list(T_list_males)){stop("T's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}

  ### Define empty lists for the accumulated kin of Focals's life-course -- each list entry will reflect a time-period
  changing_pop_struct <- list()
  Focal_array <- list()
  mom_array <- list()
  gran_array <- list()
  great_gran_array <- list()
  daughter_array <- list()
  younger_sis_array <- list()
  grand_daughter_array <-list()
  great_grand_daughter_array <- list()
  older_sister_array <- list()
  younger_aunt_array <- list()
  older_aunt_array <- list()
  younger_niece_array <- list()
  older_niece_array <- list()
  younger_cousin_array <- list()
  older_cousin_array <- list()

  ### At each time-period we: 1) -- construct the time-variant projection matrices:
  ###                               U_tilde : transfers across stage and advances age
  ###                               F_tilde : makes newborns from stage/age; puts them to stage/age
  ###                         2) -- project Focal and kin using above projection matrices

  pb <- progress::progress_bar$new(
    format = "Timescale [:bar] :percent",
    total = no_years + 1, clear = FALSE, width = 60)
  tictoc::tic()
  for(year in 1:no_years){
    pb$tick()

    T_data_f <- T_list_females[[year]] ## For each year we have na number of Transfer matrices
    T_data_m <- T_list_males[[year]]   ## which give probabilities of age-dep movement from stage to stage
    T_f_list <- list()
    T_m_list <- list()
    F_f_list <- list()
    F_m_list <- list()
    U_f_list <- list()
    U_m_list <- list()
    H_list2 <- list()

    for(stage in 1:ns){
      Uf <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      Matrix::diag(Uf[-1,-ncol(Uf)]) <- U_list_females[[year]][1:(na-1),stage]
      Uf[na,na] <- U_list_females[[year]][na,stage]
      Um <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      Matrix::diag(Um[-1,-ncol(Um)]) <- U_list_males[[year]][1:(na-1),stage]
      Um[na,na] <- U_list_males[[year]][na,stage]
      U_f_list[[(1+length(U_f_list))]] <- Uf
      U_m_list[[(1+length(U_m_list))]] <- Um
      H_mat <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      H_mat[1,] <- 1
      H_list2[[(1+length(H_list2))]] <- H_mat
    }
    for(age in 1:na){
      T_f <- T_data_f[[age]]
      T_m <- T_data_m[[age]]
      T_f_list[[(1+length(T_f_list))]] <- T_f
      T_m_list[[(1+length(T_m_list))]] <- T_m
      F_f <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
      F_m <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
      F_f[1,] <- F_list_females[[year]][age,]
      F_m[1,] <- F_list_males[[year]][age,]
      F_f_list[[(1+length(F_f_list))]] <- F_f
      F_m_list[[(1+length(F_m_list))]] <- F_m
    }
    ## create the appropriate block-diagonal matrices
    U_f_BDD <- block_diag_function(U_f_list) ## direct sum of female survivorship, independent over stage (ns diagonal blocks)
    U_m_BDD <- block_diag_function(U_m_list) ## direct sum of male survivorship, independent over stage (ns diagonal blocks)
    H_BDD <- block_diag_function(H_list2) ## direct sum of which age newborns enter, independent over stage (ns diagonal blocks)
    T_f_BDD <- block_diag_function(T_f_list) ## direct sum of female stage transitions, independent over age (na diagonal blocks)
    T_m_BDD <- block_diag_function(T_m_list) ## direct sum of male stage transitions, independent over age (na diagonal blocks)
    F_f_BDD <- block_diag_function(F_f_list) ## direct sum of female stage->stage reproductions, independent over age (na blocks)
    F_m_BDD <- block_diag_function(F_m_list) ## direct sum of male stage->stage reproductions, independent over age (na blocks)

    ## create the appropriate projection matrices
    U_tilde_females <- Matrix::t(K_perm_mat(ns, na)) %*%
      U_f_BDD %*%
      K_perm_mat(ns, na) %*%
      T_f_BDD

    ## create sex-specific age*stage projections
    U_tilde_males <- Matrix::t(K_perm_mat(ns, na)) %*%
      U_m_BDD %*%
      K_perm_mat(ns, na) %*%
      T_m_BDD

    F_tilde_females <- Matrix::t(K_perm_mat(ns, na)) %*%
      H_BDD %*%
      K_perm_mat(ns, na) %*%
      F_f_BDD

    F_tilde_males <- Matrix::t(K_perm_mat(ns, na)) %*%
      H_BDD %*%
      K_perm_mat(ns, na) %*%
      F_m_BDD

    ## if year == 1 we are at the boundary condition t=0 apply time-invariant kinship projections
    if(year == 1){
      ## Output of the static model
      kin_out_1 <- all_kin_dy(U_tilde_females,
                              U_tilde_males ,
                              F_tilde_females,
                              F_tilde_males,
                              1-birth_female,
                              na,
                              ns,
                              parity,
                              sex_Focal,
                              initial_stage_Focal)
      ### Relative lists' first entries
      Focal_array[[(1+length(Focal_array))]] <- kin_out_1[["Focal"]]
      daughter_array[[(1+length(daughter_array))]] <- kin_out_1[["d"]]
      grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out_1[["gd"]]
      great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out_1[["ggd"]]
      mom_array[[(1+length(mom_array))]] <- kin_out_1[["m"]]
      gran_array[[(1+length(gran_array))]] <- kin_out_1[["gm"]]
      great_gran_array[[(1+length(great_gran_array))]] <- kin_out_1[["ggm"]]
      younger_sis_array[[( 1+length(younger_sis_array))]] <- kin_out_1[["ys"]]
      older_sister_array[[(1+length(older_sister_array))]] <- kin_out_1[["os"]]
      younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out_1[["ya"]]
      older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out_1[["oa"]]
      younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out_1[["nys"]]
      older_niece_array[[(1+length(older_niece_array))]] <- kin_out_1[["nos"]]
      younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out_1[["cya"]]
      older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out_1[["coa"]]
      changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out_1[["ps"]]

    }
    updating_Focal <- Focal_array[[year]]
    updating_daughter <- daughter_array[[year]]
    updating_grand_daughter <- grand_daughter_array[[year]]
    updating_great_grand_daughter <- great_grand_daughter_array[[year]]
    updating_mom <- mom_array[[year]]
    updating_gran <- gran_array[[year]]
    updating_great_gran <- great_gran_array[[year]]
    updating_younger_sis <- younger_sis_array[[year]]
    updating_older_sis <- older_sister_array[[year]]
    updating_youner_aunt <- younger_aunt_array[[year]]
    updating_older_aunt <- older_aunt_array[[year]]
    updating_younger_niece <- younger_niece_array[[year]]
    updating_older_niece <- older_niece_array[[year]]
    updating_younger_cousin <- younger_cousin_array[[year]]
    updating_older_cousin <- older_cousin_array[[year]]
    updating_pop_struct <- changing_pop_struct[[year]]

    ## Output of the time-variant model
    kin_out <- all_kin_dy_TV(U_tilde_females,
                             U_tilde_males,
                             F_tilde_females,
                             F_tilde_males,
                             1-birth_female,
                             na,
                             ns,
                             parity,
                             sex_Focal,
                             initial_stage_Focal,
                             updating_Focal,
                             updating_daughter,
                             updating_grand_daughter,
                             updating_great_grand_daughter,
                             updating_mom,
                             updating_gran,
                             updating_great_gran,
                             updating_older_sis,
                             updating_younger_sis,
                             updating_older_niece,
                             updating_younger_niece,
                             updating_older_aunt,
                             updating_youner_aunt,
                             updating_older_cousin,
                             updating_younger_cousin,
                             updating_pop_struct)
    ## Relative lists entries correspond to timescale periods (each entry an kin age*stage*2 by Focal age matrix)
    Focal_array[[(1+length(Focal_array))]] <- kin_out[["Focal"]]
    daughter_array[[(1+length(daughter_array))]] <- kin_out[["d"]]
    grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out[["gd"]]
    great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out[["ggd"]]
    mom_array[[(1+length(mom_array))]] <- kin_out[["m"]]
    gran_array[[(1+length(gran_array))]] <- kin_out[["gm"]]
    great_gran_array[[(1+length(great_gran_array))]] <- kin_out[["ggm"]]
    younger_sis_array[[(1+length(younger_sis_array))]] <- kin_out[["ys"]]
    older_sister_array[[(1+length(older_sister_array))]] <- kin_out[["os"]]
    younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out[["ya"]]
    older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out[["oa"]]
    younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out[["nys"]]
    older_niece_array[[(1+length(older_niece_array))]] <- kin_out[["nos"]]
    younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out[["cya"]]
    older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out[["coa"]]
    changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out[["ps"]]
    pb$tick()
  }
  tictoc::toc()
  ## create a list of output kin -- each element a time-period specific list of matrices
  ## label the kin names to match DemoKin:
  relative_data <- list("Focal" = Focal_array,
                        "d" = daughter_array,
                        "gd" = grand_daughter_array,
                        "ggd" = great_grand_daughter_array,
                        "m" = mom_array,
                        "gm" = gran_array,
                        "ggm" = great_gran_array,
                        "ys" = younger_sis_array,
                        "os" = older_sister_array,
                        "ya" = younger_aunt_array,
                        "oa" = older_aunt_array,
                        "nys" = younger_niece_array,
                        "nos" = older_niece_array,
                        "cya" = younger_cousin_array,
                        "coa" = older_cousin_array)

  relative_names <- names(relative_data)
  ## create a nice data frame output
  kin_full <- create_full_dists_df(relative_data,
                                   relative_names,
                                   output_years,
                                   model_years,
                                   na,
                                   ns,
                                   output_kin,
                                   age_increment)
  if(summary_kin){
    kin_summary <- create_cumsum_df(relative_data,
                                    relative_names,
                                    output_years,
                                    model_years,
                                    na,
                                    ns,
                                    output_kin,
                                    age_increment)
    kin_out <- list(kin_full = kin_full, kin_summary = kin_summary)}
  else{
    kin_out <- kin_full
  }
  return(kin_out)
}


#' Title time invariant two-sex multi-state kin projections
#'
#' @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
#' @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
#' @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
#' @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
#' @param alpha scalar. birth ratio (male:female)
#' @param na scalar. number of ages.
#' @param ns scalar. number of stages.
#' @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
#' @param sex_Focal logical. Female or Male
#' @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
#'
#' @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
#' nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
#' yielding the age*stage distribution of kin for each age of Focal

all_kin_dy <- function(Uf,
                       Um,
                       Ff,
                       Fm,
                       alpha, ## alpha = sex ratio male:female (i.e., 1 - birth_female)
                       na, ## na = number of ages
                       ns, ## ns = number of stages
                       Parity,
                       sex_Focal, ## binary "F" or "M"
                       Initial_stage_Focal){

  n <- nrow(Uf) ## number of ages * stages for each sex

  ## Projection matrices:

  ## Uproj is a block diagonal matrix of block-structured Age*Stage matrices; independently over sex transfers individuals across stage and up age
  Uproj <- Matrix::Matrix(block_diag_function(list(Uf, Um)), sparse = TRUE)
  ## Fproj is a Sex-block-structured matrix of block-structured Age*Stage matrices where males and females BOTH reproduce (by stage)
  Fproj <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  Fproj[1:n, 1:n] <- (1-alpha)*Ff ## Ff is Age*Stage block structured giving rate at which females in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm ## Fm is Age*Stage block structured giving rate at which males in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm

  ## Fprojstar is a Sex-block-structured matrix of block-structured Age*Stage matrices where ONLY females reproduce
  Fprojstar <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE) ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff

  ## The stable population structure is an age*stage*sex vector:
  ##                                                            1:n gives the female age*stage structure
  ##                                                            (1+n):2n gives the male age*stage structure
  population_age_stage_structure <- SD(Uproj + Fprojstar)

  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1

    population_age_stage_of_parenting <- pi_mix_parity(Uf, Um, Ff, Fm, alpha, na, ns)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]

    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]

  }
  else{
    population_age_stage_of_parenting <- pi_mix(Uf, Um, Ff, Fm, alpha, na, ns)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]

    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]

  }

  ####################################### The dynamics of Kinship, starting with Focal who is no longer a unit vector

  ### Focal requires its own dynamic: G_tilde constructed below tracks Focal's age*stage advancement over the time-scale
  f_t <- get_G(Uf, na, ns) ## get_G function in "Functions_required.R"
  m_t <- get_G(Um, na, ns)
  G_tilde <- block_diag_function(list(f_t,m_t))
  X_Focal  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }

  ### empty kin matrices for all of Focal's kin
  X_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_sibs  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_sibs  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_niece_nephew  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_niece_nephew  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_aunt_uncle  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_aunts_uncles <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_cousins  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_cousins  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)


  ### Initial distributions for kin with non-zero deterministic initial conditions:
  # Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage

  ### projection all kin with deterministic initial conditions
  for(i in 1 : (na-1)){
    X_Focal[,i+1] <- G_tilde %*% X_Focal[,i]
    X_parents[, i+1] <- Uproj %*% X_parents[, i]
    X_younger_sibs[,i+1] <- Uproj %*% X_younger_sibs[,i] + Fprojstar %*% X_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj %*% X_younger_niece_nephew[,i] + Fproj %*% X_younger_sibs[,i]
    X_children[,i+1] <- Uproj %*% X_children[,i] + Fproj %*% X_Focal[,i]
    X_grand_children[,i+1] <- Uproj %*% X_grand_children[,i] + Fproj %*% X_children[,i]
    X_great_grand_children[,i+1] <- Uproj %*% X_great_grand_children[,i] + Fproj %*% X_grand_children[,i]
  }

  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews):
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  IC_f_great_grand_pars <- mothers_age_dist
  IC_m_great_grand_pars <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  IC_younger_aunts_uncles_f <- mothers_age_dist
  IC_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  for(ic in 1 : (na)){
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*X_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_great_grand_parents[, 1] <- X_great_grand_parents[, 1] + (IC_f_great_grand_pars[ic] + IC_m_great_grand_pars[ic])*X_grand_parents[,ic]
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*X_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*X_grand_children[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (IC_younger_aunts_uncles_f[ic] + IC_younger_aunts_uncles_m[ic])*X_younger_sibs[,ic]
  }

  ### Projections of grand parenst, older sibs, younger aunts/uncles, older nieces/nephews
  for(i in 1: (na-1)){
    X_grand_parents[, i+1] <- Uproj %*% X_grand_parents[, i]
    X_great_grand_parents[, i+1] <- Uproj %*% X_great_grand_parents[, i]
    X_older_sibs[,i+1] <- Uproj %*% X_older_sibs[,i]
    X_older_niece_nephew[,i+1] <- Uproj %*% X_older_niece_nephew[,i] + Fproj %*% X_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj %*% X_younger_aunts_uncles[,i] + Fprojstar %*% X_grand_parents[,i]
  }

  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  IC_older_aunt_uncle_f <- mothers_age_dist
  IC_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  for(ic in 1 : (na-1)){
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (IC_older_aunt_uncle_f[ic] + IC_older_aunt_uncle_m[ic])*X_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*X_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*X_younger_niece_nephew[,ic]
  }

  ## Projections of older unts/uncles, older cousins, younger cousins
  for(i in 1: (na-1)){
    X_older_aunt_uncle[,i+1] <- Uproj %*% X_older_aunt_uncle[,i]
    X_older_cousins[,i+1] <- Uproj %*% X_older_cousins[,i] + Fproj %*% X_older_aunt_uncle[,i]
    X_younger_cousins[,i+1] <- Uproj %*% X_younger_cousins[,i] + Fproj %*% X_younger_aunts_uncles[,i]
  }

  #### OUTPUT of all kin
  return(list("Focal" = X_Focal,
              "d" = X_children,
              "gd" = X_grand_children,
              "ggd" = X_great_grand_children,
              "m" = X_parents,
              "gm" = X_grand_parents,
              "ggm" = X_great_grand_parents,
              "os" = X_older_sibs,
              "ys" = X_younger_sibs,
              "nos" = X_older_niece_nephew,
              "nys" = X_younger_niece_nephew,
              "oa" = X_older_aunt_uncle,
              "ya" = X_younger_aunts_uncles,
              "coa" = X_older_cousins,
              "cya" = X_younger_cousins,
              "ps" = population_age_stage_structure
  ))
}


#' Title time-variant two-sex multi-state kin projections
#'
#' @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
#' @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
#' @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
#' @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
#' @param alpha scalar. birth ratio (male:female)
#' @param na scalar. number of ages.
#' @param ns scalar. number of stages.
#' @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
#' @param sex_Focal logical. Female or Male
#' @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
#' @param previous_kin_Focal matrix. last years kinship output.
#' @param prev_kin_children matrix. last years kinship output.
#' @param prev_kin_grandchildren matrix. last years kinship output.
#' @param prev_kin_greatgrandchildren matrix. last years kinship output.
#' @param prev_kin_parents matrix. last years kinship output.
#' @param prev_kin_grand_parents matrix. last years kinship output.
#' @param prev_kin_older_sibs matrix. last years kinship output.
#' @param prev_kin_younger_sibs matrix. last years kinship output.
#' @param prev_kin_older_niece_nephew matrix. last years kinship output.
#' @param prev_kin_younger_niece_nephew matrix. last years kinship output.
#' @param prev_kin_older_aunts_uncles matrix. last years kinship output.
#' @param prev_kin_younger_aunts_uncles matrix. last years kinship output.
#' @param prev_kin_older_cousins matrix. last years kinship output.
#' @param prev_kin_younger_cousins matrix. last years kinship output.
#' @param previous_population_age_stage_structure vector. The transient "population structure" (age*stage distributed)
#'
#' @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
#' nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
#' yielding the age*stage distribution of kin for each age of Focal
#'
all_kin_dy_TV <- function(Uf,
                          Um,
                          Ff,
                          Fm,
                          alpha, ## alpha = sex ratio male:female (i.e., 1 - birth_female)
                          na, ## number of ages
                          ns, ## number of stages
                          Parity,
                          sex_Focal,
                          Initial_stage_Focal,
                          previous_kin_Focal,
                          prev_kin_children,
                          prev_kin_grandchildren,
                          prev_kin_greatgrandchildren,
                          prev_kin_parents,
                          prev_kin_grand_parents,
                          prev_kin_great_grand_parents,
                          prev_kin_older_sibs,
                          prev_kin_younger_sibs,
                          prev_kin_older_niece_nephew,
                          prev_kin_younger_niece_nephew,
                          prev_kin_older_aunts_uncles,
                          prev_kin_younger_aunts_uncles,
                          prev_kin_older_cousins,
                          prev_kin_younger_cousins,
                          previous_population_age_stage_structure){

  n <- nrow(Uf)
  Uproj <- Matrix::Matrix(block_diag_function(list(Uf, Um)), sparse = TRUE)
  Fproj <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  Fproj[1:n, 1:n] <- (1-alpha)*Ff
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm
  Fprojstar <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE) ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff

  population_age_stage_structure <- previous_population_age_stage_structure
  population_age_stage_structure <- population_age_stage_structure/sum(population_age_stage_structure)
  population_age_stage_structure_next <- (Uproj + Fprojstar)%*%population_age_stage_structure

  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1

    population_age_stage_of_parenting <- pi_mix_TV_parity(Ff, Fm, alpha, na, ns, population_age_stage_structure)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]

    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]

  }
  else{

    population_age_stage_of_parenting <- pi_mix_TV(Ff, Fm, alpha, na, ns, population_age_stage_structure)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]

    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]

  }

  ### Focal requires its own dynamic: G_tilde constructed below tracks Focal's age*stage advancement over the time-scale
  f_t <- get_G(Uf, na, ns) ## get_G function in "Functions_required.R"
  m_t <- get_G(Um, na, ns)
  G_tilde <- block_diag_function(list(f_t,m_t))
  X_Focal  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }

  ### empty kin matrices for all of Focal's kin
  X_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_sibs <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_sibs <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_niece_nephew <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_niece_nephew <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_aunt_uncle <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_aunts_uncles <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_cousins <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_cousins <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)

  ### Initial distributions for kin with non-zero deterministic initial conditions:
  ## Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage
  ### projection all above kin with deterministic initial conditions
  for(i in 1 : (na-1)){
    X_Focal[,i+1] <- G_tilde %*% previous_kin_Focal[,i]
    X_parents[, i+1] <- Uproj %*% prev_kin_parents[, i]
    X_younger_sibs[,i+1] <- Uproj %*% prev_kin_younger_sibs[,i] + Fprojstar %*% prev_kin_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj %*% prev_kin_younger_niece_nephew[,i] + Fproj %*% prev_kin_younger_sibs[,i]
    X_children[,i+1] <- Uproj %*% prev_kin_children[,i] + Fproj %*% previous_kin_Focal[,i]
    X_grand_children[,i+1] <- Uproj %*% prev_kin_grandchildren[,i] + Fproj %*% prev_kin_children[,i]
    X_great_grand_children[,i+1] <- Uproj %*% prev_kin_greatgrandchildren[,i] + Fproj %*% prev_kin_grandchildren[,i]
  }

  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews):
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  IC_f_great_grand_pars <- mothers_age_dist
  IC_m_great_grand_pars <- fathers_age_dist
  IC_younger_aunts_uncles_f <- mothers_age_dist
  IC_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  for(ic in 1 : (na)){
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*prev_kin_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_great_grand_parents[, 1] <- X_great_grand_parents[, 1] + (IC_f_great_grand_pars[ic] + IC_m_great_grand_pars[ic])*prev_kin_grand_parents[,ic]
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*prev_kin_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*prev_kin_grandchildren[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (IC_younger_aunts_uncles_f[ic] + IC_younger_aunts_uncles_m[ic])*prev_kin_younger_sibs[,ic]
  }

  ### Projections of older sibs, younger aunts/uncles, older nieces/nephews
  for(i in 1: (na-1)){
    X_grand_parents[, i+1] <- Uproj %*% prev_kin_grand_parents[, i]
    X_great_grand_parents[, i+1] <- Uproj %*% prev_kin_great_grand_parents[, i]
    X_older_sibs[,i+1] <- Uproj %*% prev_kin_older_sibs[,i]
    X_older_niece_nephew[,i+1] <- Uproj %*% prev_kin_older_niece_nephew[,i] + Fproj %*% prev_kin_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj %*% prev_kin_younger_aunts_uncles[,i] + Fprojstar %*% prev_kin_grand_parents[,i]
  }

  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  IC_older_aunt_uncle_f <- mothers_age_dist
  IC_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  for(ic in 1 : (na-1)){
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (IC_older_aunt_uncle_f[ic] + IC_older_aunt_uncle_m[ic])*prev_kin_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*prev_kin_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*prev_kin_younger_niece_nephew[,ic]
  }

  ## Projections of older unts/uncles, older cousins, younger cousins
  for(i in 1: (na-1)){
    X_older_aunt_uncle[,i+1] <- Uproj %*% prev_kin_older_aunts_uncles[,i]
    X_older_cousins[,i+1] <- Uproj %*% prev_kin_older_cousins[,i] + Fproj %*% prev_kin_older_aunts_uncles[,i]
    X_younger_cousins[,i+1] <- Uproj %*% prev_kin_younger_cousins[,i] + Fproj %*% prev_kin_younger_aunts_uncles[,i]
  }

  return(list("Focal" = X_Focal,
              "d" = X_children,
              "gd" = X_grand_children,
              "ggd" = X_great_grand_children,
              "m" = X_parents,
              "gm" = X_grand_parents,
              "ggm" = X_great_grand_parents,
              "os" = X_older_sibs,
              "ys" = X_younger_sibs,
              "nos" = X_older_niece_nephew,
              "nys" = X_younger_niece_nephew,
              "oa" = X_older_aunt_uncle,
              "ya" = X_younger_aunts_uncles,
              "coa" = X_older_cousins,
              "cya" = X_younger_cousins,
              "ps" = population_age_stage_structure_next))
}

################## Create data frame output

## Use of "pipe" (don't understand the name, but hey)
`%>%` <- magrittr::`%>%`

#' Title Accumulated kin by each age of Focal, for each time period, and cohort of birth
#'
#' @param kin_matrix_lists list of lists of kin matrices: list( list(X_focal), list(X_parents), ... ). Outer list is length 14  = number of kin. Inner lists have lenght = timescale
#'                 so list(X_focal) = list(X_focal[year1],X_focal[year2],...,X_focal[yearlast])
#' @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
#' @param years vector. The timescale on which we implement the kinship model.
#' @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
#' @param na numeric. Number of ages.
#' @param ns numeric. Number of stages.
#' @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#'
#' @return A data frame which gives for each age of Focal at each year in the timescale, Focal's experienced number kin demarcated by stages (summed over all ages)
#'
create_cumsum_df <- function(kin_matrix_lists,
                             kin_names,
                             years,
                             model_years,
                             na,
                             ns,
                             specific_kin = NULL,
                             increment = NULL){
  if(length(years) > length(kin_matrix_lists[[1]])){stop("More years than data")}

  matrix_model_time <- model_years
  age_inc <- increment

  df_year_list <- list()
  for(j in years){
    ii <- which(matrix_model_time == j)
    df_list <- list()
    for(i in 1 : length(kin_names)){
      kin_member <- kin_names[[i]]
      kin_data <- kin_matrix_lists[[i]]
      kin_data <- kin_data[[ii]]
      df <- as.data.frame(as.matrix(kin_data))
      dims <- dim( kin_data )
      nr <- dims[1]
      nc <- dims[2]
      female_kin <- df[1:(nr/2), 1:nc]
      male_kin <- df[ (1+nr/2) : nr, 1:nc]
      colnames(female_kin) <- seq( 0 , age_inc*(na-1) , by = age_inc )
      colnames(male_kin) <- seq( 0 , age_inc*(na-1) , by = age_inc )
      female_kin$stage <- rep(seq(1, ns), na)
      male_kin$stage <- rep(seq(1, ns), na)
      female_kin$age <- rep( seq( 0 , age_inc*(na-1) , by = age_inc ) , each = ns)
      male_kin$age <- rep( seq( 0 , age_inc*(na-1) , by = age_inc ) , each = ns)
      female_kin$Sex <- "Female"
      male_kin$Sex <- "Male"
      both_kin <- rbind(female_kin, male_kin)
      both_kin <- both_kin %>% reshape2::melt(id = c("age","stage","Sex")) %>%
        dplyr::group_by(variable, stage, Sex) %>%
        dplyr::summarise(num = sum(value)) %>%
        dplyr::ungroup()
      both_kin <- both_kin %>% dplyr::transmute(age_focal = variable,
                                                stage_kin = as.factor(stage),
                                                count = num,
                                                sex_kin = Sex)
      both_kin$age_focal <- as.numeric(paste(both_kin$age_focal))
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(age_focal),
                                                 cohort_factor = as.factor(cohort))
  if(!is.null(specific_kin)){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}

#' Title joint age*stage distributions of kin by each age of Focal, for each time period, and cohort of birth
#'
#' @param kin_matrix_lists list of lists of kin matrices: list( list(X_focal), list(X_parents), ... ). Outer list is length 14  = number of kin. Inner lists have lenght = timescale
#'                 so list(X_focal) = list(X_focal[year1],X_focal[year2],...,X_focal[yearlast])
#' @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
#' @param years vector. The timescale on which we implement the kinship model.
#' @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
#' @param na numeric. Number of ages.
#' @param ns numeric. Number of stages.
#' @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#'
#' @return A data frame which gives for each age of Focal at each year in the timescale, the full age*stage dist of kin
#'
create_full_dists_df <- function(kin_matrix_lists,
                                 kin_names,
                                 years,
                                 model_years,
                                 na,
                                 ns,
                                 specific_kin = NULL,
                                 increment = NULL){
  if(length(years) > length(kin_matrix_lists[[1]])){stop("More years than data")}

  matrix_model_time <- model_years
  age_inc <- increment

  df_year_list <- list()
  for(j in years){
    df_list <- list()
    ii <- which(matrix_model_time == j)
    for(i in 1 : length(kin_names)){
      kin_member <- kin_names[[i]]
      kin_data <- kin_matrix_lists[[i]]
      kin_data <- kin_data[[ii]]
      df <- as.data.frame(as.matrix(kin_data))
      dims <- dim( kin_data)
      nr <- dims[1]
      nc <- dims[2]
      female_kin <- df[1:(nr/2), 1:nc]
      male_kin <- df[ (1+nr/2) : nr, 1:nc]
      colnames(female_kin) <- seq( 0 , age_inc*(na-1) , by = age_inc )
      colnames(male_kin) <- seq( 0 , age_inc*(na-1) , by = age_inc )
      female_kin$stage <- rep(seq(1, ns), na)
      male_kin$stage <- rep(seq(1, ns), na)
      female_kin$age <- rep( seq( 0 , age_inc*(na-1) , by = age_inc ) , each = ns)
      male_kin$age <- rep( seq( 0 , age_inc*(na-1) , by = age_inc ) , each = ns)
      female_kin$Sex <- "Female"
      male_kin$Sex <- "Male"
      both_kin <- rbind(female_kin, male_kin)
      both_kin <- both_kin %>% reshape2::melt(id = c("age","stage","Sex")) %>%
        dplyr::transmute(age_focal = variable,
                         age_kin = age,
                         stage_kin = as.factor(stage),
                         count = value,
                         sex_kin = Sex)
      both_kin$age_focal <- as.numeric(paste(both_kin$age_focal))
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list

  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(age_focal),
                                                 cohort_factor = as.factor(cohort))
  if(!is.null(specific_kin)){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}



## Construct a matrix composed as a direct sum of a list of matrices
block_diag_function <- function(mat_list){
  s = length(mat_list)
  u1 = mat_list[[1]]
  dims <- dim(u1)
  r = dims[1]
  diagmat <- Matrix::Matrix(nrow = (r*s), ncol = (r*s), data = 0, sparse = TRUE)
  for(i in 1:s){
    diagmat = diagmat + kronecker(E_matrix(i,i,s,s), mat_list[[i]])
  }
  return(diagmat)
}

## Construct a matrix which transfers Focal across stages, while ensuring Focal survives with probability 1
get_G <- function(U, na, ns){
  sig <- Matrix::t(rep(1,na*ns)) %*% U
  diag <- Matrix::diag(sig[1,])
  G <- U %*% MASS::ginv(diag)
  return(G)
}

#' Mixing distributions for the time-invariant multi-state 2-sex model: Non-parity case
#'
#' @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
#' @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
#' @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
#' @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
#' @param alpha scalar. Birth ratio male:female
#' @param na scalar. Number of age-classes
#' @param ns scalar. Number of stages
#'
#' @return list (of vectors). list[[1]] = full age*stage*sex distribution
#'                            list[[2]] = female age*stage distribution normalised
#'                            list[[3]] = male age*stage distribution normalised
#'                            list[[4]] = female marginal age distribution normalised
#'                            list[[5]] = male marginal age distribution normalised
#'
pi_mix <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

#' Mixing distributions for the time-variant multi-state 2-sex model: Non-parity case
#'
#' @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
#' @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
#' @param alpha scalar. Birth ratio male:female
#' @param na scalar. Number of age-classes
#' @param ns scalar. Number of stages
#' @param previous_age_stage_dist vector. Last years population structure (age*stage*sex full distribution)
#'
#' @return list (of vectors). list[[1]] = full age*stage*sex distribution
#'                            list[[2]] = female age*stage distribution normalised
#'                            list[[3]] = male age*stage distribution normalised
#'                            list[[4]] = female marginal age distribution normalised
#'                            list[[5]] = male marginal age distribution normalised
#'
pi_mix_TV <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

#' Mixing distributions for the time-invariant multi-state 2-sex model: Parity-specific case
#'
#' @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
#' @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
#' @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
#' @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
#' @param alpha scalar. Birth ratio male:female
#' @param na scalar. Number of age-classes
#' @param ns scalar. Number of stages
#'
#' @return list (of vectors). list[[1]] = full age*stage*sex distribution
#'                            list[[2]] = female age*stage distribution normalised
#'                            list[[3]] = male age*stage distribution normalised
#'                            list[[4]] = female marginal age distribution normalised
#'                            list[[5]] = male marginal age distribution normalised
#'
pi_mix_parity <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  for(i in 1:na){
    m_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% m_mat[,i]
    d_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% d_mat[,i]
  }
  out_mum <- m_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(m_mat)))
  out_dad <- d_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(d_mat)))
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(c(pi_f,pi_m), pi_f, pi_m, pi_F, pi_M))
}

#' Mixing distributions for the time-variant multi-state 2-sex model: Parity-specific case
#'
#' @param Uf matrix. Block-structured matrix which transfers females over stage and advances their age
#' @param Um matrix. Block-structured matrix which transfers males over stage and advances their age
#' @param Ff matrix. Block-structured matrix which counts reproduction by females and assigns newborns an age and stage
#' @param Fm matrix. Block-structured matrix which counts reproduction by males and assigns newborns an age and stage
#' @param alpha scalar. Birth ratio male:female
#' @param na scalar. Number of age-classes
#' @param ns scalar. Number of stages
#' @param previous_age_stage_dist vector. Last years population structure (age*stage*sex full distribution)
#'
#' @return list (of vectors). list[[1]] = full age*stage*sex distribution
#'                            list[[2]] = female age*stage distribution normalised
#'                            list[[3]] = male age*stage distribution normalised
#'                            list[[4]] = female marginal age distribution normalised
#'                            list[[5]] = male marginal age distribution normalised
#'
pi_mix_TV_parity <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F  <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  for(i in 1:na){
    m_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% m_mat[,i]
    d_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% d_mat[,i]
  }
  out_mum <- m_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(m_mat)))
  out_dad <- d_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(d_mat)))
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}


######################################################### Some useful utility functions required


###################################################### Eigen-decomposition of a matrix

# Calculate the spectral radius of a matrix (growth rate in Demographics)
lambda <- function(PM) {
  lead_eig <- (abs(eigen(PM, only.values = TRUE)$values))
  lead_eig <- lead_eig[which.max(lead_eig)]
  return(lead_eig)
}
# Find the column-eigenvector corresponding to the spectral radius (Stable population structure in Demographics)
SD <- function(PM) {
  spectral_stuff <- eigen(PM)
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)
}
# Find the row-eigenvector corresponding to the spectral radius (Stable reproductive values in Demographics)
RD <- function(PM) {
  spectral_stuff <- eigen(t(PM))
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)
}

###################################################### Useful matrix operations

## Constructing a unit vector with a 1 in the ith position
e_vector <- function(i, n){
  e <- rep(0, n)
  e[i] <- 1
  return(e)
}
## Creating a matrix of zeros with a 1 in the i,j-th entry
E_matrix <- function(i,j,n,m){
  E <- Matrix::Matrix(nrow = (n), ncol = (m), data = 0, sparse = TRUE)
  E[i,j] <- 1
  return(E)

}
## Creating the Vec-commutation matrix
K_perm_mat <- function(n,m){
  perm <- Matrix::Matrix(nrow = (n*m), ncol = (n*m), data = 0, sparse = TRUE)
  for(i in 1:n){
    for(j in 1:m){
      perm = perm + kronecker( E_matrix(i,j,n,m) , Matrix::t(E_matrix(i,j,n,m)) )
    }
  }
  return(perm)
}









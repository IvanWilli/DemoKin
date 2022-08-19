# replicating Caswell´s figures: choose some kin
library(DemoKin)
library(tidyverse)
library(R.matlab)
source("R/kin_time_invariant.R")

# paper input from https://www.demographic-research.org/volumes/vol45/16/45-16.pdf
input_time_variant <- readMat("tests/SWEhist_matrices.mat")

# check structure from reading mat
class(input_time_variant)
names(input_time_variant)
length(input_time_variant[["matrices"]]) # number of years
input_time_variant[["matrices"]][[128]][[1]][[1]] # U
input_time_variant[["matrices"]][[1]][[1]][[2]] # F
input_time_variant[["matrices"]][[1]][[1]][[3]] # popsize
input_time_variant[["matrices"]][[1]][[1]][[4]] # pi
length(input_time_variant_proj[["matrices"]]) # number of years

# reshape
U_hal <- f_hal <-N_hal <- pi_hal <-matrix(rep(0,111))
for(y in 1:128){
  U <- input_time_variant[["matrices"]][[y]][[1]][[1]] %>% as.matrix()
  f <- input_time_variant[["matrices"]][[y]][[1]][[2]] %>% as.matrix()
  N <- input_time_variant[["matrices"]][[y]][[1]][[3]] %>% as.matrix()
  pi <- input_time_variant[["matrices"]][[y]][[1]][[4]] %>% as.matrix()
  U_hal <- cbind(U_hal, c(U[col(U)==row(U)-1], U[ncol(U),nrow(U)]))
  f_hal <- cbind(f_hal ,f[1,])
  N_hal <- cbind(N_hal ,N)
  pi_hal <-cbind(pi_hal, pi)
}
U_hal <- U_hal[,-1]
f_hal <- f_hal[,-1]
N_hal <- N_hal[,-1]
pi_hal <-pi_hal[,-1]
colnames(U_hal) <- colnames(f_hal) <- colnames(N_hal) <- colnames(pi_hal) <-1891:2018
dim(U_hal);class(U_hal %>% as.matrix)

# output from Hal (dropbox link https://www.dropbox.com/t/3YiILmn7SpczN3oM)
output_time_variant <- readMat("tests/time-varying_sweden.mat")

# inspect the way the package reads mat
class(output_time_variant)
names(output_time_variant)
length(output_time_variant[["allkin"]]) # number of years
length(output_time_variant[["allkin"]][[1]])
length(output_time_variant[["allkin"]][[1]])
class(output_time_variant[["allkin"]][[1]][[1]]) # 1 array with kin matrix
dim(output_time_variant[["allkin"]][[1]][[1]][,,14]) # the matrix of the nth kin, 111 ages

# use own codes to interpret
codes <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
caswell_codes <- c("t", "v", "a", "b", "c", "h", "g", "d", "p", "q", "r", "s", "m", "n")

# re shape data to tidy
output_time_variant_df <- map_df(1:128, function(i){
  array_branch(output_time_variant[["allkin"]][[i]][[1]], margin = 3) %>%
    map_df(., as.data.frame)}) %>%
  setNames(as.character(0:110)) %>%
  bind_cols(crossing(year = 1891+(0:127), # years
                     kin_index = 1:14, # number of possible kin
                     age_kin = 0:110) # ages
            ) %>%
  inner_join(tibble(kin = codes, caswell_codes) %>%
               arrange(caswell_codes) %>% mutate(kin_index = 1:14))

# check dimension: 128 years, 14 types of kin, 111 ages
nrow(output_time_variant_df); 128*14*111

# own calculation for first year
out_first_year <- kin_time_invariant(
                      U = U_hal[,"1891"],
                      f = f_hal[,"1891"],
                      pi = pi_hal[,"1891"],
                      birth_female = 1)

# check first visually demokin
out_first_year %>%
  filter(alive == "yes") %>%
  group_by(age_focal, kin) %>%
  summarise(count = sum(count, na.rm=T)) %>%
  ggplot(aes(age_focal, count)) +
  geom_line() +
  facet_wrap(~kin, scales="free_y")

# compare with paper results
comparison <- out_first_year %>%
  filter(alive == "yes") %>%
  group_by(age_focal, kin) %>%
  summarise(count = sum(count, na.rm=T)) %>%
  mutate(source = "demokin") %>%
  bind_rows(
    output_time_variant_df %>%
      filter(year %in% 1891) %>%
      pivot_longer(`0`:`110`, names_to = "age", values_to = "count") %>%
      mutate(age = as.integer(age)) %>%
      group_by(age_focal=age, kin) %>%
      summarise(count = sum(count)) %>%
      mutate(source = "paper"))

# comparison visually
comparison %>%
  ggplot() +
  geom_line(aes(age_focal, count, color=source, linetype=source)) +
  facet_wrap(~kin, scales="free_y") +
  theme_bw()

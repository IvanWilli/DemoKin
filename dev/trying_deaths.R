library(devtools)
library(tidyverse)
load_all()
# debugonce(kins)

ego_age = 35
year = 2015
P = swe_surv
asfr = swe_asfr
N = swe_pop
age = 0:100
birth_female = 1/2.04

start_time = Sys.time()
swe35_2015_stable <- kins(ego_age = 35, year = 2015,
                          P = swe_surv, asfr = swe_asfr, N=swe_pop,
                          stable = FALSE, alive = c("asdfsd"))
end_time = Sys.time()
end_time-start_time
(end_time-start_time)/60

D <- swe35_2015_stable[["kins_death"]]
L <- swe35_2015_stable[["kins_living"]]

L[["kins_by_age_kin"]] %>%
  select(-x) %>%
  gather(kin, count, -x_kin) %>%
  ggplot() +
  geom_line(aes(x_kin, count))  +
  geom_vline(xintercept = 35, color=2)+
  theme_bw() +
  facet_wrap(~kin)

L[["kins_by_age_ego"]] %>%
  gather(kin, count, -x) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = 35, color=2)+
  theme_bw() +
  facet_wrap(~kin)


D$kins_cum_death_total

D$kins_death_by_age_ego %>%
  gather(kin, count, -x,) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = 35, color=2)+
  theme_bw() +
  facet_wrap(~kin,scales="free")

D$kins_cum_death_by_age_ego %>%
  gather(kin, count, -x,) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = 35, color=2)+
  theme_bw() +
  facet_wrap(~kin,scales="free")

D$lost_mean_age

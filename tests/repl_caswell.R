# replicating Caswell´s figures: choose some kin

library(devtools)
load_all()
library(DemoKin)
library(tidyverse)

#time invariant

# Caswell's assumption on time invariant !!!!!!!!!!!!!!!!!!!!!!
ft[1,1:ages] = f * U * birth_female

swe_surv_2010 <- swe_surv %>% pull(`2011`)
swe_asfr_2010 <- swe_asfr %>% pull(`2011`)
debugonce(kin)
swe50_2015_stable <- kin(U = swe_surv_2010, f = swe_asfr_2010, output_cohort = c(1911,1930),
                         output_kin = c("d","m"))

swe_kin_cohorts <- kin(U = U_caswell_2021, f = f_caswell_2021, time_invariant = F,
                       birth_female = 1,
                       output_cohort = c(1911),
                       output_kin = c("d"))

U = U_caswell_2021; f = f_caswell_2021; pi = pi_caswell_2021; birth_female = 1;
output_cohort = c(1911);output_period = NULL; output_kin = c("d")

# FIGURE 5

# cohort
swe_kin_cohorts <- kin(U = U_caswell_2021, f = f_caswell_2021, pi = pi_caswell_2021, time_invariant = F,
                       birth_female = 1,
                output_cohort = c(1891,1911,1931,1951,1971,2011,2041),
                output_kin = c("d","gd","ggd","m","gm","ggm"))

swe_kin_cohorts$kin_summary %>%
  ggplot(aes(year,count,color=factor(cohort))) +
  scale_y_continuous(name = "",labels = seq(0,3,.2),breaks = seq(0,3,.2))+
  geom_line(size=1)+
  geom_vline(xintercept = 2019, linetype=2)+
  facet_wrap(~kin,scales = "free")+
  theme_bw()

# period
swe_kin_period <- kin(U = U_caswell_2021, f = f_caswell_2021, pi = pi_caswell_2021, stable = F, birth_female = 1,
                 focal_year = c(1891,1921,1951,2010,2050,2080,2120),
                 selected_kin = c("d","gd","ggd","m","gm","ggm","os","ys","oa","ya"))

swe_kin_period$kin_summary %>%
  ggplot(aes(age_focal,count,color=factor(year))) +
  geom_line(size=1)+
  scale_y_continuous(name = "",labels = seq(0,3,.2),breaks = seq(0,3,.2))+
  facet_wrap(~kin, scales = "free")+
  theme_bw()

# ADDITIONAL PLOTS cohrot and period
ggplot(swe_kin_cohorts$kin_summary %>% filter(cohort == 1911),
       aes(year,mean_age)) +
  geom_point(aes(size=count,color=kin)) +
  geom_line(aes(color=kin)) +
  scale_y_continuous(name = "Edad", breaks = seq(0,110,10), labels = seq(0,110,10), limits = c(0,110))+
  geom_segment(x = 1911, y = 0, xend = 2025, yend = 110, color = 1)+
  geom_vline(xintercept = 1911, linetype=2)+
  theme_light()+ coord_fixed()+
  labs(title = "Kin cohort 1911")

swe_kin_period$kin_summary %>%
  filter(age_focal==50) %>%
  ggplot(aes(year, mean_age, color=kin)) +
  geom_point(aes(size=count)) +
  geom_line() +
  geom_hline(yintercept = 50, color=1, linetype=1)+
  theme_light()+
  coord_fixed()+
  labs(title = "Kin period")


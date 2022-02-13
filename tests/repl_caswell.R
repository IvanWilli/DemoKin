# replicating Caswell´s figures: choose some kins

# load input data
# devtools::load_all()
# load("data/swe_caswell.rda")
library(DemoKin)

# FIGURE 5

# cohort
swe_kins_cohorts <- kins(U = U, f = f, pi = pi, stable = F,birth_female = 1,
                ego_cohort = c(1891,1911,1931,1951,1971,2011,2041),
                selected_kins = c("d","gd","ggd","m","gm","ggm"))
swe_kins_cohorts$kins_by_age_ego %>%
  ggplot(aes(year,total,color=factor(cohort))) +
  scale_y_continuous(name = "",labels = seq(0,3,.2),breaks = seq(0,3,.2))+
  geom_line(size=1)+
  geom_vline(xintercept = 2019, linetype=2)+
  facet_wrap(~kin,scales = "free")+
  theme_bw()

# period
swe_kins_period <- kins(U = U, f = f, pi = pi, stable = F, birth_female = 1,
                 ego_year = c(1891,1921,1951,2010,2050,2080,2120),
                 selected_kins = c("d","gd","ggd","m","gm","ggm","os","ys","oa","ya"))

swe_kins_period$kins_by_age_ego %>%
  ggplot(aes(age_ego,total,color=factor(year))) +
  geom_line(size=1)+
  scale_y_continuous(name = "",labels = seq(0,3,.2),breaks = seq(0,3,.2))+
  facet_wrap(~kin, scales = "free")+
  theme_bw()

# ADDITIONAL PLOTS cohrot and period
ggplot(swe_kins_cohorts$kins_by_age_ego %>% filter(cohort == 1911),
       aes(year,mean_age)) +
  geom_point(aes(size=total,color=kin)) +
  geom_line(aes(color=kin)) +
  scale_y_continuous(name = "Edad", breaks = seq(0,110,10), labels = seq(0,110,10), limits = c(0,110))+
  geom_segment(x = 1911, y = 0, xend = 2025, yend = 110, color = 1)+
  geom_vline(xintercept = 1911, linetype=2)+
  theme_light()+ coord_fixed()+
  labs(title = "Kin cohort 1911")

swe_kins_period$kins_by_age_ego %>%
  filter(age_ego==50) %>%
  ggplot(aes(year, mean_age, color=kin)) +
  geom_point(aes(size=total)) +
  geom_line() +
  geom_hline(yintercept = 50, color=1, linetype=1)+
  theme_light()+
  coord_fixed()+
  labs(title = "Kin period")


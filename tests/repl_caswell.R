# replicating Caswell´s figures: choose some kin

library(devtools)
load_all()
library(DemoKin)
library(tidyverse)
library(progress)
library(R.matlab)
load("tests/test.RData")

# basic
swe_kin_period_pack <- kin(U = swe_surv,
                      f = swe_asfr,
                      N = swe_pop,
                      time_invariant = F,
                      birth_female = 1,
                      output_period = c(1900, 1950, 2010),
                      output_kin = c("d","gd","m","gm","oa", "os"))

swe_kin_period_pack$kin_full %>%
  filter(alive == "yes") %>%
  group_by(age_focal, kin, year) %>%
  summarise(count = sum(count, na.rm=T)) %>%
  ggplot(aes(age_focal, count, color=factor(year))) +
  geom_line() +
  facet_wrap(~kin, scales="free_y")

# time variant ------------------------------------------------------------

# inputs
input_time_variant <- readMat("tests/SWEhist_matrices.mat")
input_time_variant_proj <- readMat("tests/SWEproj_matrices.mat")
# class(input_time_variant)
# names(input_time_variant)
# length(input_time_variant[["matrices"]]) # number of years
# input_time_variant[["matrices"]][[128]][[1]][[1]] # U
# input_time_variant[["matrices"]][[1]][[1]][[2]] # F
# input_time_variant[["matrices"]][[1]][[1]][[3]] # popsize
# input_time_variant[["matrices"]][[1]][[1]][[4]] # pi
# length(input_time_variant_proj[["matrices"]]) # number of years

U_hal <- f_hal <-N_hal <- pi_hal <-matrix(rep(0,111))
for(y in 1:128){
  # y = 1
  U <- input_time_variant[["matrices"]][[y]][[1]][[1]] %>% as.matrix()
  f <- input_time_variant[["matrices"]][[y]][[1]][[2]] %>% as.matrix()
  N <- input_time_variant[["matrices"]][[y]][[1]][[3]] %>% as.matrix()
  pi <- input_time_variant[["matrices"]][[y]][[1]][[4]] %>% as.matrix()
  U_hal <- cbind(U_hal, c(U[col(U)==row(U)-1], U[ncol(U),nrow(U)]))
  f_hal <- cbind(f_hal ,f[1,])
  N_hal <- cbind(N_hal ,N)
  pi_hal <-cbind(pi_hal, pi)
}
U_hal_end <- U_hal[,-1]
f_hal_end <- f_hal[,-1]
N_hal_end <- N_hal[,-1]
pi_hal_end <-pi_hal[,-1]
colnames(U_hal_end) <- colnames(f_hal_end) <- colnames(N_hal_end) <- colnames(pi_hal_end) <-1891:2018
dim(U_hal_end);class(U_hal_end %>% as.matrix)

# period
swe_kin_period <- kin(U = U_hal_end %>% as.matrix(),
                      f = f_hal_end %>% as.matrix(),
                      pi = pi_hal_end %>% as.matrix(),
                      time_invariant = F,
                      birth_female = 1,
                      output_period = c(1891,1921,1951,2010),
                      output_kin = c("d","gd","m","gm","oa", "os"))

# check first-row plots from figures 5-A  and 5-B from https://www.demographic-research.org/volumes/vol45/16/45-16.pdf
swe_kin_period$kin_full %>%
  filter(alive == "yes") %>%
  group_by(age_focal, kin, year) %>%
  summarise(count = sum(count, na.rm=T)) %>%
  ggplot(aes(age_focal, count, color=factor(year))) +
  geom_line() +
  facet_wrap(~kin, scales="free_y")

# read from https://www.dropbox.com/t/3YiILmn7SpczN3oM
output_time_variant <- readMat("tests/time-varying_sweden.mat")

# inspect the way the package reads
# class(output_time_variant)
# names(output_time_variant)
# length(output_time_variant[["allkin"]]) # number of years
# length(output_time_variant[["allkin"]][[1]])
# length(output_time_variant[["allkin"]][[1]])
# class(output_time_variant[["allkin"]][[1]][[1]]) # 1 array with kin matrix
# dim(output_time_variant[["allkin"]][[1]][[1]][,,14]) # the matrix of the nth kin, 111 ages

# use own codes to interpret
codes <- c("coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
caswell_codes <- c("t", "v", "a", "b", "c", "h", "g", "d", "p", "q", "r", "s", "m", "n")

# re arrange all data to a dataframe
output_time_variant_df <- map_df(1:128, function(i){
  array_branch(output_time_variant[["allkin"]][[i]][[1]], margin = 3) %>%
    map_df(., as.data.frame)}) %>%
  setNames(as.character(0:110)) %>%
  bind_cols(crossing(year = 1891+(0:127),
                     kin_index = 1:14,
                     age_kin = 0:110)) %>%
  inner_join(tibble(kin = codes, caswell_codes) %>%
               arrange(caswell_codes) %>% mutate(kin_index = 1:14))

# check dimension: 128 years, 14 types of kin, 111 ages
nrow(output_time_variant_df); 128*14*111

# check first-row plots from figures 5-A  and 5-B from https://www.demographic-research.org/volumes/vol45/16/45-16.pdf
output_time_variant_df %>%
  filter(year %in% c(1891,1921,1951,2010), kin %in% c("d","gd", "m", "gm", "oa", "os")) %>%
  pivot_longer(`0`:`110`, names_to = "age", values_to = "count") %>%
  mutate(age = as.integer(age)) %>%
  group_by(age, kin, year) %>%
  summarise(count = sum(count)) %>%
  ggplot(aes(age, count, color=factor(year))) +
  geom_line() +
  facet_wrap(~kin, scales="free_y")

# differences - look d, gd, in 1891 and 1951
swe_period_together <- swe_kin_period$kin_full %>%
  filter(alive == "yes") %>%
  filter(year %in% c(1891,1921,1951,2010), kin %in% c("d","gd","m","gm","oa", "os")) %>%
  group_by(age_focal, kin, year) %>% summarise(count_demokin = sum(count, na.rm=T)) %>%
  inner_join(
    output_time_variant_df %>%
      filter(year %in% c(1891,1921,1951,2010), kin %in% c("d","gd", "m", "gm","oa", "os")) %>%
      pivot_longer(`0`:`110`, names_to = "age", values_to = "count") %>%
      mutate(age = as.integer(age)) %>%
      group_by(age_focal=age, kin, year) %>%
      summarise(count_paper = sum(count)))

swe_period_together %>%
  filter(year == 1891) %>%
  ggplot() +
  geom_line(aes(age_focal, count_demokin, color=factor(year)), linetype=1) +
  geom_line(aes(age_focal, count_paper, color=factor(year)), linetype=2) +
  facet_wrap(~kin, scales="free_y")

swe_period_rel_dif <- swe_period_together %>%
  mutate(rel_dif = round(100*(count_paper/count_demokin-1),3)) %>%
  arrange(year, kin) %>%
  as.data.frame() %>%
  group_by(year, kin) %>% summarise(sum(rel_dif, na.rm=T))














# to bind projected
# U_hal <- U_hal[1:106,]
# f_hal <- f_hal[1:106,]
# N_hal <- N_hal[1:106,]
# pi_hal <-pi_hal[1:106,]
# for(y in 1:102){
#   # y = 1
#   U <- input_time_variant_proj[["matrices"]][[y]][[1]][[1]]
#   f <- input_time_variant_proj[["matrices"]][[y]][[1]][[2]]
#   N <- input_time_variant_proj[["matrices"]][[y]][[1]][[3]]
#   pi <- input_time_variant_proj[["matrices"]][[y]][[1]][[4]]
#   U_hal <- U_hal %>% bind_cols(c(U[col(U)==row(U)-1], U[ncol(U),nrow(U)]))
#   f_hal <- f_hal %>% bind_cols(f[1,])
#   N_hal <- N_hal %>% bind_cols(N)
#   pi_hal <-pi_hal%>% bind_cols(as.numeric(pi))
# }
# dim(U_hal[,-1])
# U_hal_end <- U_hal[,-1] %>% setNames(as.character(1891:2120))
# f_hal_end <- f_hal[,-1] %>% setNames(as.character(1891:2120))
# N_hal_end <- N_hal[,-1] %>% setNames(as.character(1891:2120))
# pi_hal_end <-pi_hal[,-1] %>% setNames(as.character(1891:2120))
# dim(U_hal_end);names(U_hal_end)

# time invariant ----------------------------------------------------------

### data: survival probability and fertility by age for Japan
# available at https://www.demographic-research.org/volumes/vol41/24/default.htm

p_1947 <- 1 - read.csv("tests/qx_years.csv", header = F, sep = " ")[[4]]
f_1947 <- read.csv("tests/fx_years.csv", header = F, sep = " ")[[4]]
p_2014 <- 1 - read.csv("tests/qx_years.csv", header = F, sep = " ")[[205]]
f_2014 <- read.csv("tests/fx_years.csv", header = F, sep = " ")[[205]]

# Caswell assumption on first age
f_1947 <- f_1947 * p_1947
f_2014 <- f_2014 * p_2014

kins_japan_1947 <- kin(p_1947, f_1947, living = F)$kin_full
kins_japan_1947 %>%
  filter(alive=="yes", kin=="ggm") %>%
  group_by(age_focal) %>% summarise(sum(count))


### results
kins_japan <- rbind(tibble(Year = 1947, kin(p_1947, f_1947, living = F)$kin_full),
                    tibble(Year = 2014, kin(p_2014, f_2014, living = F)$kin_full))

# kins alive by age when ego is aged 30 or 70
kins_japan %>%
  filter(age_focal %in% c(30,70), alive=="yes") %>%
  ggplot() +
  geom_line(aes(x=age_kin, y=count,
                color=factor(age_focal), linetype=factor(Year)))  +
  facet_wrap(~kin,scales = "free_y") +
  theme_classic() +
  facet_wrap(~kin,scales = "free_y")

kins_japan %>%
  filter(age_focal %in% 30, alive=="yes", kin == "m", Year==2014)

### get paper results: done with https://plotdigitizer.com/app

m_30_2014 <- c(48.124993716677295, 0.0068724848600042006,
               52.13541022398433, 0.022765097394085585,
               56.14582673129136, 0.056697985757917374,
               60.04166165822103, 0.07398657677157613,
               64.16665974590543, 0.054765100671140945,
               68.17707625321246, 0.02330201014576342,
               71.95832959976478, 0.0035436192454907658) %>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "m", year = 2014, age = 30)
m_30_1947 <- c(47.89583055592257, 0.010630874121749168,
               51.791665482852245, 0.045100671140939595,
               56.37499863406029, 0.05111409232120386,
               59.92708007784368, 0.03908724586435613,
               63.82291500477335, 0.02577181208053692,
               68.17707625321246, 0.012671136024014266,
               71.84374801938742, 0.0025771828465813683) %>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "m", year = 1947, age = 30)
s_30_1947 <- c(117.30467373388943, 0.028030307190039253,
               3.8671896129379024, 0.0011363654972982584,
               8.05663990468026, 0.005681817853417949,
               12.031249743886312, 0.01792929767285178,
               16.220700035628674, 0.036111111913867226,
               19.873045098211307, 0.05782828333912762,
               23.84765903523633, 0.07626262618964844,
               28.037105229159724, 0.08244949644554042,
               32.119139392452894, 0.06717171906914071,
               36.09374513383999, 0.044696973856705915,
               39.853510422690775, 0.024621210698144477,
               43.93554458598395, 0.010858592937435215,
               47.80273010110288, 0.00303030478177092) %>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "s", year = 1947, age = 30)
s_70_1947 <- c(43.93554458598395, 0.0011363654972982584,
               47.80273010110288, 0.00441919166376952,
               52.20702494320051, 0.013383845316732092,
               56.074218653957374, 0.026388889290266948,
               59.94140416907631, 0.03952020358922534,
               64.02343013673156, 0.048358586916764375,
               68.10546430002474, 0.045959600046354354,
               71.97264981514367, 0.03143939886539736,
               76.05468397843684, 0.015404045293554923,
               80.02928971982392, 0.005429296468717607,
               84.00390365684893, 0.0011363654972982584) %>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "s", year = 1947, age = 70)
s_30_2014 <- c(7.949219678412107, 0.0005050524024740438,
               12.031249743886312, 0.0032828357995446046,
               16.005859583092366, 0.009217175037662914,
               19.873045098211307, 0.020328289359798468,
               23.6328103870621, 0.03194444645133473,
               27.822264776623417, 0.03888889049440111,
               32.01171916618474, 0.0337121250434572,
               35.77148445503553, 0.0215909155494469,
               39.96093884459685, 0.010732322612011683,
               44.04296481225208, 0.0032828357995446046)%>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "s", year = 2014, age = 30)
s_70_2014 <- c(51.88476426439605, 0.002904044089420749,
               56.18163888022552, 0.009343435730013084,
               59.83398394280815, 0.01944444524720057,
               64.13085855863763, 0.03068182026168629,
               68.21288452629288, 0.035479798819043,
               71.75780936260736, 0.03068182026168629,
               76.16210420470499, 0.019318184554850383,
               79.92186949355577, 0.008459601250488509,
               84.00390365684893, 0.00265152270472039) %>% matrix(ncol=2, byrow = T) %>%
  as.data.frame() %>% setNames(c("age_kin","count")) %>%
  mutate(kin = "s", year = 2014, age = 70)

output_time_invariant <- m_30_2014 %>%
  bind_rows(m_30_1947, s_30_1947, s_70_1947, s_30_2014, s_70_2014) %>%
  mutate(age_kin = trunc(age_kin), count=round(count,7))

compare_time_invariant <- kins_japan %>%
  filter(kin %in% c("os", "ys"), alive == "yes") %>%
  group_by(Year, age_focal, age_kin, alive) %>%
  summarise(count = sum(count)) %>%
  mutate(kin = "s") %>%
  bind_rows(kins_japan %>% filter(alive == "yes")) %>%
  select(-year, -cohort, -alive) %>%
  rename(count_demokin = count, year = Year) %>%
  mutate(count_demokin = round(count_demokin,7)) %>%
  right_join(output_time_invariant %>%
               rename(age_focal=age, count_paper = count))

compare_time_invariant %>%
  ggplot() +
  geom_line(aes(age_kin, count_demokin, linetype=factor(year)), col=1)+
  geom_line(aes(age_kin, count_paper, linetype=factor(year)), col=2) +
  facet_grid(~kin+age_focal)+
  theme_bw()


### compare values










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

### plots
# kins alive by age when ego is aged 30 or 70
kins_japan %>%
  filter(age_focal %in% c(30,70), alive=="yes") %>%
  ggplot() +
  geom_line(aes(x=age_kin, y=count,
                color=factor(age_focal), linetype=factor(Year)))  +
  facet_wrap(~kin,scales = "free_y") +
  theme_classic() +
  facet_wrap(~kin,scales = "free_y")
# kins alive during ego´s life
kins_japan %>%
  filter(alive=="yes") %>%
  group_by(Year, kin, age_focal) %>% summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_focal, count, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# experienced deaths
kins_japan %>%
  filter(alive=="no") %>%
  group_by(Year, kin, age_focal) %>% summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_focal, count, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# variation coefficient of age by kin
kins_japan %>%
  filter(alive=="yes") %>%
  group_by(Year, kin, age_focal) %>%
  summarise(mean_age = sum(count*age_kin)/sum(count),
            var_age  = sum(count*age_kin^2)/sum(count) - mean_age^2,
            cv_age  = round(sqrt(var_age)/mean_age*100,1)) %>%
  ggplot() +
  geom_line(aes(age_focal, cv_age, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# dependency ages
kins_japan %>%
  filter(alive=="yes") %>%
  mutate(age_kin_dep = ifelse(age_kin<15,"0-14",
                              ifelse(age_kin<65,"15-64","65+"))) %>%
  group_by(Year, kin, age_focal, age_kin_dep) %>%
  summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_focal, count,
                color = age_kin_dep, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")

















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




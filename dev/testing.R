### expected kin counts function
### based on Caswell (2019)
### functions
# count on ego´s trip   | -
# kins by age today     | -
# kins today, moments   | -
# sum third degree      |
# agreggate degree      |
# death acumulation     |
# plot - diagramm       | -
# prevalence and stages |
# forecast/interpolate  |

library(tidyverse)
library(devtools)
library(HMDHFDplus)
load_all()
SWE_data <- get_HMDHFD(country = "SWE", min_year = 1900, max_year = 2015,
                       pass = rstudioapi::askForPassword(),
                       user = rstudioapi::askForSecret(name="user"))
                       # user = "ivanwilliams1985@gmail.com", pass = "AD2")
debugonce(kins)

#stable
ego_SWE35_stable <- kins(ego_age = 30, year = 1950, P = SWE_data[["P"]],
                         asfr = SWE_data[["asfr"]], stable = TRUE)

ego_SWE35_stable[["kins_by_age_ego"]] %>%
  gather(kin, count, -x) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = 30, color=2)+
  theme_bw() +
  facet_wrap(~kin)


# non-stable
ego_SWE35_nonstable <- kins(ego_age = 30, year = 1990,
                         P = SWE_data[["P"]], asfr = SWE_data[["asfr"]],
                         N = SWE_data[["N"]])

ego_SWE35_nonstable[["kins_by_age_ego"]] %>%
  gather(kin, count, -x) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = 30, color=2)+
  theme_bw() +
  facet_wrap(~kin)

plot_diagramm(ego_SWE35_stable[["kins_total"]],ego_age = 30)


ego_SWE35_nonstable[["kins_total"]]
ego_SWE35_stable[["kins_total"]]









source("Code/R/formal/data.R")
source("Code/R/formal/kins_gral.R")
data <- get_data_HMD(country = "SWE",
                     user = "ivanwilliams1985@gmail.com", pass = "AD2")
U = data[["U"]]
f = data[["f"]]
W = data[["mac"]] %>% as.data.frame()
ego_age_example = 20; year_example = 2015
kin_example <- kins(ego_age = ego_age_example,
                    year = year_example,
                    U = U, f = f, W = W)

kins_by_age_ego <- kin_example[["kins_by_age_ego"]]
kins_by_age_kin <- kin_example[["kins_by_age_kin"]]
kins_total <- kin_example[["kins_total"]]
kins_by_age_ego %>%
  gather(kin, count, -x) %>%
  ggplot() +
  geom_line(aes(x, count))  +
  geom_vline(xintercept = ego_age, color=2)+
  theme_bw() +
  facet_wrap(~kin)
kins_by_age_kin %>%
  select(-x) %>% gather(kin, count, -x_kin) %>%
  ggplot() +
  geom_line(aes(x_kin, count))  +
  geom_vline(xintercept = ego_age, color=2)+
  theme_bw() +
  facet_wrap(~kin)


# load library and data
library(tidyverse)
library(DiagrammeR)


# replicate caswell -------------------------------------------------------

### data: survival probability and fertility by age for Japan
# available at https://www.demographic-research.org/volumes/vol41/24/default.htm
p_1947 = 1 - read.csv("Code/R/formal/caswell/qx_years.csv",
                      header = F, sep = " ")[["V4"]]
f_1947 = read.csv("Code/R/formal/caswell/fx_years.csv",
                  header = F, sep = " ")[["V4"]]
p_2014 = 1 - read.csv("Code/R/formal/caswell/qx_years.csv",
                      header = F, sep = " ")[["V205"]]
f_2014 = read.csv("Code/R/formal/caswell/fx_years.csv",
                  header = F, sep = " ")[["V205"]]


### apply function
kins_japan <- rbind(tibble(Year = 1947, kins_stable(p_1947, f_1947)),
                    tibble(Year = 2014, kins_stable(p_2014, f_2014)))

### diagnostics
# kins alive by age when ego is aged 30 or 70
kins_japan %>%
  filter(age_ego %in% c(30,70), alive=="yes") %>%
  ggplot() +
  geom_line(aes(x=age_kin, y=count,
                color=factor(age_ego), linetype=factor(Year)))  +
  facet_wrap(~kin,scales = "free_y") +
  theme_classic() +
  facet_wrap(~kin,scales = "free_y")
# kins alive duting ego´s life
kins_japan %>%
  filter(alive=="yes") %>%
  group_by(Year, kin, age_ego) %>% summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_ego, count, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# experienced deaths
kins_japan %>%
  filter(alive=="no") %>%
  group_by(Year, kin, age_ego) %>% summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_ego, count, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# variation coefficient of age by kin
kins_japan %>%
  filter(alive=="yes") %>%
  group_by(Year, kin, age_ego) %>%
  summarise(mean_age = sum(count*age)/sum(count),
            var_age  = sum(count*age^2)/sum(count) - mean_age^2,
            cv_age  = round(sqrt(var_age)/mean_age*100,1)) %>%
  ggplot() +
  geom_line(aes(age_ego, cv_age, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")
# dependency ages
kins_japan %>%
  filter(alive=="yes") %>%
  mutate(age_kin_dep = ifelse(age_kin<15,"0-14",
                              ifelse(age_kin<65,"15-64","65+"))) %>%
  group_by(Year, kin, age_ego, age_kin_dep) %>%
  summarise(count = sum(count)) %>%
  ggplot() +
  geom_line(aes(age_ego, count,
                color = age_kin_dep, linetype=factor(Year)))  +
  theme_classic() +
  facet_wrap(~kin, scales = "free_y")

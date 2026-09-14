library(deSolve)
library(tidyverse)
library(zoo)

## Attn: The complementary setup file for NC is commented in detail

#### Phase 1 Target data ####
start <- as.Date("5-1-2020", "%m-%d-%Y") # start of the calibration period

df_CA <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  filter(location != 37) %>%
  mutate(calib_start = difftime(start, min(date), units = "days")) %>%
  mutate(day = difftime(date, min(date), units = "days")) %>% 
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>% 
  mutate(observed = if_else(target == "inc death", value + min_suppressed, value)) %>% 
  group_by(day, target, race_ethnicity, calib_start) %>% 
  summarise(observed = sum(observed, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(target = gsub("inc death", "Deaths", target),
         target = gsub("inc case", "Cases", target)
  )

## check case fatality rate for each group
cfr <- df_CA %>%
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

# CALIBRATION DATA
calibration_CA <- df_CA %>%
  filter(day > as.numeric(gsub("days","",calib_start))) %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 
# 
# test <- calibration_CA %>%
#   ggplot(aes(x = day, y = observed_cum)) + 
#   geom_line() +
#   facet_wrap(race_ethnicity~target, scales = "free")
# test

# Get cumulative reported cases and deaths by day 119
cum_init_CA <- df_CA %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

test_IFR_CA <- cum_init_CA %>% 
  filter(day == 315) %>%
  filter(target == "Deaths") %>%
  mutate(Infections_Estimate_72 = observed_cum/0.0072,
         Infections_Estimate_65 = observed_cum/0.0065
  ) %>%
  mutate(pop = case_when(
    race_ethnicity == "asian" ~ Pop_Asian_CA,
    race_ethnicity == "black" ~ Pop_Black_CA,
    race_ethnicity == "white" ~ Pop_White_CA,
    race_ethnicity == "other" ~ Pop_Other_CA,
    race_ethnicity == "latino" ~ Pop_Latino_CA,
  )) %>%
  group_by(race_ethnicity) %>% 
  reframe(Infs_perc_72 = 100*Infections_Estimate_72/pop,
         Infs_perc_65 = 100*Infections_Estimate_65/pop
  ) %>%
  ungroup()

#### Phase 2 Target data ####
start <- as.Date("5-1-2020", "%m-%d-%Y") # start of the calibration period

df_CA_phase2 <- read_csv("../data/fit/target_data_phase2.csv") %>% 
  filter(location != 37) %>%
  mutate(calib_start = difftime(start, min(date), units = "days")) %>%
  mutate(day = difftime(date, min(date), units = "days")) %>% 
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>% 
  mutate(observed = if_else(target == "inc death", value + min_suppressed, value)) %>% 
  group_by(day, target, race_ethnicity, calib_start) %>% 
  summarise(observed = sum(observed, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(target = gsub("inc death", "Deaths", target),
         target = gsub("inc case", "Cases", target)
  )

## check case fatality rate for first wave
cfr1 <- df_CA_phase2 %>%
  filter(day <= 285) %>%
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

## check case fatality rate for each group - just on 2nd wave
cfr2 <- df_CA_phase2 %>%
  filter(day > 285) %>% 
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

# CALIBRATION DATA
calibration_CA_phase2 <- df_CA_phase2 %>%
  filter(day > as.numeric(gsub("days","",calib_start))) %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

# Get cumulative reported cases and deaths by day 119 - same as in phase 1
cum_init_CA <- df_CA %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

#### Phase 1 Contacts ####

mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
df_contact_ca <- read_csv("../data/mobility_for_spline_ca.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>%
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix for day 311
df_contact_ca_constant <- df_contact_ca %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  mutate(day = 311) %>%
  select(day, alter, ego, value)
## For holidays, allow a 10% increase in contacts starting after day 311 and peaking day 356
contact_max_ca <- df_contact_ca_constant %>%
  mutate(value = 1.1*value) %>%
  mutate(day = 356) %>% 
  select(day, alter, ego, value)
## Then decline back down to end of calibration period
df_contact_ca_constant2 <- df_contact_ca %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(alter, ego, value)
df_contact_ca_constant2 <- expand_grid(day = seq(370, 455, by = 1), df_contact_ca_constant2)

df_contact_ca <- rbind(df_contact_ca, df_contact_ca_constant) %>%
  rbind(contact_max_ca) %>%
  rbind(df_contact_ca_constant2) %>%
  arrange(day)

aca <- contact(df_contact_ca, "asian", "asian")
acb <- contact(df_contact_ca, "asian", "black")
acl <- contact(df_contact_ca, "asian", "latino")
aco <- contact(df_contact_ca, "asian", "other")
acw <- contact(df_contact_ca, "asian", "white")

bca <- contact(df_contact_ca, "black", "asian")
bcb <- contact(df_contact_ca, "black", "black")
bcl <- contact(df_contact_ca, "black", "latino")
bco <- contact(df_contact_ca, "black", "other")
bcw <- contact(df_contact_ca, "black", "white")

lca <- contact(df_contact_ca, "latino", "asian")
lcb <- contact(df_contact_ca, "latino", "black")
lcl <- contact(df_contact_ca, "latino", "latino")
lco <- contact(df_contact_ca, "latino", "other")
lcw <- contact(df_contact_ca, "latino", "white")

oca <- contact(df_contact_ca, "other", "asian")
ocb <- contact(df_contact_ca, "other", "black")
ocl <- contact(df_contact_ca, "other", "latino")
oco <- contact(df_contact_ca, "other", "other")
ocw <- contact(df_contact_ca, "other", "white")

wca <- contact(df_contact_ca, "white", "asian")
wcb <- contact(df_contact_ca, "white", "black")
wcl <- contact(df_contact_ca, "white", "latino")
wco <- contact(df_contact_ca, "white", "other")
wcw <- contact(df_contact_ca, "white", "white")


datnull <- tibble(day = 119:456, cont =  wcw(119:456))

testplot <- datnull %>% 
  ggplot(aes(x = day, y = cont)) +
  geom_line() 
testplot

#### Vaccination ####
vax_CA <- read_csv("../data/vaccination_data.csv") %>% 
  filter(demographic_value %in% c("white", "asian", "other", "latino", "black")) %>% 
  filter(location == "06") %>%
  mutate(day = difftime(date, mindate, units = "days")) %>% 
  mutate(day = gsub(" days", "", day)) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(day, ego = demographic_value, value = full_vax)

## add days before vaccination with zeroes for spline - 2020-12-15 vaccination started, so 2020-12-14

dates <- seq(from = mindate, to = as.Date("2020-12-14", "%Y-%m-%d"), by = 1)

df_add <- expand_grid(
  date = dates, 
  ego = c("white", "asian", "other", "latino", "black"),
  value = c(0)
) %>% 
  mutate(day = difftime(date, mindate, units = "days")) %>% 
  mutate(day = gsub(" days", "", day)) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(day, ego, value)

## get all data 

vax_CA <- rbind(df_add, vax_CA) %>% 
  mutate(pop = case_when(
    ego == "asian" ~ Pop_Asian_CA,
    ego == "black" ~ Pop_Black_CA,
    ego == "latino" ~ Pop_Latino_CA,
    ego == "other" ~ Pop_Other_CA,
    ego == "white" ~ Pop_White_CA
  )) %>%
  mutate(value = value/pop) %>%
  group_by(ego) %>% 
  arrange(day) %>% 
  mutate(delta_value = value - lag(value),
         delta_day = day - lag(day) 
  ) %>% 
  mutate(daily_value = delta_value/delta_day) %>%
  ungroup() %>%
  select(day, ego, daily_value)

vax_a <- vaccinate(vax_CA, "asian")
vax_b <- vaccinate(vax_CA, "black")
vax_l <- vaccinate(vax_CA, "latino")
vax_o <- vaccinate(vax_CA, "other")
vax_w <- vaccinate(vax_CA, "white")


### var parameter ###
date <- c(as.Date("05-01-2020", "%m-%d-%Y"),
          as.Date("06-01-2020", "%m-%d-%Y"),
          as.Date("10-01-2020", "%m-%d-%Y"),
          as.Date("11-15-2020", "%m-%d-%Y"), 
          as.Date("03-15-2021", "%m-%d-%Y"),
          as.Date("04-15-2021", "%m-%d-%Y")
)
var <- c(
  1,
  1,
  1,
  (0.99*1) + (0.01*1.5),
  (0.5*1) + (0.5*1.5),
  (0.44*1) + (0.66*1.5)
)
var_frame <- tibble(date, var) %>%
  mutate(day = difftime(date, mindate, "days")) %>%
  mutate(day = gsub(" days", "", day)) %>%
  mutate(day = as.numeric(day)) %>%
  select(day, var) %>% arrange(day)

variant <- splinefun(x = var_frame$day, y = var_frame$var)


#### Phase 2 contacts - make a step function ####
# Redefine fxns in the phase 2 codes
mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
df_contact_ca2 <- read_csv("../data/mobility_for_spline_ca.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>%
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix from final day of contact
df_contact_ca_constant2a <- df_contact_ca2 %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(day, alter, ego, value)
## For holidays, allow a 10% increase in contacts
contact_max_ca2 <- df_contact_ca_constant2a %>%
  mutate(value = 1.15*value) %>% 
  select(-day)
contact_max_ca2 <- expand_grid(day = seq(306,363, by = 1), contact_max_ca2)
## outside of holidays, make flat contact from final matrix
df_contact_ca_constant2a <- expand_grid(day = c(seq(304,305, by  = 1), seq(364, 455, by = 1)), df_contact_ca_constant2a %>% select(-day))
df_contact_ca2 <- rbind(df_contact_ca2, df_contact_ca_constant2a) %>%
  rbind(contact_max_ca2) %>%
  arrange(day) %>% 
  filter(!is.na(value))


#### Phase 2 counter contacts ####
# Redefine fxns in the phase 2 codes
mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
df_contact_ca2_counter <- read_csv("../data/mobility_for_spline_ca_counterfactual.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>%
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix from final day of contact
df_contact_ca_constant2a_counter <- df_contact_ca2_counter %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(day, alter, ego, value)
## For holidays, allow a 10% increase in contacts
contact_max_ca2_counter <- df_contact_ca_constant2a_counter %>%
  mutate(value = 1.15*value) %>% 
  select(-day)
contact_max_ca2_counter <- expand_grid(day = seq(306,363, by = 1), contact_max_ca2_counter)
## outside of holidays, make flat contact from final matrix
df_contact_ca_constant2a_counter <- expand_grid(day = c(seq(304,305, by  = 1), seq(364, 455, by = 1)), df_contact_ca_constant2a_counter %>% select(-day))
df_contact_ca2_counter <- rbind(df_contact_ca2_counter, df_contact_ca_constant2a_counter) %>%
  rbind(contact_max_ca2_counter) %>%
  arrange(day) %>% 
  filter(!is.na(value))


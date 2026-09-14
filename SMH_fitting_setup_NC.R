library(deSolve)
library(tidyverse)

#### Phase 1 Target data ####
start <- as.Date("5-1-2020", "%m-%d-%Y") # start of the calibration period

df_NC <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  filter(location == 37) %>%
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


## check case fatality rate for each group, this goes in to SMH_params.R
cfr <- df_NC %>%
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

# CALIBRATION DATA - we are fitting on cumulative
calibration_NC <- df_NC %>%
  filter(day > as.numeric(gsub("days","",calib_start))) %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

# grabbing some cumulative counts for initial conditions
cum_init_NC <- df_NC %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

## test IFR vs. reported deaths
test_IFR_NC <- cum_init_NC %>% 
  filter(day == 315) %>%
  filter(target == "Deaths") %>%
  mutate(Infections_Estimate_72 = observed_cum/0.0072,
         Infections_Estimate_65 = observed_cum/0.0065
         ) %>%
  mutate(pop = case_when(
    race_ethnicity == "asian" ~ Pop_Asian_NC,
    race_ethnicity == "black" ~ Pop_Black_NC,
    race_ethnicity == "white" ~ Pop_White_NC,
    race_ethnicity == "other" ~ Pop_Other_NC,
  )) %>%
  group_by(race_ethnicity) %>% 
  reframe(Infs_prop_72 = 100*Infections_Estimate_72/pop,
          Infs_prop_65 = 100*Infections_Estimate_65/pop
          ) %>%
  ungroup()

#### Phase 2 Target data ####
start <- as.Date("5-1-2020", "%m-%d-%Y") # start of the calibration period

df_NC_phase2 <- read_csv("../data/fit/target_data_phase2.csv") %>% 
  filter(location == 37) %>%
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

## check case fatality rate for each group - just on 1st wave
cfr1 <- df_NC_phase2 %>%
  filter(day <= 285) %>%
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

## check case fatality rate for each group - just on 2nd wave
cfr2 <- df_NC_phase2 %>%
  filter(day > 285) %>% 
  pivot_wider(names_from = "target", values_from = "observed") %>% 
  group_by(race_ethnicity) %>% 
  mutate(Deaths = sum(Deaths, na.rm = T), Cases = sum(Cases, na.rm = T)) %>% 
  summarise(CFR = Deaths/Cases) %>%
  ungroup() %>% distinct()

# CALIBRATION DATA
calibration_NC_phase2 <- df_NC_phase2 %>%
  filter(day > as.numeric(gsub("days","",calib_start))) %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

## initial conditions are the same as in phase 1
cum_init_NC <- df_NC %>%
  select(day, target, race_ethnicity, observed) %>%
  group_by(race_ethnicity, target) %>%
  arrange(day) %>%
  mutate(observed_sum = cumsum(observed)) %>%
  ungroup() %>%
  select(day, race_ethnicity, target, observed_cum = observed_sum) 

#### Phase 1 Contacts ####
## Grabbing the minimum date available in the target data, as this is what days are indexed off of
mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
## Here is the mobility data, and we convert dates into days since mindate
df_contact_nc <- read_csv("../data/mobility_for_spline_nc.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>% ## rolling mean of contact data
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix for day 311
df_contact_nc_constant <- df_contact_nc %>% 
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  mutate(day = 311) %>%
  select(day, alter, ego, value)
## For holidays, allow a 20% increase in contacts
contact_max_nc <- df_contact_nc_constant %>%
  mutate(value = 1.1*value) %>%
  mutate(day = 356) %>% 
  select(day, alter, ego, value)
## Then decline back down to end of calibration period
df_contact_nc_constant2 <- df_contact_nc %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(alter, ego, value)
df_contact_nc_constant2 <- expand_grid(day = seq(370, 455, by = 1), df_contact_nc_constant2)

# make a full dataframe with all these contact values
df_contact_nc <- rbind(df_contact_nc, df_contact_nc_constant) %>%
  rbind(contact_max_nc) %>%
  rbind(df_contact_nc_constant2) %>%
  arrange(day)

# construct the splines
aca <- contact(df_contact_nc, "asian", "asian")
acb <- contact(df_contact_nc, "asian", "black")
aco <- contact(df_contact_nc, "asian", "other")
acw <- contact(df_contact_nc, "asian", "white")

bca <- contact(df_contact_nc, "black", "asian")
bcb <- contact(df_contact_nc, "black", "black")
bco <- contact(df_contact_nc, "black", "other")
bcw <- contact(df_contact_nc, "black", "white")

oca <- contact(df_contact_nc, "other", "asian")
ocb <- contact(df_contact_nc, "other", "black")
oco <- contact(df_contact_nc, "other", "other")
ocw <- contact(df_contact_nc, "other", "white")

wca <- contact(df_contact_nc, "white", "asian")
wcb <- contact(df_contact_nc, "white", "black")
wco <- contact(df_contact_nc, "white", "other")
wcw <- contact(df_contact_nc, "white", "white")

# contact visualization
datnull <- tibble(day = 119:456, cont =  wcw(119:456))

testplot <- datnull %>% 
  #filter(day < 400) %>%
  ggplot(aes(x = day, y = cont)) +
  geom_line() 
testplot

#### Vaccination ####
# read the vaccination data
vax_NC <- read_csv("../data/vaccination_data.csv") %>% 
  filter(demographic_value %in% c("white", "asian", "other", "black")) %>% 
  filter(location == "37") %>%
  mutate(day = difftime(date, mindate, units = "days")) %>% 
  mutate(day = gsub(" days", "", day)) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(day, ego = demographic_value, value = full_vax) %>% 
  group_by(day, ego) %>% 
  filter(value == max(value, na.rm = T)) %>% ## ensure only one value reported per day
  ungroup()

## add days before vaccination with zeroes for spline - Jan 20, 2021 vaccination started, so 2021-01-19 but this is 381 days from start and vaccinations are already reported so 
## will end the sequence at 12-13-2020, 1 day before first reported data

dates <- seq(from = mindate, to = as.Date("2020-12-13", "%Y-%m-%d"), by = 1)

df_add <- expand_grid(
  date = dates, 
  ego = c("white", "asian", "other","black"),
  value = c(0)
) %>% 
  mutate(day = difftime(date, mindate, units = "days")) %>% 
  mutate(day = gsub(" days", "", day)) %>% 
  mutate(day = as.numeric(day)) %>% 
  select(day, ego, value)

## get all data 

vax_NC <- rbind(df_add, vax_NC) %>% 
  mutate(pop = case_when(
    ego == "asian" ~ Pop_Asian_NC,
    ego == "black" ~ Pop_Black_NC,
    ego == "other" ~ Pop_Other_NC,
    ego == "white" ~ Pop_White_NC
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

# make splines
vax_a <- vaccinate(vax_NC, "asian")
vax_b <- vaccinate(vax_NC, "black")
vax_o <- vaccinate(vax_NC, "other")
vax_w <- vaccinate(vax_NC, "white")


#### Variant (var parameter) ####

# series of dates with anticipated weighted proportion of new variant
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

# make spline function
variant <- splinefun(x = var_frame$day, y = var_frame$var)

#### Phase 2 contacts - make a step function ####
# Redefine fxns in the phase 2 codes
mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
df_contact_nc2 <- read_csv("../data/mobility_for_spline_nc.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>%
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix from final day of contact
df_contact_nc_constant2a <- df_contact_nc2 %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(day, alter, ego, value)
## For holidays, allow a 10% increase in contacts
contact_max_nc2 <- df_contact_nc_constant2a %>%
  mutate(value = 1.15*value) %>% 
  select(-day)
contact_max_nc2 <- expand_grid(day = seq(306,363, by = 1), contact_max_nc2) # could change back to 315-370
## outside of holidays, make flat contact from final matrix
df_contact_nc_constant2a <- expand_grid(day = c(seq(304,305, by  = 1), seq(364, 455, by = 1)), df_contact_nc_constant2a %>% select(-day))
df_contact_nc2 <- rbind(df_contact_nc2, df_contact_nc_constant2a) %>%
  rbind(contact_max_nc2) %>%
  arrange(day) %>% 
  filter(!is.na(value))


#### Phase 2 counter-contact  ####
# counterfactual contact where disparities are reduced
# Redefine fxns in the phase 2 codes
mindate <- read_csv("../data/fit/target_data_phase1.csv") %>% 
  select(date) %>%
  filter(date == min(date)) %>%
  distinct()
mindate <- mindate$date[1]
df_contact_nc2_counter <- read_csv("../data/mobility_for_spline_nc_counterfactual.csv") %>% 
  mutate(day = difftime(start_date, mindate, units = "days")) %>%
  mutate(day = gsub(" days", "", day), day = as.numeric(day)) %>%
  group_by(alter, ego) %>%
  arrange(day) %>%
  mutate(value = rollmean(value, 3, NA)) %>%
  ungroup() %>%
  select(day, alter, ego, value)
## grab constant matrix from final day of contact
df_contact_nc_constant2a_counter <- df_contact_nc2_counter %>%
  filter(!is.na(value)) %>%
  filter(day == max(day)) %>%
  select(day, alter, ego, value)
## For holidays, allow a 10% increase in contacts
contact_max_nc2_counter <- df_contact_nc_constant2a_counter %>%
  mutate(value = 1.15*value) %>% 
  select(-day)
contact_max_nc2_counter <- expand_grid(day = seq(306,363, by = 1), contact_max_nc2_counter) # could change back to 315-370
## outside of holidays, make flat contact from final matrix
df_contact_nc_constant2a_counter <- expand_grid(day = c(seq(304,305, by  = 1), seq(364, 455, by = 1)), df_contact_nc_constant2a_counter %>% select(-day))
df_contact_nc2_counter <- rbind(df_contact_nc2_counter, df_contact_nc_constant2a_counter) %>%
  rbind(contact_max_nc2_counter) %>%
  arrange(day) %>% 
  filter(!is.na(value))

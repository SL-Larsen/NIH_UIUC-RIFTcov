#### load packages ####
library(tidyverse)
library(deSolve)
library(MetBrewer)
library(directlabels)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(MetBrewer)
library(ggrepel)
library(spatstat.utils)
library(bbmle)
library(Metrics)
library(lhs)
library(scales)
library(Hmisc)
library(zoo)
library(lubridate)
library(parallel)

#### source ####
source("SMH_functions.R")
source("SMH_params.R")
source("SMH_fitting_setup_NC.R")

#### Phase 1 ####
##### initial conditions #####

init_NC <- rep(0,length(names_var_nolatino))
names(init_NC) <- names_var_nolatino
initial_I_NC <- df_NC %>%
  filter(day == 119)
initial_DR_NC <- cum_init_NC %>%
  filter(day == 119)
## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_NC and the cumulative number of reported cases up to 119
init_NC["RepCases_asian"] <- 106
init_NC["D_asian"] <- 3
init_NC["R1_asian"] <- 285

init_NC["RepCases_black"] <- 559
init_NC["D_black"] <- 143
init_NC["R1_black"] <- 3570

init_NC["RepCases_other"] <- 428
init_NC["D_other"] <- 9
init_NC["R1_other"] <- 1346

init_NC["RepCases_white"] <- 1109
init_NC["D_white"] <- 270
init_NC["R1_white"] <- 5445

##### Run maximum likelihood fitting #####
t_start1 = 119
t_end1 = 315

data = calibration_NC
data_calib = calibration_NC
run_parms = NC_parms
run_init = init_NC

start = list(
  #mu = 0.022, 
  rr_asian_scale =  0.9, rr_black_scale =  0.9, rr_other_scale =  0.9, rr_white =  0.15,
  npi_asian_scale = 0.9, npi_black_scale = 0.9, 
  npi_other_scale = 0.9, npi_white = 0.4,
  npi_out = 1/100,
  sig_dist = .187)

upper_limit = list(
  #mu = 0.1, 
  rr_asian_scale = 1, rr_black_scale = 1, 
  rr_other_scale = 1, 
  rr_white = 0.4, 
  npi_asian_scale = 1, npi_black_scale = 1, 
  npi_other_scale = 1, 
  npi_white = 0.9,
  npi_out = 1/3,
  sig_dist = 10
)

lower_limit = list(
  #mu = 0.01, 
  rr_asian_scale = 0.9, rr_black_scale = 0.9, 
  rr_other_scale = 0.9, 
  rr_white = 0.3, 
  npi_asian_scale = 0, npi_black_scale = 0, 
  npi_other_scale = 0, 
  npi_white = 0.1,
  npi_out = 0,
  sig_dist = 0.001
)

m1_NC = mle2(minuslogl = loglik_nc2,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             upper= upper_limit,
             lower = lower_limit,
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_NC)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word")

write_csv(dt_prepan, "../data/out_data/param_estimate_NC.csv")

##### Profile parameters #####
dt_run <-  read_csv("../data/out_data/param_estimate_NC.csv")

param_names <- c(#"mu",
                 "rr_asian_scale", "rr_black_scale", "rr_other_scale", "rr_white",
                 "npi_asian_scale", "npi_black_scale", "npi_other_scale", "npi_white",
                 "npi_out",
                 "sig_dist")
dt_profile <- dt_run %>%
  mutate(prof_param = param_names) %>%
  select(prof_param, Estimate) %>%
  mutate(upper = unlist(unname(upper_limit)),
         lower = unlist(unname(lower_limit))
  )

t_start1 = 119
t_end1 = 315

data = calibration_NC
data_calib = calibration_NC
run_parms = NC_parms
run_init = init_NC


mclapply(1:10, profile_function_NC)

#### Alternate method of fitting the data ####
source("SMH_functions.R")
source("SMH_params.R")
source("SMH_fitting_setup_NC.R")

## Replace contact functions for phase 2 fitting 
aca <- contact(df_contact_nc2, "asian", "asian")
acb <- contact(df_contact_nc2, "asian", "black")
aco <- contact(df_contact_nc2, "asian", "other")
acw <- contact(df_contact_nc2, "asian", "white")

bca <- contact(df_contact_nc2, "black", "asian")
bcb <- contact(df_contact_nc2, "black", "black")
bco <- contact(df_contact_nc2, "black", "other")
bcw <- contact(df_contact_nc2, "black", "white")

oca <- contact(df_contact_nc2, "other", "asian")
ocb <- contact(df_contact_nc2, "other", "black")
oco <- contact(df_contact_nc2, "other", "other")
ocw <- contact(df_contact_nc2, "other", "white")

wca <- contact(df_contact_nc2, "white", "asian")
wcb <- contact(df_contact_nc2, "white", "black")
wco <- contact(df_contact_nc2, "white", "other")
wcw <- contact(df_contact_nc2, "white", "white")

##### run phase 1 mle #####
##### initial conditions #####

init_NC <- rep(0,length(names_var_nolatino))
names(init_NC) <- names_var_nolatino
initial_I_NC <- df_NC %>%
  filter(day == 119)
initial_DR_NC <- cum_init_NC %>%
  filter(day == 119)
## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_NC and the cumulative number of reported cases up to 119
init_NC["RepCases_asian"] <- 106
init_NC["D_asian"] <- 3
init_NC["R1_asian"] <- 285

init_NC["RepCases_black"] <- 559
init_NC["D_black"] <- 143
init_NC["R1_black"] <- 3570

init_NC["RepCases_other"] <- 428
init_NC["D_other"] <- 9
init_NC["R1_other"] <- 1346

init_NC["RepCases_white"] <- 1109
init_NC["D_white"] <- 270
init_NC["R1_white"] <- 5445

##### Run maximum likelihood fitting #####
t_start1 = 119
t_end1 = 285

data = calibration_NC_phase2
data_calib = calibration_NC_phase2
run_parms = NC_parms
run_init = init_NC

start = list(
  rr_asian_scale =  0.95, rr_black_scale =  0.95, 
  rr_other_scale =  0.95, rr_white =  0.95,
  npi_asian_scale = 0.9, npi_black_scale = 0.9, 
  npi_other_scale = 0.9, npi_white = 0.4,
  npi_out = 1/100,
  sig_dist = .187)

upper_limit1 = list(
  rr_asian_scale =  1, rr_black_scale =  1, 
  rr_other_scale =  1, rr_white =  1,
  npi_asian_scale = 1, npi_black_scale = 1, 
  npi_other_scale = 1, 
  npi_white = 0.9,
  npi_out = 1/3,
  sig_dist = 10
)

lower_limit1 = list(
  rr_asian_scale =  0.9, rr_black_scale =  0.9, 
  rr_other_scale =  0.9, rr_white =  0.35,  
  npi_asian_scale = 0.5, npi_black_scale = 0.5, 
  npi_other_scale = 0.5, 
  npi_white = 0.1,
  npi_out = 0,
  sig_dist = 0.001
)

m1_NC = mle2(minuslogl = loglik_nc2_phase1_short,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             upper= upper_limit1,
             lower = lower_limit1,
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_NC)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word") %>%
  mutate(loglik = -m1_NC@min)

write_csv(dt_prepan, "../data/out_data/param_estimate_NC_phase1_for_phase2.csv")

##### Profile Wave 1 #####

dt_run <-  read_csv("../data/out_data/param_estimate_NC_phase1_for_phase2.csv")
out_ll_A <- dt_run$loglik[1]

param_names1 <- c(
  "rr_asian_scale", "rr_black_scale", "rr_other_scale", "rr_white",
  "npi_asian_scale", "npi_black_scale", "npi_other_scale", "npi_white",
  "npi_out",
  "sig_dist"
)

dt_profile <- dt_run %>%
  mutate(prof_param = param_names1) %>%
  select(prof_param, Estimate) %>%
  mutate(upper = unlist(unname(c(upper_limit1))),
         lower = unlist(unname(c(lower_limit1)))
  )

t_start1 = 119
t_end1 = 285

data = calibration_NC
data_calib = calibration_NC
run_parms = NC_parms
run_init = init_NC

mclapply(1:9, profile_function_NC2_ph1)

##### Run first wave, based on params estimated in Phase 1 #####
init_NC <- rep(0,length(names_var_nolatino))
names(init_NC) <- names_var_nolatino
initial_I_NC <- df_NC %>%
  filter(day == 119)
initial_DR_NC <- cum_init_NC %>%
  filter(day == 119)
## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_NC and the cumulative number of reported cases up to 119
init_NC["RepCases_asian"] <- 106
init_NC["D_asian"] <- 3
init_NC["R1_asian"] <- 285

init_NC["RepCases_black"] <- 559
init_NC["D_black"] <- 143
init_NC["R1_black"] <- 3570

init_NC["RepCases_other"] <- 428
init_NC["D_other"] <- 9
init_NC["R1_other"] <- 1346

init_NC["RepCases_white"] <- 1109
init_NC["D_white"] <- 270
init_NC["R1_white"] <- 5445

dt_run <- read_csv("../data/out_data/param_estimate_NC_phase1_for_phase2.csv")

t_start1 = 119
t_end1 = 285

data = calibration_NC
data_calib = calibration_NC
run_parms = NC_parms
run_init = init_NC

paras1 <- run_parms

# reporting rate
paras1["rr_asian"] <- dt_run$Estimate[1] * dt_run$Estimate[4]
paras1["rr_black"] <- dt_run$Estimate[2] * dt_run$Estimate[4]
paras1["rr_other"] <- dt_run$Estimate[3] * dt_run$Estimate[4]
paras1["rr_white"] <- dt_run$Estimate[4]

paras1["npi_out_asian"] <- dt_run$Estimate[9] #npi_out_asian_scale*npi_out_white
paras1["npi_out_black"] <- dt_run$Estimate[9] #npi_out_black_scale*npi_out_white
paras1["npi_out_other"] <- dt_run$Estimate[9] #npi_out_other_scale*npi_out_white
paras1["npi_out_white"] <- dt_run$Estimate[9] #npi_out_white


## set to IFRs over the whole period for each group (from seroprevalence data)
paras1["sigma1_asian"] <- 0.002192918
paras1["sigma1_black"] <- 0.008161487
paras1["sigma1_other"] <- 0.001260678
paras1["sigma1_white"] <- 0.007533983

npi_asian = dt_run$Estimate[5]*dt_run$Estimate[8]
npi_black = dt_run$Estimate[6]*dt_run$Estimate[8]
npi_other = dt_run$Estimate[7]*dt_run$Estimate[8]
npi_white = dt_run$Estimate[8]

paras1["npi_asian"] <- npi_asian
paras1["npi_black"] <- npi_black
paras1["npi_other"] <- npi_other
paras1["npi_white"] <- npi_white

init <- run_init

#To build out initial conditions
# Currently infected: the reported cases for the week, minus anyone who has died this week
# Currently dead: all deaths up to this point
# Currently recovered: all cases up to this point, minus active infections (I), and all deaths up to this point
# Susceptible: everyone who is left

## adjust initial infected based on reporting rate - the number infected is the number reported that day divided by reporting rate
## remove deaths from past week
init["I1_asian"] <- init["RepCases_asian"]/paras1["rr_asian"] - 1
init["I1_black"] <- init["RepCases_black"]/paras1["rr_black"] - 29
init["I1_other"] <- init["RepCases_other"]/paras1["rr_other"] - 2
init["I1_white"] <- init["RepCases_white"]/paras1["rr_white"] - 67

# adjust initial recovered based on reporting rate and deaths ever
init["R1_asian"] <- init["R1_asian"]/paras1["rr_asian"] - init["D_asian"] - init["I1_asian"]
init["R1_black"] <- init["R1_black"]/paras1["rr_black"] - init["D_black"] - init["I1_black"]
init["R1_other"] <- init["R1_other"]/paras1["rr_other"] - init["D_other"] - init["I1_other"]
init["R1_white"] <- init["R1_white"]/paras1["rr_white"] - init["D_white"] - init["I1_white"]

# create susceptibles
init["S1_asian"] <- Pop_Asian_NC - init["I1_asian"] - init["R1_asian"] - init["D_asian"]
init["S1_black"] <- Pop_Black_NC - init["I1_black"] - init["R1_black"] - init["D_black"]
init["S1_other"] <- Pop_Other_NC - init["I1_other"] - init["R1_other"] - init["D_other"]
init["S1_white"] <- Pop_White_NC - init["I1_white"] - init["R1_white"] - init["D_white"]

## set initial stay at home based on npi1 - just S1, I1, R1

init["S1_asian_L"] <- npi_asian*init["S1_asian"]
init["S1_asian"] <- (1-npi_asian)*init["S1_asian"]
init["I1_asian_L"] <- npi_asian*init["I1_asian"]
init["I1_asian"] <- (1-npi_asian)*init["I1_asian"]
init["R1_asian_L"] <- npi_asian*init["R1_asian"]
init["R1_asian"] <- (1-npi_asian)*init["R1_asian"]

init["S1_black_L"] <- npi_black*init["S1_black"]
init["S1_black"] <- (1-npi_black)*init["S1_black"]
init["I1_black_L"] <- npi_black*init["I1_black"]
init["I1_black"] <- (1-npi_black)*init["I1_black"]
init["R1_black_L"] <- npi_black*init["R1_black"]
init["R1_black"] <- (1-npi_black)*init["R1_black"]

init["S1_other_L"] <- npi_other*init["S1_other"]
init["S1_other"] <- (1-npi_other)*init["S1_other"]
init["I1_other_L"] <- npi_other*init["I1_other"]
init["I1_other"] <- (1-npi_other)*init["I1_other"]
init["R1_other_L"] <- npi_other*init["R1_other"]
init["R1_other"] <- (1-npi_other)*init["R1_other"]

init["S1_white_L"] <- npi_white*init["S1_white"]
init["S1_white"] <- (1-npi_white)*init["S1_white"]
init["I1_white_L"] <- npi_white*init["I1_white"]
init["I1_white"] <- (1-npi_white)*init["I1_white"]
init["R1_white_L"] <- npi_white*init["R1_white"]
init["R1_white"] <- (1-npi_white)*init["R1_white"]

times1 <- seq(t_start1,t_end1, by = 1)
#times2 <- seq(t_start2,t_end2, by = 1)

out_calib1 <- ode(y=init, func = seir2_nolatino_vax, times=times1, parms = paras1, method = "euler") %>%
  as.data.frame() %>%
  as_tibble()

##### Run phase 2 MLE #####
t_start2 = 285
t_end2 = 315 #315
t_start3 = 315 #315
t_end3 = 455
  
data = calibration_NC_phase2
data_calib = calibration_NC_phase2

start = list(
  npi_out2 = 0.001065,
  rr_secondwave = 1.1,
  npi_secondwave = 1.3,
  sig_dist = .07749)

upper_limit2 = list(
  npi_out2 = 1/3,
  rr_secondwave = 2,
  npi_secondwave = 4,
  sig_dist = 10
)

lower_limit2 = list(
  npi_out2 = 0,
  rr_secondwave = 1,
  npi_secondwave = 0,
  sig_dist = 0.001
)

m1_NC = mle2(minuslogl = loglik_nc2_phase2_short,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             upper= upper_limit2,
             lower = lower_limit2,
             #fixed = list(npi_halfway = 10),
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_NC)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word") %>%
  mutate(loglik = -m1_NC@min)

write_csv(dt_prepan, "../data/out_data/param_estimate_NC_phase2_alternate.csv")

##### Profile Wave 2 #####

##### Run first wave, based on params estimated in Phase 1 #####
init_NC <- rep(0,length(names_var_nolatino))
names(init_NC) <- names_var_nolatino
initial_I_NC <- df_NC %>%
  filter(day == 119)
initial_DR_NC <- cum_init_NC %>%
  filter(day == 119)
## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_NC and the cumulative number of reported cases up to 119
init_NC["RepCases_asian"] <- 106
init_NC["D_asian"] <- 3
init_NC["R1_asian"] <- 285

init_NC["RepCases_black"] <- 559
init_NC["D_black"] <- 143
init_NC["R1_black"] <- 3570

init_NC["RepCases_other"] <- 428
init_NC["D_other"] <- 9
init_NC["R1_other"] <- 1346

init_NC["RepCases_white"] <- 1109
init_NC["D_white"] <- 270
init_NC["R1_white"] <- 5445

dt_run <- read_csv("../data/out_data/param_estimate_NC_phase1_for_phase2.csv")

t_start1 = 119
t_end1 = 285

data = calibration_NC
data_calib = calibration_NC
run_parms = NC_parms
run_init = init_NC

paras1 <- run_parms
# reporting rate
paras1["rr_asian"] <- dt_run$Estimate[1] * dt_run$Estimate[4]
paras1["rr_black"] <- dt_run$Estimate[2] * dt_run$Estimate[4]
paras1["rr_other"] <- dt_run$Estimate[3] * dt_run$Estimate[4]
paras1["rr_white"] <- dt_run$Estimate[4]

paras1["npi_out_asian"] <- dt_run$Estimate[9] #npi_out_asian_scale*npi_out_white
paras1["npi_out_black"] <- dt_run$Estimate[9] #npi_out_black_scale*npi_out_white
paras1["npi_out_other"] <- dt_run$Estimate[9] #npi_out_other_scale*npi_out_white
paras1["npi_out_white"] <- dt_run$Estimate[9] #npi_out_white


## set to IFRs over the whole period for each group (from seroprevalence data)
paras1["sigma1_asian"] <- 0.002192918
paras1["sigma1_black"] <- 0.008161487
paras1["sigma1_other"] <- 0.001260678
paras1["sigma1_white"] <- 0.007533983

npi_asian = dt_run$Estimate[5]*dt_run$Estimate[8]
npi_black = dt_run$Estimate[6]*dt_run$Estimate[8]
npi_other = dt_run$Estimate[7]*dt_run$Estimate[8]
npi_white = dt_run$Estimate[8]

paras1["npi_asian"] <- npi_asian
paras1["npi_black"] <- npi_black
paras1["npi_other"] <- npi_other
paras1["npi_white"] <- npi_white

init <- run_init

#To build out initial conditions
# Currently infected: the reported cases for the week, minus anyone who has died this week
# Currently dead: all deaths up to this point
# Currently recovered: all cases up to this point, minus active infections (I), and all deaths up to this point
# Susceptible: everyone who is left

## adjust initial infected based on reporting rate - the number infected is the number reported that day divided by reporting rate
## remove deaths from past week
init["I1_asian"] <- init["RepCases_asian"]/paras1["rr_asian"] - 1
init["I1_black"] <- init["RepCases_black"]/paras1["rr_black"] - 29
init["I1_other"] <- init["RepCases_other"]/paras1["rr_other"] - 2
init["I1_white"] <- init["RepCases_white"]/paras1["rr_white"] - 67

# adjust initial recovered based on reporting rate and deaths ever
init["R1_asian"] <- init["R1_asian"]/paras1["rr_asian"] - init["D_asian"] - init["I1_asian"]
init["R1_black"] <- init["R1_black"]/paras1["rr_black"] - init["D_black"] - init["I1_black"]
init["R1_other"] <- init["R1_other"]/paras1["rr_other"] - init["D_other"] - init["I1_other"]
init["R1_white"] <- init["R1_white"]/paras1["rr_white"] - init["D_white"] - init["I1_white"]

# create susceptibles
init["S1_asian"] <- Pop_Asian_NC - init["I1_asian"] - init["R1_asian"] - init["D_asian"]
init["S1_black"] <- Pop_Black_NC - init["I1_black"] - init["R1_black"] - init["D_black"]
init["S1_other"] <- Pop_Other_NC - init["I1_other"] - init["R1_other"] - init["D_other"]
init["S1_white"] <- Pop_White_NC - init["I1_white"] - init["R1_white"] - init["D_white"]

## set initial stay at home based on npi1 - just S1, I1, R1

init["S1_asian_L"] <- npi_asian*init["S1_asian"]
init["S1_asian"] <- (1-npi_asian)*init["S1_asian"]
init["I1_asian_L"] <- npi_asian*init["I1_asian"]
init["I1_asian"] <- (1-npi_asian)*init["I1_asian"]
init["R1_asian_L"] <- npi_asian*init["R1_asian"]
init["R1_asian"] <- (1-npi_asian)*init["R1_asian"]

init["S1_black_L"] <- npi_black*init["S1_black"]
init["S1_black"] <- (1-npi_black)*init["S1_black"]
init["I1_black_L"] <- npi_black*init["I1_black"]
init["I1_black"] <- (1-npi_black)*init["I1_black"]
init["R1_black_L"] <- npi_black*init["R1_black"]
init["R1_black"] <- (1-npi_black)*init["R1_black"]

init["S1_other_L"] <- npi_other*init["S1_other"]
init["S1_other"] <- (1-npi_other)*init["S1_other"]
init["I1_other_L"] <- npi_other*init["I1_other"]
init["I1_other"] <- (1-npi_other)*init["I1_other"]
init["R1_other_L"] <- npi_other*init["R1_other"]
init["R1_other"] <- (1-npi_other)*init["R1_other"]

init["S1_white_L"] <- npi_white*init["S1_white"]
init["S1_white"] <- (1-npi_white)*init["S1_white"]
init["I1_white_L"] <- npi_white*init["I1_white"]
init["I1_white"] <- (1-npi_white)*init["I1_white"]
init["R1_white_L"] <- npi_white*init["R1_white"]
init["R1_white"] <- (1-npi_white)*init["R1_white"]

times1 <- seq(t_start1,t_end1, by = 1)
#times2 <- seq(t_start2,t_end2, by = 1)

out_calib1 <- ode(y=init, func = seir2_nolatino_vax, times=times1, parms = paras1, method = "euler") %>%
  as.data.frame() %>%
  as_tibble()

dt_run2 <-  read_csv("../data/out_data/param_estimate_NC_phase2_alternate.csv")
out_ll_B <- dt_run2$loglik[1]

param_names2 <- c(
  "npi_out2",
  "rr_secondwave",
  "npi_secondwave",
  "sig_dist"
  )

dt_profile <- dt_run2 %>%
  mutate(prof_param = param_names2) %>%
  select(prof_param, Estimate) %>%
  mutate(upper = unlist(unname(c(upper_limit2))),
         lower = unlist(unname(c(lower_limit2)))
  )

t_start2 = 285
t_end2 = 315
t_start3 = 315
t_end3 = 455

data = calibration_NC_phase2
data_calib = calibration_NC_phase2
run_parms = NC_parms
run_init = init_NC

#profile_function_NC2(1)

mclapply(1:3, profile_function_NC2_ph2)

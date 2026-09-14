#### Packages ####
library(deSolve)
library(tidyverse)
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

#### Source ####
source("SMH_functions.R")
source("SMH_params.R")
source("SMH_fitting_setup_CA.R")

#### Phase 1 ####

##### Init #####
init_CA <- rep(0,length(names_var))
names(init_CA) <- names_var

initial_I_CA <- df_CA %>%
  filter(day == 119)
initial_DR_CA <- cum_init_CA %>%
  filter(day == 119)

## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_CA and the cumulative number of reported cases up to 119
init_CA["RepCases_asian"] <- 806
init_CA["D_asian"] <- 404
init_CA["R1_asian"] <- 4244

init_CA["RepCases_black"] <- 450
init_CA["D_black"] <- 265
init_CA["R1_black"] <- 2289

init_CA["RepCases_latino"] <- 5035
init_CA["D_latino"] <- 841
init_CA["R1_latino"] <- 17127

init_CA["RepCases_other"] <- 476
init_CA["D_other"] <- 14
init_CA["R1_other"] <- 2942

init_CA["RepCases_white"] <- 1710
init_CA["D_white"] <- 820
init_CA["R1_white"] <- 9449

##### Run MLE #####

t_start1 = 119
t_end1 = 315

data = calibration_CA
data_calib = calibration_CA
run_parms = CA_parms
run_init = init_CA

start = list(
  #mu = 0.022,
  rr_asian_scale =  0.9, rr_black_scale =  0.9, rr_latino_scale =  0.9, rr_other_scale =  0.9, rr_white =  0.3,
  npi_asian_scale = 0.9, npi_black_scale = 0.9,  npi_latino_scale = 0.9, 
  npi_other_scale = 0.9, npi_white = 0.4,
  npi_out = 1/100,
  # npi_out_asian_scale = 1/30, npi_out_black_scale = 1/30, 
  #          npi_out_other_scale = 1/30, npi_out_white = 1/30,
  sig_dist = .187)

upper_limit = list(
  #mu = 0.1, 
  rr_asian_scale = 1, rr_black_scale = 1, 
  rr_latino_scale = 1,
  rr_other_scale = 1, 
  rr_white = 0.4, 
  npi_asian_scale = 1, npi_black_scale = 1, 
  npi_latino_scale = 1,
  npi_other_scale = 1, 
  npi_white = 0.9,
  npi_out = 1/3,
  sig_dist = 10
)

lower_limit = list(
  #mu = 0.01, 
  rr_asian_scale = 0.9, rr_black_scale = 0.9, 
  rr_latino_scale = 0.9,
  rr_other_scale = 0.9, 
  rr_white = 0.3, 
  npi_asian_scale = 0.1, npi_black_scale = 0.1, 
  npi_latino_scale = 0.1,
  npi_other_scale = 0.1, 
  npi_white = 0.1,
  npi_out = 0,
  sig_dist = 0.001
)

m1_CA = mle2(minuslogl = loglik_ca2,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             #fixed = list(mu = 0.03),
             upper= upper_limit,
             lower = lower_limit,
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_CA)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word")

write_csv(dt_prepan, "../data/out_data/param_estimate_CA.csv")

##### Profile #####

dt_run <- read_csv("../data/out_data/param_estimate_CA.csv")

param_names <- c(#"mu",
                 "rr_asian_scale", "rr_black_scale", "rr_latino_scale", "rr_other_scale", "rr_white",
                 "npi_asian_scale", "npi_black_scale", "npi_latino_scale", "npi_other_scale", "npi_white",
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

data = calibration_CA
data_calib = calibration_CA
run_parms = CA_parms
run_init = init_CA

mclapply(1:12, profile_function_CA)


#### Phase 2 alternate ####

source("SMH_functions.R")
source("SMH_params.R")
source("SMH_fitting_setup_CA.R")

## Reset contacts to step fxn for phase 2
aca <- contact(df_contact_ca2, "asian", "asian")
acb <- contact(df_contact_ca2, "asian", "black")
acl <- contact(df_contact_ca2, "asian", "latino")
aco <- contact(df_contact_ca2, "asian", "other")
acw <- contact(df_contact_ca2, "asian", "white")

bca <- contact(df_contact_ca2, "black", "asian")
bcb <- contact(df_contact_ca2, "black", "black")
bcl <- contact(df_contact_ca2, "black", "latino")
bco <- contact(df_contact_ca2, "black", "other")
bcw <- contact(df_contact_ca2, "black", "white")

lca <- contact(df_contact_ca2, "latino", "asian")
lcb <- contact(df_contact_ca2, "latino", "black")
lcl <- contact(df_contact_ca2, "latino", "latino")
lco <- contact(df_contact_ca2, "latino", "other")
lcw <- contact(df_contact_ca2, "latino", "white")

oca <- contact(df_contact_ca2, "other", "asian")
ocb <- contact(df_contact_ca2, "other", "black")
ocl <- contact(df_contact_ca2, "other", "latino")
oco <- contact(df_contact_ca2, "other", "other")
ocw <- contact(df_contact_ca2, "other", "white")

wca <- contact(df_contact_ca2, "white", "asian")
wcb <- contact(df_contact_ca2, "white", "black")
wcl <- contact(df_contact_ca2, "white", "latino")
wco <- contact(df_contact_ca2, "white", "other")
wcw <- contact(df_contact_ca2, "white", "white")

##### Phase 1 MLE #####
###### Init ######
init_CA <- rep(0,length(names_var))
names(init_CA) <- names_var

initial_I_CA <- df_CA %>%
  filter(day == 119)
initial_DR_CA <- cum_init_CA %>%
  filter(day == 119)

## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_CA and the cumulative number of reported cases up to 119
init_CA["RepCases_asian"] <- 806
init_CA["D_asian"] <- 404
init_CA["R1_asian"] <- 4244

init_CA["RepCases_black"] <- 450
init_CA["D_black"] <- 265
init_CA["R1_black"] <- 2289

init_CA["RepCases_latino"] <- 5035
init_CA["D_latino"] <- 841
init_CA["R1_latino"] <- 17127

init_CA["RepCases_other"] <- 476
init_CA["D_other"] <- 14
init_CA["R1_other"] <- 2942

init_CA["RepCases_white"] <- 1710
init_CA["D_white"] <- 820
init_CA["R1_white"] <- 9449

###### Run MLE ######

t_start1 = 119
t_end1 = 285

data = calibration_CA_phase2
data_calib = calibration_CA_phase2
run_parms = CA_parms
run_init = init_CA

start = list(
  rr_asian_scale =  0.95, rr_black_scale =  0.95, 
  rr_latino_scale =  0.95, rr_other_scale =  0.95, rr_white =  0.95,
  npi_asian_scale = 0.9, npi_black_scale = 0.9,  npi_latino_scale = 0.9, 
  npi_other_scale = 0.9, npi_white = 0.7,
  npi_out = 1/100,
  sig_dist = .187)

upper_limit1 = list(
  rr_asian_scale =  1, rr_black_scale =  1, 
  rr_latino_scale =  1, rr_other_scale =  1, rr_white =  1,
  npi_asian_scale = 1, npi_black_scale = 1, 
  npi_latino_scale = 1,
  npi_other_scale = 1, 
  npi_white = 0.9,
  npi_out = 1/3,
  sig_dist = 10
)

lower_limit1 = list(
  rr_asian_scale =  0.9, rr_black_scale =  0.9, 
  rr_latino_scale =  0.9, rr_other_scale =  0.9, rr_white =  0.35,  
  npi_asian_scale = 0.5, npi_black_scale = 0.5, 
  npi_latino_scale = 0.5,
  npi_other_scale = 0.5, 
  npi_white = 0.1,
  npi_out = 0,
  sig_dist = 0.001
)

m1_CA = mle2(minuslogl = loglik_ca2_phase1_short,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             #fixed = list(mu = 0.03),
             upper= upper_limit1,
             lower = lower_limit1,
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_CA)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word") %>%
  mutate(loglik = -m1_CA@min)

write_csv(dt_prepan, "../data/out_data/param_estimate_CA_phase1_for_phase2.csv")


###### Profile Wave 1 ######

dt_run <-  read_csv("../data/out_data/param_estimate_CA_phase1_for_phase2.csv")
out_ll_A <- dt_run$loglik[1]

param_names1 <- c(
  "rr_asian_scale", "rr_black_scale", "rr_latino_scale", "rr_other_scale", "rr_white",
  "npi_asian_scale", "npi_black_scale", "npi_latino_scale", "npi_other_scale", "npi_white",
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

data = calibration_CA
data_calib = calibration_CA
run_parms = CA_parms
run_init = init_CA

mclapply(1:11, profile_function_CA2_ph1)


##### Phase 2 MLE #####
###### Run first wave, based on params estimated for first wave ######

t_start1 = 119
t_end1 = 285

init_CA <- rep(0,length(names_var))
names(init_CA) <- names_var

initial_I_CA <- df_CA %>%
  filter(day == 119)
initial_DR_CA <- cum_init_CA %>%
  filter(day == 119)

## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_CA and the cumulative number of reported cases up to 119
init_CA["RepCases_asian"] <- 806
init_CA["D_asian"] <- 404
init_CA["R1_asian"] <- 4244

init_CA["RepCases_black"] <- 450
init_CA["D_black"] <- 265
init_CA["R1_black"] <- 2289

init_CA["RepCases_latino"] <- 5035
init_CA["D_latino"] <- 841
init_CA["R1_latino"] <- 17127

init_CA["RepCases_other"] <- 476
init_CA["D_other"] <- 14
init_CA["R1_other"] <- 2942

init_CA["RepCases_white"] <- 1710
init_CA["D_white"] <- 820
init_CA["R1_white"] <- 9449

dt_run <- read_csv("../data/out_data/param_estimate_CA_phase1_for_phase2.csv")

data = calibration_CA_phase2
data_calib = calibration_CA_phase2
run_parms = CA_parms
run_init = init_CA

###### first parameter set ######
paras1 <- run_parms

# reporting rate
paras1["rr_asian"] <- dt_run$Estimate[1] * dt_run$Estimate[5]
paras1["rr_black"] <- dt_run$Estimate[1] * dt_run$Estimate[5]
paras1["rr_latino"] <- dt_run$Estimate[1] * dt_run$Estimate[5]
paras1["rr_other"] <- dt_run$Estimate[1] * dt_run$Estimate[5]
paras1["rr_white"] <- dt_run$Estimate[5]

paras1["npi_out_asian"] <-dt_run$Estimate[11] #dt_run$Estimate[10]*dt_run$Estimate[13]
paras1["npi_out_black"] <-dt_run$Estimate[11] #dt_run$Estimate[11]*dt_run$Estimate[13]
paras1["npi_out_latino"] <-dt_run$Estimate[11] #dt_run$Estimate[11]*dt_run$Estimate[13]
paras1["npi_out_other"] <-dt_run$Estimate[11] #dt_run$Estimate[12]*dt_run$Estimate[13]
paras1["npi_out_white"] <- dt_run$Estimate[11] #dt_run$Estimate[13]

## set to IFRs over the whole period for each group (from seroprevalence data)
paras1["sigma1_asian"] <- 0.013369570
paras1["sigma1_black"] <- 0.014932698
paras1["sigma1_latino"] <- 0.008206442
paras1["sigma1_other"] <- 0.003015897
paras1["sigma1_white"] <- 0.012578309

npi_asian <- dt_run$Estimate[6]*dt_run$Estimate[10]
npi_black <- dt_run$Estimate[7]*dt_run$Estimate[10]
npi_latino <- dt_run$Estimate[8]*dt_run$Estimate[10]
npi_other <- dt_run$Estimate[9]*dt_run$Estimate[10]
npi_white <- dt_run$Estimate[10]

paras1["npi_asian"] <- npi_asian
paras1["npi_black"] <- npi_black
paras1["npi_latino"] <- npi_latino
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
init["I1_asian"] <- init["RepCases_asian"]/paras1["rr_asian"] - 72
init["I1_black"] <- init["RepCases_black"]/paras1["rr_black"] - 54
init["I1_latino"] <- init["RepCases_latino"]/paras1["rr_latino"] - 213
init["I1_other"] <- init["RepCases_other"]/paras1["rr_other"] - 2
init["I1_white"] <- init["RepCases_white"]/paras1["rr_white"] - 174

# adjust initial recovered based on reporting rate and deaths ever
init["R1_asian"] <- init["R1_asian"]/paras1["rr_asian"] - init["D_asian"] - init["I1_asian"]
init["R1_black"] <- init["R1_black"]/paras1["rr_black"] - init["D_black"] - init["I1_black"]
init["R1_latino"] <- init["R1_latino"]/paras1["rr_latino"] - init["D_latino"] - init["I1_latino"]
init["R1_other"] <- init["R1_other"]/paras1["rr_other"] - init["D_other"] - init["I1_other"]
init["R1_white"] <- init["R1_white"]/paras1["rr_white"] - init["D_white"] - init["I1_white"]

# create susceptibles
init["S1_asian"] <- Pop_Asian_CA - init["I1_asian"] - init["R1_asian"] - init["D_asian"]
init["S1_black"] <- Pop_Black_CA - init["I1_black"] - init["R1_black"] - init["D_black"]
init["S1_latino"] <- Pop_Latino_CA - init["I1_latino"] - init["R1_latino"] - init["D_latino"]
init["S1_other"] <- Pop_Other_CA - init["I1_other"] - init["R1_other"] - init["D_other"]
init["S1_white"] <- Pop_White_CA - init["I1_white"] - init["R1_white"] - init["D_white"]

## important to do the lockdown calculation first so that the compartment it's scaled off of does not change before rescaling *that* compartment
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

init["S1_latino_L"] <- npi_latino*init["S1_latino"]
init["S1_latino"] <- (1-npi_latino)*init["S1_latino"]
init["I1_latino_L"] <- npi_latino*init["I1_latino"]
init["I1_latino"] <- (1-npi_latino)*init["I1_latino"]
init["R1_latino_L"] <- npi_latino*init["R1_latino"]
init["R1_latino"] <- (1-npi_latino)*init["R1_latino"]

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

###### Run sim 1 ######
times1 <- seq(t_start1,t_end1, by = 1)

out_calib1 <- ode(y=init, func = seir2_vax, times=times1, parms = paras1, method = "euler") %>%
  as.data.frame() %>%
  as_tibble() 

###### Run MLE ######

t_start2 = 285
t_end2 = 315 #315
t_start3 = 315 #315
t_end3 = 455

data = calibration_CA_phase2
data_calib = calibration_CA_phase2

# 0.8-0.9 optimal for npi_secondwave
start = list(
  npi_out2 = 1/300,
  rr_secondwave = 1.1,
  npi_secondwave = 0.9,
  sig_dist = .130616)

## wave one fits at 73% npi coverage with 0 npi_out, so cap at 25% added npi coverage - meaning can only add 25%
upper_limit2 = list(
  npi_out2 = 1/3,
  rr_secondwave = 2,
  npi_secondwave = 1.4,
  sig_dist = 10
)

lower_limit2 = list(
  npi_out2 = 0,
  rr_secondwave = 1,
  npi_secondwave = 0,
  sig_dist = 0.001
)

m1_CA = mle2(minuslogl = loglik_ca2_phase2_short,
             start = start,
             data = data_calib,
             method = "L-BFGS-B",
             upper= upper_limit2,
             lower = lower_limit2,
             control=list(maxit=20000, trace=TRUE, parscale=abs(unlist(start))))

full <- summary(m1_CA)
dt_prepan <- data.table::as.data.table(coef(full), .keep.rownames = "word") %>%
  mutate(loglik = -m1_CA@min)

write_csv(dt_prepan, "../data/out_data/param_estimate_CA_phase2_alternate.csv")

##### Profile Wave 2 #####

# Run first wave, based on params estimated for first wave

init_CA <- rep(0,length(names_var))
names(init_CA) <- names_var

initial_I_CA <- df_CA %>%
  filter(day == 119)
initial_DR_CA <- cum_init_CA %>%
  filter(day == 119)

## RepCases and D are based on the reported infections/deaths on day 119, since both observed and model data are counted cumulatively from day 119. 
## R is calculated off of calibration_CA and the cumulative number of reported cases up to 119
init_CA["RepCases_asian"] <- 806
init_CA["D_asian"] <- 404
init_CA["R1_asian"] <- 4244

init_CA["RepCases_black"] <- 450
init_CA["D_black"] <- 265
init_CA["R1_black"] <- 2289

init_CA["RepCases_latino"] <- 5035
init_CA["D_latino"] <- 841
init_CA["R1_latino"] <- 17127

init_CA["RepCases_other"] <- 476
init_CA["D_other"] <- 14
init_CA["R1_other"] <- 2942

init_CA["RepCases_white"] <- 1710
init_CA["D_white"] <- 820
init_CA["R1_white"] <- 9449

dt_run <- read_csv("../data/out_data/param_estimate_CA_phase1_for_phase2.csv")

data = calibration_CA_phase2
data_calib = calibration_CA_phase2
run_parms = CA_parms
run_init = init_CA

###### first parameter set ######
paras1 <- run_parms

# reporting rate
paras1["rr_asian"] <- dt_run$Estimate[1] * dt_run$Estimate[5]
paras1["rr_black"] <- dt_run$Estimate[2] * dt_run$Estimate[5]
paras1["rr_latino"] <- dt_run$Estimate[3] * dt_run$Estimate[5]
paras1["rr_other"] <- dt_run$Estimate[4] * dt_run$Estimate[5]
paras1["rr_white"] <- dt_run$Estimate[5]

paras1["npi_out_asian"] <-dt_run$Estimate[11] #dt_run$Estimate[10]*dt_run$Estimate[13]
paras1["npi_out_black"] <-dt_run$Estimate[11] #dt_run$Estimate[11]*dt_run$Estimate[13]
paras1["npi_out_latino"] <-dt_run$Estimate[11] #dt_run$Estimate[11]*dt_run$Estimate[13]
paras1["npi_out_other"] <-dt_run$Estimate[11] #dt_run$Estimate[12]*dt_run$Estimate[13]
paras1["npi_out_white"] <- dt_run$Estimate[11] #dt_run$Estimate[13]

## set to IFRs over the whole period for each group (from seroprevalence data)
paras1["sigma1_asian"] <- 0.013369570
paras1["sigma1_black"] <- 0.014932698
paras1["sigma1_latino"] <- 0.008206442
paras1["sigma1_other"] <- 0.003015897
paras1["sigma1_white"] <- 0.012578309

npi_asian <- dt_run$Estimate[6]*dt_run$Estimate[10]
npi_black <- dt_run$Estimate[7]*dt_run$Estimate[10]
npi_latino <- dt_run$Estimate[8]*dt_run$Estimate[10]
npi_other <- dt_run$Estimate[9]*dt_run$Estimate[10]
npi_white <- dt_run$Estimate[10]

paras1["npi_asian"] <- npi_asian
paras1["npi_black"] <- npi_black
paras1["npi_latino"] <- npi_latino
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
init["I1_asian"] <- init["RepCases_asian"]/paras1["rr_asian"] - 72
init["I1_black"] <- init["RepCases_black"]/paras1["rr_black"] - 54
init["I1_latino"] <- init["RepCases_latino"]/paras1["rr_latino"] - 213
init["I1_other"] <- init["RepCases_other"]/paras1["rr_other"] - 2
init["I1_white"] <- init["RepCases_white"]/paras1["rr_white"] - 174

# adjust initial recovered based on reporting rate and deaths ever
init["R1_asian"] <- init["R1_asian"]/paras1["rr_asian"] - init["D_asian"] - init["I1_asian"]
init["R1_black"] <- init["R1_black"]/paras1["rr_black"] - init["D_black"] - init["I1_black"]
init["R1_latino"] <- init["R1_latino"]/paras1["rr_latino"] - init["D_latino"] - init["I1_latino"]
init["R1_other"] <- init["R1_other"]/paras1["rr_other"] - init["D_other"] - init["I1_other"]
init["R1_white"] <- init["R1_white"]/paras1["rr_white"] - init["D_white"] - init["I1_white"]

# create susceptibles
init["S1_asian"] <- Pop_Asian_CA - init["I1_asian"] - init["R1_asian"] - init["D_asian"]
init["S1_black"] <- Pop_Black_CA - init["I1_black"] - init["R1_black"] - init["D_black"]
init["S1_latino"] <- Pop_Latino_CA - init["I1_latino"] - init["R1_latino"] - init["D_latino"]
init["S1_other"] <- Pop_Other_CA - init["I1_other"] - init["R1_other"] - init["D_other"]
init["S1_white"] <- Pop_White_CA - init["I1_white"] - init["R1_white"] - init["D_white"]

## important to do the lockdown calculation first so that the compartment it's scaled off of does not change before rescaling *that* compartment
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

init["S1_latino_L"] <- npi_latino*init["S1_latino"]
init["S1_latino"] <- (1-npi_latino)*init["S1_latino"]
init["I1_latino_L"] <- npi_latino*init["I1_latino"]
init["I1_latino"] <- (1-npi_latino)*init["I1_latino"]
init["R1_latino_L"] <- npi_latino*init["R1_latino"]
init["R1_latino"] <- (1-npi_latino)*init["R1_latino"]

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

out_calib1 <- ode(y=init, func = seir2_vax, times=times1, parms = paras1, method = "euler") %>%
  as.data.frame() %>%
  as_tibble() 

dt_run2 <-  read_csv("../data/out_data/param_estimate_CA_phase2_alternate.csv")
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
t_end2 = 315 #315
t_start3 = 315 #315
t_end3 = 455

data = calibration_CA_phase2
data_calib = calibration_CA_phase2
run_parms = CA_parms
run_init = init_CA

mclapply(1:3, profile_function_CA2_ph2)
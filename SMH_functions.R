library(deSolve)
library(tidyverse)

#### Contact/npi splines ####
# This function takes in a set of contact data, an ego (race/ethnicity) and an alter (race/ethnicity)
# It creates a spline function using the contact values for each day
contact <- function(data_contact, ego_string, alter_string){
  dfc <- data_contact %>%
    filter(ego == ego_string, alter == alter_string) %>% 
    arrange(day)
  contF <- splinefun(x = dfc$day, y = dfc$value)
  return(contF)
}

#### Vaccination ####
# This function takes in a set of vaccination data, and an ego (race/ethnicity). 
# It creates a spline function using the daily vaccination rate for each day (NOT cumulative)
vaccinate <- function(data_vax, ego_string){
  dfv <- data_vax %>% 
    filter(ego == ego_string) %>%
    arrange(day)
  vaxF <- splinefun(x = dfv$day, y = dfv$daily_value)
  return(vaxF)
}


#### NPI ####
type_III_resp <- function (t, Smax, WVH, k) {
  
  return(Smax*(t^k)/((WVH^k) + (t^k)))
}

npi_weekly <- function (Smax, WVH, k){
  datOut <- tibble(Week = 1:60, Npi = type_III_resp(1:60, Smax,WVH, k)) %>%
    mutate(dailyNpi =  c(Npi[1],diff(Npi)),
           NpiCumSum = cumsum(dailyNpi)) %>% 
    arrange(Week) %>%
    mutate(Day = Week*7-7+1, dailyNpi=dailyNpi/7)
  
  npiF <- splinefun(x = datOut$Day, y = datOut$dailyNpi)
  return(npiF)
}


#### Phase 1 Log likelihood ####

# These are the log likelihood functions for fitting Phase 1 data. There are separate functions for NC and CA. 
# The function for NC is commented with explanations. 

loglik_nc2 <- function(
    rr_asian_scale, rr_black_scale, rr_other_scale, rr_white, 
    npi_asian_scale, npi_black_scale, npi_other_scale, npi_white,
    npi_out,
    sig_dist
){
  
  # First, we feed in the run_parms from SMH_params.R
  paras1 <- run_parms

  # Next, any params that are going to get fit need to be adjusted. We pull these from the arguments to the log likelihood function
  paras1["rr_asian"] <- rr_asian_scale*rr_white
  paras1["rr_black"] <- rr_black_scale*rr_white
  paras1["rr_other"] <- rr_other_scale*rr_white
  paras1["rr_white"] <- rr_white
  paras1["npi_out_asian"] <- npi_out
  paras1["npi_out_black"] <- npi_out
  paras1["npi_out_other"] <- npi_out
  paras1["npi_out_white"] <- npi_out
  
  # Using the reporting rate that is being tested, combined with the case fatality rate for phase 1 which we know,
  # we can infer the infection fatality rate for each group
  paras1["sigma1_asian"] <- paras1["sigma1_asian"]*paras1["rr_asian"]
  paras1["sigma1_black"] <- paras1["sigma1_black"]*paras1["rr_black"]
  paras1["sigma1_other"] <- paras1["sigma1_other"]*paras1["rr_other"]
  paras1["sigma1_white"] <- paras1["sigma1_white"]*paras1["rr_white"]
  
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
  
  ## set initial stay at home based on npi's - S1, I1, R1
  
  npi_asian = npi_asian_scale*npi_white
  npi_black = npi_black_scale*npi_white
  npi_other = npi_other_scale*npi_white
  
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
  
  # Set times
  times1 <- seq(t_start1,t_end1, by = 1)
  
  # Run the model 
  out_calib1 <- ode(y=init, func = seir2_nolatino_vax, times=times1, parms = paras1, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  # Reformat the model output, and inner join with the observed data 
  out_calib1 <- out_calib1 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 

  # Normalize the data. Here, I'm using the maximum in the observed data for each target and 
  # race/ethnicity combination. This step is very important. Without it, the MLE weights targets and race/ethnicity differently,
  # and doesn't converge. 
  out_calib <- out_calib1 %>% 
    group_by(target, race_ethnicity) %>%
    mutate(max_normal = max(observed_cum, na.rm = T),
           observed_cum = observed_cum/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    ungroup()
  
  # Now get the negative log likelihood. Note we are concurrently fitting the sd
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed_cum,sd=sig_dist,log=TRUE))
  
  # This code ensures that if the output of the log likelihood is NA, the MLE won't stop. Instead it will return
  # an extremely large, positive value of the negative log likelihood and keep iterating.
  ll_final = if_else(is.na(ll), 10^6, ll)
  
  return(ll_final)
}


loglik_ca2 <- function(
    rr_asian_scale, rr_black_scale, rr_latino_scale, rr_other_scale, rr_white, 
    npi_asian_scale, npi_black_scale, npi_latino_scale, npi_other_scale, npi_white,
    npi_out,
    sig_dist
){
  
  paras1 <- run_parms
  paras1["rr_asian"] <- rr_asian_scale*rr_white
  paras1["rr_black"] <- rr_black_scale*rr_white
  paras1["rr_latino"] <- rr_latino_scale*rr_white
  paras1["rr_other"] <- rr_other_scale*rr_white
  paras1["rr_white"] <- rr_white
  paras1["npi_out_asian"] <- npi_out
  paras1["npi_out_black"] <- npi_out
  paras1["npi_out_latino"] <- npi_out
  paras1["npi_out_other"] <- npi_out
  paras1["npi_out_white"] <- npi_out
  paras1["sigma1_asian"] <- paras1["sigma1_asian"]*paras1["rr_asian"]
  paras1["sigma1_black"] <- paras1["sigma1_black"]*paras1["rr_black"]
  paras1["sigma1_latino"] <- paras1["sigma1_latino"]*paras1["rr_latino"]
  paras1["sigma1_other"] <- paras1["sigma1_other"]*paras1["rr_other"]
  paras1["sigma1_white"] <- paras1["sigma1_white"]*paras1["rr_white"]

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
  
  
  ## set initial stay at home based on npi1 - just S1, I1, R1
  
  npi_asian = npi_asian_scale*npi_white
  npi_black = npi_black_scale*npi_white
  npi_latino = npi_latino_scale*npi_white
  npi_other = npi_other_scale*npi_white
  
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
  out_calib1 <- out_calib1 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 

  out_calib <- out_calib1 %>% 
    group_by(target, race_ethnicity) %>%
    mutate(max_normal = max(observed_cum, na.rm = T),
           observed_cum = observed_cum/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    ungroup()
  
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed_cum,sd=sig_dist,log=TRUE))
  
  ll_final = if_else(is.na(ll), 10^6, ll)
  
  return(ll_final)
}

#### Phase 1 Profiling functions ####

# These are separate profiling functions for NC and CA. The function works by taking in an index i
# corresponding to a parameter that was fit during the MLE process. 
# We take the initial estimate for this parameter and the others, and run the MLE to get a maximum likelihood. 
# If params are on the boundary, we need to add
# or subtract a super small value because the MLE gets mad if we directly input boundary values. 
# Once we get the likelihood out, we subtract 1.92 to get the confidence boundary.
# Then, we do a right (+) and left (-) profile by iterating out from the original parameter estimate:
# We add some starting amount to the parameter estimate (1, 0.1, etc)
# Anytime a +value exceeds the parameter boundary, we divide it by 10 without running the MLE
# and try again until reaching a predefined threshold of accuracy
# If the +value does not push the parameter across the boundary, we run the MLE, checking if the LL is within
# or outside of the confidence interval. If it's outside, we go back and divide by 10, and try re-adding this smaller value. 
# If it's inside, we add the small value again. (e.g. param + 0.1 --> param + 0.1 + 0.1)
# Using this procedure, coupled with our threshold value, we can approximate the confidence interval.

profile_function_NC<- function(i){
  
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms. 
  start = list(
    #mu = dt_run$Estimate[1],
    rr_asian_scale =  dt_run$Estimate[1]+0.001, rr_black_scale =  dt_run$Estimate[2],
    rr_other_scale =  dt_run$Estimate[3]-0.001, rr_white =  dt_run$Estimate[4]+0.001,
    npi_asian_scale = dt_run$Estimate[5]-0.01, npi_black_scale = dt_run$Estimate[6]+0.001,  
    npi_other_scale = dt_run$Estimate[7]+0.001, npi_white = dt_run$Estimate[8]+0.001,
    npi_out = dt_run$Estimate[9]+0.001,
    sig_dist = dt_run$Estimate[10])
  upper_list = upper_limit
  lower_list = lower_limit
  fixed_list = list(#mu = parm_test,
                    rr_asian_scale= parm_test, rr_black_scale= parm_test, 
                    rr_other_scale = parm_test, rr_white = parm_test,
                    npi_asian_scale = parm_test, npi_black_scale = parm_test, 
                    npi_other_scale = parm_test, npi_white = parm_test,
                    npi_out = parm_test,
                    sig_dist = parm_test
  )
  start2 <- start[-i] ## remove the one we are fixing
  upper_list2 <- upper_list[-i]
  lower_list2 <- lower_list[-i]
  fixed_list2 <- fixed_list[i] ## make a list containing only the param we are fixing
  
  print("Running start MLE")
  ## run m1_A for starting param
  m1_NC = mle2(minuslogl = loglik_nc2,
               start = start2,
               data = data_calib,
               fixed = fixed_list2,
               method="L-BFGS-B",
               upper= upper_list2,
               lower = lower_list2,
               control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
  )
  
  out_ll <- -m1_NC@min ## grab the base likelihood (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
  conf_threshold <- out_ll - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] ## start at the true parameter value
  upper_limit <- dt_profile$upper[i] ## This is bad coding practice but works because R is an unintuitive language with weird environment rules. I will fix this at some point to not look like a "rewrite"
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    param_test_right = 10^3
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_right >= upper_limit){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_list2 = list(#mu = parm_test_right,
                           rr_asian_scale = parm_test_right, rr_black_scale= parm_test_right,
                           rr_other_scale = parm_test_right, rr_white = parm_test_right,
                           npi_asian_scale = parm_test_right, npi_black_scale = parm_test_right, 
                           npi_other_scale = parm_test_right, npi_white = parm_test_right,
                           npi_out = parm_test_right,
                           sig_dist = parm_test_right
        )
        fixed_list2 <- fixed_list2[i]
        ##
        
        ## run mle and grab the new log likelihood
        m1_NC = mle2(minuslogl = loglik_nc2,
                     start = start2,
                     data = data_calib,
                     fixed = fixed_list2,
                     method="L-BFGS-B",
                     upper= upper_list2,
                     lower = lower_list2,
                     control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
        )
        new_ll_right <- -m1_NC@min
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  ## get left threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i] ## start at the true parameter value
  lower_limit <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    param_test_left = 10^3
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_list2 <- list(#mu = parm_test_left,
                            rr_asian_scale= parm_test_left, rr_black_scale= parm_test_left, 
                            rr_other_scale = parm_test_left, rr_white = parm_test_left,
                            npi_asian_scale = parm_test_left, npi_black_scale = parm_test_left, 
                            npi_other_scale = parm_test_left, npi_white = parm_test_left,
                            npi_out = parm_test_left, 
                            sig_dist = parm_test_left
        )
        fixed_list2 <- fixed_list2[i]
        ##
        
        ## run mle and grab the new log likelihood
        ## run mle and grab the new log likelihood
        m1_NC = mle2(minuslogl = loglik_nc2,
                     start = start2,
                     data = data_calib,
                     fixed = fixed_list2,
                     method="L-BFGS-B",
                     upper= upper_list2,
                     lower = lower_list2,
                     control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
        )
        new_ll_left <- -m1_NC@min
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_out_data_NC/", dt_profile$prof_param[i], ".csv"))
  
}

profile_function_CA<- function(i){
  
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the mle
  start = list(
    #mu = dt_run$Estimate[1],
    rr_asian_scale =  dt_run$Estimate[1]+0.001, rr_black_scale =  dt_run$Estimate[2]-0.001, rr_latino_scale =  dt_run$Estimate[3]+0.001,
    rr_other_scale =  dt_run$Estimate[4]-0.001, rr_white = dt_run$Estimate[5]+0.001,
    npi_asian_scale = dt_run$Estimate[6]-0.001, npi_black_scale = dt_run$Estimate[7],  npi_latino_scale = dt_run$Estimate[8], 
    npi_other_scale = dt_run$Estimate[9]+0.001, npi_white = dt_run$Estimate[10]+0.001,
    npi_out = dt_run$Estimate[11]+0.001,
    # npi_out_asian_scale = 1/30, npi_out_black_scale = 1/30, 
    #          npi_out_other_scale = 1/30, npi_out_white = 1/30,
    sig_dist = dt_run$Estimate[12])
  upper_list = upper_limit
  lower_list = lower_limit
  fixed_list = list(#mu = parm_test,
                    rr_asian_scale= parm_test, rr_black_scale= parm_test, rr_latino_scale = parm_test,
                    rr_other_scale = parm_test, rr_white = parm_test,
                    npi_asian_scale = parm_test, npi_black_scale = parm_test, npi_latino_scale = parm_test,
                    npi_other_scale = parm_test, npi_white = parm_test,
                    npi_out = parm_test,
                    sig_dist = parm_test
  )
  start2 <- start[-i] ## remove the one we are fixing
  upper_list2 <- upper_list[-i]
  lower_list2 <- lower_list[-i]
  fixed_list2 <- fixed_list[i] ## make a list containing only the param we are fixing
  
  print("Running start MLE")
  ## run m1_A for starting param
  m1_CA = mle2(minuslogl = loglik_ca2,
               start = start2,
               data = data_calib,
               fixed = fixed_list2,
               method="L-BFGS-B",
               upper= upper_list2,
               lower = lower_list2,
               control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
  )
  
  out_ll <- -m1_CA@min ## grab the base likelihood (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
  conf_threshold <- out_ll - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] ## start at the true parameter value
  upper_limit <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    param_test_right = 10^3
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_right >= upper_limit){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_list2 = list(#mu = parm_test_right,
                           rr_asian_scale = parm_test_right, rr_black_scale= parm_test_right, rr_latino_scale = parm_test_right,
                           rr_other_scale = parm_test_right, rr_white = parm_test_right,
                           npi_asian_scale = parm_test_right, npi_black_scale = parm_test_right, npi_latino_scale = parm_test_right,
                           npi_other_scale = parm_test_right, npi_white = parm_test_right,
                           npi_out = parm_test_right,
                           sig_dist = parm_test_right
        )
        fixed_list2 <- fixed_list2[i]
        ##
        
        ## run mle and grab the new log likelihood
        m1_CA = mle2(minuslogl = loglik_ca2,
                     start = start2,
                     data = data_calib,
                     fixed = fixed_list2,
                     method="L-BFGS-B",
                     upper= upper_list2,
                     lower = lower_list2,
                     control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
        )
        new_ll_right <- -m1_CA@min
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  ## get left threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i] ## start at the true parameter value
  lower_limit <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    param_test_left = 10^3
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_list2 <- list(#mu = parm_test_left,
                            rr_asian_scale= parm_test_left, rr_black_scale= parm_test_left, rr_latino_scale = parm_test_left,
                            rr_other_scale = parm_test_left, rr_white = parm_test_left,
                            npi_asian_scale = parm_test_left, npi_black_scale = parm_test_left, npi_latino_scale = parm_test_left,
                            npi_other_scale = parm_test_left, npi_white = parm_test_left,
                            npi_out = parm_test_left, 
                            sig_dist = parm_test_left
        )
        fixed_list2 <- fixed_list2[i]
        ##
        
        ## run mle and grab the new log likelihood
        ## run mle and grab the new log likelihood
        m1_CA = mle2(minuslogl = loglik_ca2,
                     start = start2,
                     data = data_calib,
                     fixed = fixed_list2,
                     method="L-BFGS-B",
                     upper= upper_list2,
                     lower = lower_list2,
                     control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(start2)))
        )
        new_ll_left <- -m1_CA@min
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_out_data_CA/", dt_profile$prof_param[i], ".csv"))
  
}

#### Phase 2 log likelihood alternate ####

# This is a second pass at the log likelihood functions for phase 2. In this method, we use an adapted
# phase 1 function to fit phase 1, then use this function to fit any scaling parameters on those
# params moving into phase 2

# These are the log likelihood functions for fitting Phase 1 data. There are separate functions for NC and CA. 
# The function for NC is commented with explanations. 

loglik_nc2_phase1_short <- function(
    rr_asian_scale, rr_black_scale, rr_other_scale, rr_white, 
    npi_asian_scale, npi_black_scale, npi_other_scale, npi_white,
    npi_out,
    sig_dist
){
  
  # First, we feed in the run_parms from SMH_params.R
  paras1 <- run_parms
  
  # reporting rate
  paras1["rr_asian"] <- rr_asian_scale * rr_white
  paras1["rr_black"] <- rr_black_scale * rr_white
  paras1["rr_other"] <- rr_other_scale * rr_white
  paras1["rr_white"] <- rr_white
  
  # Next, any params that are going to get fit need to be adjusted. We pull these from the arguments to the log likelihood function
  paras1["npi_out_asian"] <- npi_out
  paras1["npi_out_black"] <- npi_out
  paras1["npi_out_other"] <- npi_out
  paras1["npi_out_white"] <- npi_out
  
  ## set to IFRs over the whole period for each group (from seroprevalence data)
  paras1["sigma1_asian"] <- 0.002192918
  paras1["sigma1_black"] <- 0.008161487
  paras1["sigma1_other"] <- 0.001260678
  paras1["sigma1_white"] <- 0.007533983
  
  npi_asian = npi_asian_scale*npi_white
  npi_black = npi_black_scale*npi_white
  npi_other = npi_other_scale*npi_white
  
  ## NPIs
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
  
  ## set initial stay at home based on npi's - S1, I1, R1
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
  
  # Set times
  times1 <- seq(t_start1,t_end1, by = 1)
  
  # Run the model 
  out_calib1 <- ode(y=init, func = seir2_nolatino_vax, times=times1, parms = paras1, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  # Reformat the model output, and inner join with the observed data 
  out_calib1 <- out_calib1 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 
  
  # Normalize the data. Here, I'm using the maximum in the observed data for each target and 
  # race/ethnicity combination. This step is very important. Without it, the MLE weights targets and race/ethnicity differently,
  # and doesn't converge. 
  out_calib <- out_calib1 %>% 
    group_by(target, race_ethnicity) %>%
    arrange(day) %>%
    mutate(observed = observed_cum - lag(observed_cum),
           model_predicted = model_predicted - lag(model_predicted)
           ) %>%
    mutate(max_normal = max(observed, na.rm = T),
           observed = observed/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    filter(!is.na(observed), !is.na(model_predicted)) %>%
    ungroup()
  
  # Now get the negative log likelihood. Note we are concurrently fitting the sd
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed,sd=sig_dist,log=TRUE))
  
  # This code ensures that if the output of the log likelihood is NA, the MLE won't stop. Instead it will return
  # an extremely large, positive value of the negative log likelihood and keep iterating.
  ll_final = if_else(is.na(ll), 10^6, ll)
  
  return(ll_final)
}


loglik_ca2_phase1_short <- function(
    rr_asian_scale, rr_black_scale, rr_latino_scale, rr_other_scale, rr_white, 
    npi_asian_scale, npi_black_scale, npi_latino_scale, npi_other_scale, npi_white,
    npi_out,
    sig_dist
){
  
  paras1 <- run_parms
  
  # reporting rate
  paras1["rr_asian"] <- rr_asian_scale * rr_white
  paras1["rr_black"] <- rr_black_scale * rr_white
  paras1["rr_latino"] <- rr_latino_scale * rr_white
  paras1["rr_other"] <- rr_other_scale * rr_white
  paras1["rr_white"] <- rr_white
  
  paras1["npi_out_asian"] <- npi_out
  paras1["npi_out_black"] <- npi_out
  paras1["npi_out_latino"] <- npi_out
  paras1["npi_out_other"] <- npi_out
  paras1["npi_out_white"] <- npi_out
  
  ## set to IFRs over the whole period for each group (from seroprevalence data)
  paras1["sigma1_asian"] <- 0.013369570
  paras1["sigma1_black"] <- 0.014932698
  paras1["sigma1_latino"] <- 0.008206442
  paras1["sigma1_other"] <- 0.003015897
  paras1["sigma1_white"] <- 0.012578309
  
  npi_asian = npi_asian_scale*npi_white
  npi_black = npi_black_scale*npi_white
  npi_latino = npi_latino_scale*npi_white
  npi_other = npi_other_scale*npi_white
  
  ## NPIs
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
  
  
  ## set initial stay at home based on npi1 - just S1, I1, R1
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
  out_calib1 <- out_calib1 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 
  
  out_calib <- out_calib1 %>% 
    group_by(target, race_ethnicity) %>%
    arrange(day) %>%
    mutate(observed = observed_cum - lag(observed_cum),
           model_predicted = model_predicted - lag(model_predicted)
    ) %>%
    mutate(max_normal = max(observed, na.rm = T),
           observed = observed/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    filter(!is.na(observed), !is.na(model_predicted)) %>%
    ungroup()
  
    # group_by(target, race_ethnicity) %>%
    # mutate(max_normal = max(observed_cum, na.rm = T),
    #        observed_cum = observed_cum/max_normal,
    #        model_predicted = model_predicted/max_normal) %>%
    # arrange(day) %>%
    # mutate(observed = observed_cum - lag(observed_cum),
    #        model_predicted = model_predicted - lag(model_predicted)
    # ) %>%
    # filter(!is.na(observed), !is.na(model_predicted)) %>%
    # ungroup()
  
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed,sd=sig_dist,log=TRUE))
  
  ll_final = if_else(is.na(ll), 10^6, ll)
  
  return(ll_final)
}

loglik_nc2_phase2_short <- function(
    npi_out2,
    rr_secondwave,
    npi_secondwave,
    #npi_halfway, 
    sig_dist
){

  ## first parm set - here we're going to grab the params from phase 1, which will have been set up 
  ## while running the phase 1 simulation 
  paras2 <- paras1
  ## The rate of leaving npi's from phase 1 is changed for phase 2
  paras2["npi_out_asian"] <- npi_out2
  paras2["npi_out_black"] <- npi_out2
  paras2["npi_out_other"] <- npi_out2 
  paras2["npi_out_white"] <- npi_out2 
  ## Reporting rate change
  paras2["rr_asian"] <- rr_secondwave*paras2["rr_asian"]
  paras2["rr_black"] <- rr_secondwave*paras2["rr_black"]
  paras2["rr_other"] <- rr_secondwave*paras2["rr_other"]
  paras2["rr_white"] <- rr_secondwave*paras2["rr_white"]
  
  paras2["npi_secondwave"] <- npi_secondwave
  #paras2["npi_halfway"] <- npi_halfway
  
  ## set time for simulations
  times2 <- seq(t_start2,t_end2, by = 1)
  times3 <- seq(t_start3, t_end3, by = 1)
  
  # init for second run 
  init2 <- out_calib1 %>% 
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist() 

  # second run
  out_calib2 <- ode(y=init2, func = seir2_nolatino_vax, times=times2, parms = paras2, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  
  # init for third run 
  init3 <- out_calib2 %>%
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist()
  
  # format 2nd run

  out_calib2 <- out_calib2 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) %>%
    filter(day != t_end2) # remove end day - this is init for third run

  ## third run - no param changes for this run, just init changes
  out_calib3 <- ode(y=init3, func = seir2_nolatino_vax2, times=times3, parms = paras2, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  
  out_calib3 <- out_calib3 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 
  
  ## grab the full output and normalize it for comparison, the same as phase 1
  out_calib <- rbind(out_calib2, out_calib3) %>%
    group_by(target, race_ethnicity) %>%
    arrange(day) %>%
    mutate(observed = observed_cum - lag(observed_cum),
           model_predicted = model_predicted - lag(model_predicted)
    ) %>%
    mutate(max_normal = max(observed, na.rm = T),
           observed = observed/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    filter(!is.na(observed), !is.na(model_predicted)) %>%
    ungroup()
  
    # group_by(target, race_ethnicity) %>%
    # mutate(max_normal = max(observed_cum, na.rm = T),
    #        observed_cum = observed_cum/max_normal,
    #        model_predicted = model_predicted/max_normal) %>%
    # arrange(day) %>%
    # mutate(observed = observed_cum - lag(observed_cum),
    #        model_predicted = model_predicted - lag(model_predicted)
    # ) %>%
    # filter(!is.na(observed), !is.na(model_predicted)) %>%
    # ungroup()
  
  ## calculate the negative loglikelihood
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed,sd=sig_dist,log=TRUE))
  
  ## NA insurance policy so that simulations don't break
  ll_final = if_else(is.na(ll) | is.infinite(ll), 10^6, ll)
  
  return(ll_final)
}

loglik_ca2_phase2_short <- function(
    npi_out2,
    rr_secondwave,
    npi_secondwave,
    #npi_halfway,
    sig_dist
){
  
  #print("running")
  paras2 <- paras1
  paras2["npi_out_asian"] <- npi_out2
  paras2["npi_out_black"] <- npi_out2
  paras2["npi_out_latino"] <- npi_out2
  paras2["npi_out_other"] <- npi_out2
  paras2["npi_out_white"] <- npi_out2
  ## Reporting rate change
  paras2["rr_asian"] <- rr_secondwave*paras2["rr_asian"]
  paras2["rr_black"] <- rr_secondwave*paras2["rr_black"]
  paras2["rr_latino"] <- rr_secondwave*paras2["rr_latino"]
  paras2["rr_other"] <- rr_secondwave*paras2["rr_other"]
  paras2["rr_white"] <- rr_secondwave*paras2["rr_white"]
  
  paras2["npi_secondwave"] <- npi_secondwave
  #paras2["npi_halfway"] <- npi_halfway

  init <- run_init
  
  times2 <- seq(t_start2,t_end2, by = 1)
  times3 <- seq(t_start3,t_end3, by = 1)
  
  # init for second run 
  init2 <- out_calib1 %>% 
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist() 
  
  # second run 
  out_calib2 <- ode(y=init2, func = seir2_vax, times=times2, parms = paras2, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  
  # init for third run 
  init3 <- out_calib2 %>%
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist()
  
  # format 2nd run 
  out_calib2 <- out_calib2 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) %>%
    filter(day != t_end2) # filter final day - init for run 3
  
  # third run 
  out_calib3 <- ode(y=init3, func = seir2_vax2, times=times3, parms = paras2, method = "euler") %>% # no parm change
    as.data.frame() %>%
    as_tibble()
  
  out_calib3 <- out_calib3 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 
  
  out_calib <- rbind(out_calib2, out_calib3) %>%
    group_by(target, race_ethnicity) %>%
    arrange(day) %>%
    mutate(observed = observed_cum - lag(observed_cum),
           model_predicted = model_predicted - lag(model_predicted)
    ) %>%
    mutate(max_normal = max(observed, na.rm = T),
           observed = observed/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    filter(!is.na(observed), !is.na(model_predicted)) %>%
    ungroup()
    # # mutate(phase = if_else(day < t_start2, "Phase 1", "Phase 2")) %>%
    # # group_by(target, race_ethnicity, phase) %>%
    # group_by(target, race_ethnicity) %>%
    # mutate(max_normal = max(observed_cum, na.rm = T),
    #        observed_cum = observed_cum/max_normal,
    #        model_predicted = model_predicted/max_normal) %>%
    # arrange(day) %>%
    # mutate(observed = observed_cum - lag(observed_cum),
    #        model_predicted = model_predicted - lag(model_predicted)
    # ) %>%
    # filter(!is.na(observed), !is.na(model_predicted)) %>%
    # ungroup()
  
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed,sd=sig_dist,log=TRUE))
  
  ll_final = if_else(is.na(ll) | is.infinite(ll), 10^6, ll)
  
  return(ll_final)
}

#### Debug ####
loglik_test <- function(
    npi_out2,
    rr_secondwave,
    npi_secondwave,
    npi_halfway,
    sig_dist
){
  
  #print(c(npi_out2, rr_secondwave, npi_secondwave, npi_halfway, sig_dist))
  
  #print("running")
  paras2 <- paras1
  paras2["npi_out_asian"] <- npi_out2
  paras2["npi_out_black"] <- npi_out2
  paras2["npi_out_latino"] <- npi_out2
  paras2["npi_out_other"] <- npi_out2
  paras2["npi_out_white"] <- npi_out2
  ## Reporting rate change
  paras2["rr_asian"] <- rr_secondwave*paras2["rr_asian"]
  paras2["rr_black"] <- rr_secondwave*paras2["rr_black"]
  paras2["rr_latino"] <- rr_secondwave*paras2["rr_latino"]
  paras2["rr_other"] <- rr_secondwave*paras2["rr_other"]
  paras2["rr_white"] <- rr_secondwave*paras2["rr_white"]
  
  ## set to IFRs over the whole period for each group (from seroprevalence data)
  paras1["sigma1_asian"] <- 0.013369570
  paras1["sigma1_black"] <- 0.014932698
  paras1["sigma1_latino"] <- 0.008206442
  paras1["sigma1_other"] <- 0.003015897
  paras1["sigma1_white"] <- 0.012578309
  
  init <- run_init
  
  times2 <- seq(t_start2,t_end2, by = 1)
  times3 <- seq(t_start3,t_end3, by = 1)
  
  # init for second run 
  init2 <- get %>% #out_calib1 %>% 
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist() 
  
  # second run 
  out_calib2 <- ode(y=init2, func = seir2_vax, times=times2, parms = paras2, method = "euler") %>%
    as.data.frame() %>%
    as_tibble()
  
  # init for third run
  # all out of lockdown to start
  # init3 <- out_calib2 %>%
  #   filter(time == max(time, na.rm = T)) %>%
  #   pivot_longer(-time) %>%
  #   separate(name, into = c("compartment", "race_ethnicity", "compliance"), "_") %>%
  #   mutate(compliance = if_else(is.na(compliance), "NL", compliance)) %>%
  #   pivot_wider(names_from = "compliance", values_from = "value") %>%
  #   mutate(total = NL + L,
  #          npi = case_when(
  #            race_ethnicity == "asian" ~ 0,
  #            race_ethnicity == "black" ~ 0,
  #            race_ethnicity == "latino" ~ 0,
  #            race_ethnicity == "other" ~ 0,
  #            race_ethnicity == "white" ~ 0
  #          ),
  #          new_NL = if_else(compartment %in% c("D", "RepCases", "Inc"), NL, (1-npi)*total),
  #          new_L = if_else(compartment %in% c("D", "RepCases", "Inc"), L, npi*total)
  #   ) %>%
  #   select(compartment, race_ethnicity,NL = new_NL, L = new_L) %>%
  #   pivot_longer(NL:L) %>%
  #   mutate(class = paste(compartment, race_ethnicity, name, sep = "_")) %>%
  #   mutate(class = gsub("_NL", "", class)) %>%
  #   filter(!is.na(value)) %>%
  #   select(class, value) %>%
  #   arrange(factor(class, levels = names_var)) %>%
  #   pivot_wider(names_from = "class", values_from = "value") %>%
  #   unlist()
  
  npi_a <- npi_weekly(npi_secondwave*npi_asian, npi_halfway, 2)
  npi_b <- npi_weekly(npi_secondwave*npi_black, npi_halfway, 2)
  npi_l <- npi_weekly(npi_secondwave*npi_latino, npi_halfway, 2)
  npi_o <- npi_weekly(npi_secondwave*npi_other, npi_halfway, 2)
  npi_w <- npi_weekly(npi_secondwave*npi_white, npi_halfway, 2)
  
  # init for third run 
  init3 <- out_calib2 %>%
    filter(time == max(time, na.rm = T)) %>%
    pivot_longer(-time) %>%
    select(name, value) %>%
    arrange(factor(name, levels = names_var)) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    unlist()
  
  # format 2nd run 
  out_calib2 <- out_calib2 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) %>%
    filter(day != t_end2) # filter final day - init for run 3
  
  # third run 
  out_calib3 <- ode(y=init3, func = seir2_vax, times=times3, parms = paras2, method = "euler") %>% # no parm change
    as.data.frame() %>%
    as_tibble()
  
  print(out_calib3)
  
  get3 <- out_calib3 %>%
    pivot_longer(-time) %>% 
    filter(substr(name, 1,1) == "D" | substr(name, 1,2) == "Re") %>% 
    separate(name, into = c("target", "race_ethnicity")) %>%
    mutate(target = case_when(
      target == "D" ~ "Deaths",
      target == "RepCases" ~ "Cases"
    )) %>% 
    select(day = time, target, race_ethnicity, model_predicted = value) %>% 
    inner_join(data, by = c("target", "race_ethnicity", "day")) 
  
  out_calib <- rbind(out_calib2, get3) %>%
    group_by(target, race_ethnicity) %>%
    arrange(day) %>%
    mutate(observed = observed_cum - lag(observed_cum),
           model_predicted = model_predicted - lag(model_predicted)
    ) %>%
    mutate(max_normal = max(observed, na.rm = T),
           observed = observed/max_normal,
           model_predicted = model_predicted/max_normal) %>%
    filter(!is.na(observed), !is.na(model_predicted)) %>%
    ungroup()
  # # mutate(phase = if_else(day < t_start2, "Phase 1", "Phase 2")) %>%
  # # group_by(target, race_ethnicity, phase) %>%
  # group_by(target, race_ethnicity) %>%
  # mutate(max_normal = max(observed_cum, na.rm = T),
  #        observed_cum = observed_cum/max_normal,
  #        model_predicted = model_predicted/max_normal) %>%
  # arrange(day) %>%
  # mutate(observed = observed_cum - lag(observed_cum),
  #        model_predicted = model_predicted - lag(model_predicted)
  # ) %>%
  # filter(!is.na(observed), !is.na(model_predicted)) %>%
  # ungroup()
  
  ll <- -sum(dnorm(x=out_calib$model_predicted,mean=out_calib$observed,sd=sig_dist,log=TRUE))
  
  print(ll)
  
  ll_final = if_else(is.na(ll), 10^6, ll)
  
  testing3 <- out_calib3 %>% filter(day == max(day))
  print(testing3)
  
  return(out_calib)
}


#### Phase 2 profiling ####


profile_function_NC2_ph1 <- function(i){
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms. 
  start1 = list(
    #mu = dt_run$Estimate[1],
    rr_asian_scale =  dt_run$Estimate[1]+0.001, rr_black_scale =  dt_run$Estimate[2],
    rr_other_scale =  dt_run$Estimate[3]-0.001, rr_white =  dt_run$Estimate[4]+0.001,
    npi_asian_scale = dt_run$Estimate[5]-0.001, npi_black_scale = dt_run$Estimate[6]+0.001,  
    npi_other_scale = dt_run$Estimate[7]+0.001, npi_white = dt_run$Estimate[8]+0.001,
    npi_out = dt_run$Estimate[9]+0.001,
    sig_dist = dt_run$Estimate[10])
  upper_list1 = upper_limit1
  lower_list1 = lower_limit1
  fixed_list1 = list(
    rr_asian_scale = parm_test, rr_black_scale = parm_test, 
    rr_other_scale = parm_test, rr_white = parm_test,
    npi_asian_scale = parm_test, npi_black_scale = parm_test, 
    npi_other_scale = parm_test, npi_white = parm_test,
    npi_out = parm_test,
    sig_dist = parm_test
  )
  
  startA <- start1[-i] ## remove the one we are fixing
  upper_listA <- upper_list1[-i]
  lower_listA <- lower_list1[-i]
  fixed_listA <- fixed_list1[i] ## make a list containing only the param we are fixing
  
  conf_threshold <- out_ll_A - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] + precision ## start at the true parameter value plus a step
  upper_limit_var <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    parm_test_right = dt_profile$Estimate[i]
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_right >= upper_limit_var){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listA = list(#mu = parm_test_right,
          rr_asian_scale = parm_test_right, rr_black_scale = parm_test_right,
          rr_other_scale = parm_test_right, rr_white = parm_test_right,
          npi_asian_scale = parm_test_right, npi_black_scale = parm_test_right, 
          npi_other_scale = parm_test_right, npi_white = parm_test_right,
          npi_out = parm_test_right,
          sig_dist = parm_test_right
        )
        fixed_listA <- fixed_listA[i]
        
        ## run mle and grab the new log likelihood
        m1_NC_A = mle2(minuslogl = loglik_nc2_phase1_short,
                       start = startA,
                       data = data_calib,
                       fixed = fixed_listA,
                       method="L-BFGS-B",
                       upper= upper_listA,
                       lower = lower_listA,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startA)))
        )

        
        right_ll_A <- -m1_NC_A@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        
        new_ll_right <- right_ll_A
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i] - precision ## start at the true parameter value minus a step
  lower_limit_var <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    parm_test_left = dt_profile$Estimate[i]
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit_var){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        
        
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listA = list(#mu = parm_test_right,
          rr_asian_scale = parm_test_left, rr_black_scale = parm_test_left,
          rr_other_scale = parm_test_left, rr_white = parm_test_left,
          npi_asian_scale = parm_test_left, npi_black_scale = parm_test_left, 
          npi_other_scale = parm_test_left, npi_white = parm_test_left,
          npi_out = parm_test_left,
          sig_dist = parm_test_left
        )
        fixed_listA <- fixed_listA[i]
        
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        m1_NC_A = mle2(minuslogl = loglik_nc2_phase1_short,
                       start = startA,
                       data = data_calib,
                       fixed = fixed_listA,
                       method="L-BFGS-B",
                       upper= upper_listA,
                       lower = lower_listA,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startA)))
        )
        
        left_ll_A <- -m1_NC_A@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        new_ll_left <- left_ll_A
        
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll_A,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_NC_2/", dt_profile$prof_param[i], "ph1.csv"))
  
}

profile_function_NC2_ph2 <- function(i){
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms. 
  start2 = list(
    npi_out2 =  dt_run2$Estimate[1] + 0.001,
    rr_secondwave = dt_run2$Estimate[2],
    npi_secondwave = dt_run2$Estimate[3],
    #npi_halfway = dt_run$Estimate[4],
    sig_dist = dt_run2$Estimate[4])
  upper_list2 = upper_limit2
  lower_list2 = lower_limit2
  fixed_list2 = list(
    npi_out2 = parm_test,
    rr_secondwave = parm_test,
    npi_secondwave = parm_test,
    #npi_halfway = parm_test,
    sig_dist = parm_test
  )
  
  startB <- start2[-i] ## remove the one we are fixing
  upper_listB <- upper_list2[-i]
  lower_listB <- lower_list2[-i]
  fixed_listB <- fixed_list2[i] ## make a list containing only the param we are fixing

  conf_threshold <- out_ll_B - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] + precision ## start close to the true parameter value
  upper_limit_var <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    parm_test_right = dt_profile$Estimate[i]
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_right >= upper_limit_var){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          npi_out2 = parm_test_right,
          rr_secondwave = parm_test_right,
          npi_secondwave = parm_test_right,
          #npi_halfway = parm_test_right,
          sig_dist = parm_test_right
        )
        fixed_listB <- fixed_listB[i]
        
        
        ## run mle and grab the new log likelihood
        m1_NC_B = mle2(minuslogl = loglik_nc2_phase2_short,
                       start = startB,
                       data = data_calib,
                       fixed = fixed_listB,
                       method="L-BFGS-B",
                       upper= upper_listB,
                       lower = lower_listB,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        right_ll_B <- -m1_NC_B@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        
        new_ll_right <- right_ll_B
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i] - precision ## start at the true parameter value minus a step
  lower_limit_var <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    parm_test_left = dt_profile$Estimate[i]
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit_var){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        
        
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          npi_out2 = parm_test_left,
          rr_secondwave = parm_test_left,
          npi_secondwave = parm_test_left,
          #npi_halfway = parm_test_left,
          sig_dist = parm_test_left
        )
        fixed_listB <- fixed_listB[i]
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        m1_NC_B = mle2(minuslogl = loglik_nc2_phase2_short,
                       start = startB,
                       data = data_calib,
                       fixed = fixed_listB,
                       method="L-BFGS-B",
                       upper= upper_listB,
                       lower = lower_listB,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        left_ll_B <- -m1_NC_B@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        new_ll_left <- left_ll_B
        
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll_B,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_NC_2/", dt_profile$prof_param[i], "ph2.csv"))
  
}



profile_function_CA2_ph1 <- function(i){
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms. 
  start1 = list(
    #mu = dt_run$Estimate[1],
    rr_asian_scale =  dt_run$Estimate[1]+0.001, rr_black_scale =  dt_run$Estimate[2],
    rr_latino_scale = dt_run$Estimate[3]+0.001,
    rr_other_scale =  dt_run$Estimate[4]-0.001, rr_white =  dt_run$Estimate[5]+0.001,
    npi_asian_scale = dt_run$Estimate[6]-0.001, npi_black_scale = dt_run$Estimate[7]+0.001,  
    npi_latino_scale = dt_run$Estimate[8] + 0.001,
    npi_other_scale = dt_run$Estimate[9]+0.001, npi_white = dt_run$Estimate[10]+0.001,
    npi_out = dt_run$Estimate[11]+0.001,
    sig_dist = dt_run$Estimate[12])
  upper_list1 = upper_limit1
  lower_list1 = lower_limit1
  fixed_list1 = list(
    rr_asian_scale = parm_test, rr_black_scale = parm_test, 
    rr_latino_scale = parm_test,
    rr_other_scale = parm_test, rr_white = parm_test,
    npi_asian_scale = parm_test, npi_black_scale = parm_test, 
    npi_latino_scale = parm_test,
    npi_other_scale = parm_test, npi_white = parm_test,
    npi_out = parm_test,
    sig_dist = parm_test
  )
  
  startA <- start1[-i] ## remove the one we are fixing
  upper_listA <- upper_list1[-i]
  lower_listA <- lower_list1[-i]
  fixed_listA <- fixed_list1[i] ## make a list containing only the param we are fixing
  
  conf_threshold <- out_ll_A - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] + precision ## start at the true parameter value plus a step
  upper_limit_var <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    parm_test_right = dt_profile$Estimate[i]
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_right >= upper_limit_var){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listA = list(#mu = parm_test_right,
          rr_asian_scale = parm_test_right, rr_black_scale = parm_test_right,
          rr_latino_scale = parm_test_right,
          rr_other_scale = parm_test_right, rr_white = parm_test_right,
          npi_asian_scale = parm_test_right, npi_black_scale = parm_test_right, 
          npi_latino_scale = parm_test_right,
          npi_other_scale = parm_test_right, npi_white = parm_test_right,
          npi_out = parm_test_right,
          sig_dist = parm_test_right
        )
        fixed_listA <- fixed_listA[i]
        
        
        ## run mle and grab the new log likelihood
        m1_CA_A = mle2(minuslogl = loglik_ca2_phase1_short,
                       start = startA,
                       data = data_calib,
                       fixed = fixed_listA,
                       method="L-BFGS-B",
                       upper= upper_listA,
                       lower = lower_listA,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startA)))
        )
        
        right_ll_A <- -m1_CA_A@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        
        new_ll_right <- right_ll_A
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i] - precision ## start at the true parameter value minus a step
  lower_limit_var <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    parm_test_left = dt_profile$Estimate[i]
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit_var){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        
        
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listA = list(#mu = parm_test_right,
          rr_asian_scale = parm_test_left, rr_black_scale = parm_test_left,
          rr_latino_scale = parm_test_left,
          rr_other_scale = parm_test_left, rr_white = parm_test_left,
          npi_asian_scale = parm_test_left, npi_black_scale = parm_test_left, 
          npi_latino_scale = parm_test_left,
          npi_other_scale = parm_test_left, npi_white = parm_test_left,
          npi_out = parm_test_left,
          sig_dist = parm_test_left
        )
        fixed_listA <- fixed_listA[i]
        
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        m1_CA_A = mle2(minuslogl = loglik_ca2_phase1_short,
                       start = startA,
                       data = data_calib,
                       fixed = fixed_listA,
                       method="L-BFGS-B",
                       upper= upper_listA,
                       lower = lower_listA,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startA)))
        )
        
        left_ll_A <- -m1_CA_A@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        new_ll_left <- left_ll_A
        
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll_A,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_CA_2/", dt_profile$prof_param[i], "ph1.csv"))
  
}

profile_function_CA2_ph2 <- function(i){
  ## grab the parameter estimate from the maximum likelihood
  parm_test = dt_profile$Estimate[i]
  
  ## set up the start and upper/lowers for the first mle. If any parameters are on the boundary, shift slightly so that the start parms are not boundary parms. 
  start2 = list(
    npi_out2 =  dt_run2$Estimate[1],
    rr_secondwave = dt_run2$Estimate[2],
    npi_secondwave = dt_run2$Estimate[3],
    sig_dist = dt_run2$Estimate[4])
  upper_list2 = upper_limit2
  lower_list2 = lower_limit2
  fixed_list2 = list(
    npi_out2 = parm_test,
    rr_secondwave = parm_test,
    npi_secondwave = parm_test,
    sig_dist = parm_test
  )
  
  startB <- start2[-i] ## remove the one we are fixing
  upper_listB <- upper_list2[-i]
  lower_listB <- lower_list2[-i]
  fixed_listB <- fixed_list2[i] ## make a list containing only the param we are fixing
  
  conf_threshold <- out_ll_B - 1.92 ## grab the confidence threshold for log likelihood
  
  ## get right threshold
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_right = dt_profile$Estimate[i] +  precision ## start at the true parameter value plus a step
  upper_limit_var <- dt_profile$upper[i]
  
  print("Running right MLE")
  
  if (round(dt_profile$Estimate[i],4) == round(dt_profile$upper[i],4)){
    parm_test_right = dt_profile$Estimate[i]
    new_ll_right = -10^3
  }
  else{
    while (precision > 10e-5) {
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      if (parm_test_right >= upper_limit_var){
        new_ll_right <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          npi_out2 = parm_test_right,
          rr_secondwave = parm_test_right,
          npi_secondwave = parm_test_right,
          #npi_halfway = parm_test_right,
          sig_dist = parm_test_right
        )
        fixed_listB <- fixed_listB[i]
        
        
        ## run mle and grab the new log likelihood
        m1_CA_B = mle2(minuslogl = loglik_ca2_phase2_short,
                       start = startB,
                       data = data_calib,
                       fixed = fixed_listB,
                       method="L-BFGS-B",
                       upper= upper_listB,
                       lower = lower_listB,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        right_ll_B <- -m1_CA_B@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        
        new_ll_right <- right_ll_B
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_right == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_right > conf_threshold){
        parm_test_right = parm_test_right + precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_right < conf_threshold){
        parm_test_right = parm_test_right - precision
        precision = precision/10
        parm_test_right = parm_test_right + precision
      }
    }
    
  }
  
  
  ## First, set the starting steps away from the parameter
  precision = 0.1
  parm_test_left = dt_profile$Estimate[i]  - precision ## start at the true parameter value minus first step
  lower_limit_var <- dt_profile$lower[i]
  
  print("Running left MLE")
  
  if (round(dt_profile$Estimate[i], 4) == round(dt_profile$lower[i], 4)){
    parm_test_left = dt_profile$Estimate[i]
    new_ll_left = -10^3
  }
  else{
    while (precision > 10e-5) {
      
      ## IF the test param has exceeded the allowed window for the parameter, we will not run the mle, and instead go back a step and adjust precision (treat it like it is < conf_threshold)
      
      if (parm_test_left <= lower_limit_var){
        new_ll_left <- -10^3 ## return dummy log likelihood that will always be outside the threshold
      } ## OTHERWISE we will run the mle
      else{
        
        
        #parm_test <- parm_test + precision
        ## update the list of fixed parms
        fixed_listB = list(
          npi_out2 = parm_test_left,
          rr_secondwave = parm_test_left,
          npi_secondwave = parm_test_left,
          #npi_halfway = parm_test_left,
          sig_dist = parm_test_left
        )
        fixed_listB <- fixed_listB[i]
        
       #print(parm_test_left)
        
        ## run mle and grab the new log likelihood
        ## run m1_A for starting param
        m1_CA_B = mle2(minuslogl = loglik_ca2_phase2_short,
                       start = startB,
                       data = data_calib,
                       fixed = fixed_listB,
                       method="L-BFGS-B",
                       upper= upper_listB,
                       lower = lower_listB,
                       control=list(maxit=1000, trace=TRUE, parscale=abs(unlist(startB)))
        )
        
        left_ll_B <- -m1_CA_B@min ## grab the base likelihood for phase 1 (remember, we fit on the negative log likelihood, so this should have a minus sign in front)
        new_ll_left <- left_ll_B
        
      }
      
      ## If it equals the confidence threshold, stop
      if (new_ll_left == conf_threshold){
        break
      }
      ## if it's above the threshold, shift the param again to the right
      else if(new_ll_left > conf_threshold){
        parm_test_left = parm_test_left - precision}
      ## if it's below the threshold, shift the param to the right but by less
      else if(new_ll_left < conf_threshold){
        parm_test_left = parm_test_left + precision
        precision = precision/10
        parm_test_left = parm_test_left - precision
      }
    }
  }
  
  out_data <- tibble(param = dt_profile$prof_param[i],
                     loglik = out_ll_B,
                     ll_05 = new_ll_left,
                     ll_95 = new_ll_right,
                     conf_05 = parm_test_left,
                     conf_95 = parm_test_right)
  
  write_csv(out_data, paste("../data/profile_CA_2/", dt_profile$prof_param[i], "ph2.csv"))
  
}


#### SEIR2 with vaccines ####
# For california
seir2_vax <- function(t, y, pars){
  # Parameters (fixed)
  epsilon <- pars["epsilon"] ## latent period
  omega <- pars["omega"] ## waning of immunity
  gamma <- pars["gamma"] ## infectious period
  ir_2 <- pars["ir_2"] ## scale infection risk for prior immunity
  ifr_2 <- pars["ifr_2"] ## scale infection fatality rate for prior immunity
  npi_efficacy <- pars["npi_efficacy"]
  
  ## Parameters (fit)
  mu <- pars["mu"] ## risk of infection
  
  rr_asian <- pars["rr_asian"] # reporting rate
  rr_black <- pars["rr_black"]
  rr_latino <- pars["rr_latino"]
  rr_other <- pars["rr_other"]
  rr_white <- pars["rr_white"]
  
  npi_out_asian <- pars["npi_out_asian"] ## rate of leaving npi's
  npi_out_black <- pars["npi_out_black"]
  npi_out_latino <- pars["npi_out_latino"]
  npi_out_other <- pars["npi_out_other"]
  npi_out_white <- pars["npi_out_white"]
  
  sigma1_asian <- pars["sigma1_asian"] ## infection fatality rate
  sigma1_black <- pars["sigma1_black"] 
  sigma1_latino <- pars["sigma1_latino"] 
  sigma1_other <- pars["sigma1_other"] 
  sigma1_white <- pars["sigma1_white"] 
  
  # State variables
  ## Asian
  S1_asian <- y[1]
  E1_asian <- y[2]
  I1_asian <- y[3]
  R1_asian <- y[4]
  S2_asian <- y[5]
  E2_asian <- y[6]
  I2_asian <- y[7]
  R2_asian <- y[8]
  
  S1_asian_L <- y[9]
  E1_asian_L <- y[10]
  I1_asian_L <- y[11]
  R1_asian_L <- y[12]
  S2_asian_L <- y[13]
  E2_asian_L <- y[14]
  I2_asian_L <- y[15]
  R2_asian_L <- y[16]
  
  D_asian <- y[17]
  Inc_asian <- y[18]
  RepCases_asian <- y[19]
  
  ## Black
  S1_black <- y[20]
  E1_black <- y[21]
  I1_black <- y[22]
  R1_black <- y[23]
  S2_black <- y[24]
  E2_black <- y[25]
  I2_black <- y[26]
  R2_black <- y[27]
  
  S1_black_L <- y[28]
  E1_black_L <- y[29]
  I1_black_L <- y[30]
  R1_black_L <- y[31]
  S2_black_L <- y[32]
  E2_black_L <- y[33]
  I2_black_L <- y[34]
  R2_black_L <- y[35]
  
  D_black <- y[36]
  Inc_black <- y[37]
  RepCases_black <- y[38]
  
  ## Latino
  S1_latino <- y[39]
  E1_latino <- y[40]
  I1_latino <- y[41]
  R1_latino <- y[42]
  S2_latino <- y[43]
  E2_latino <- y[44]
  I2_latino <- y[45]
  R2_latino <- y[46]
  
  S1_latino_L <- y[47]
  E1_latino_L <- y[48]
  I1_latino_L <- y[49]
  R1_latino_L <- y[50]
  S2_latino_L <- y[51]
  E2_latino_L <- y[52]
  I2_latino_L <- y[53]
  R2_latino_L <- y[54]
  
  D_latino <- y[55]
  Inc_latino <- y[56]
  RepCases_latino <- y[57]
  
  ## Other
  S1_other <- y[58]
  E1_other <- y[59]
  I1_other <- y[60]
  R1_other <- y[61]
  S2_other <- y[62]
  E2_other <- y[63]
  I2_other <- y[64]
  R2_other <- y[65]
  
  S1_other_L <- y[66]
  E1_other_L <- y[67]
  I1_other_L <- y[68]
  R1_other_L <- y[69]
  S2_other_L <- y[70]
  E2_other_L <- y[71]
  I2_other_L <- y[72]
  R2_other_L <- y[73]
  
  D_other <- y[74]
  Inc_other <- y[75]
  RepCases_other <- y[76]
  
  ## White
  S1_white <- y[77]
  E1_white <- y[78]
  I1_white <- y[79]
  R1_white <- y[80]
  S2_white <- y[81]
  E2_white <- y[82]
  I2_white <- y[83]
  R2_white <- y[84]
  
  S1_white_L <- y[85]
  E1_white_L <- y[86]
  I1_white_L <- y[87]
  R1_white_L <- y[88]
  S2_white_L <- y[89]
  E2_white_L <- y[90]
  I2_white_L <- y[91]
  R2_white_L <- y[92]
  
  D_white <- y[93]
  Inc_white <- y[94]
  RepCases_white <- y[95]
  
  ## calculations
  I_asian_eff <- I1_asian + I2_asian + (1-npi_efficacy)*(I1_asian_L + I2_asian_L)
  I_black_eff <- I1_black + I2_black + (1-npi_efficacy)*(I1_black_L + I2_black_L)
  I_latino_eff <- I1_latino + I1_latino + (1-npi_efficacy)*(I1_latino_L + I2_latino_L)
  I_other_eff <- I1_other + I2_other + (1-npi_efficacy)*(I1_other_L + I2_other_L)
  I_white_eff <- I1_white + I2_white + (1-npi_efficacy)*(I1_white_L + I2_white_L)
  
  P_asian_eff <- S1_asian + E1_asian + I1_asian + R1_asian + S2_asian + E2_asian + I2_asian + R2_asian  
  P_black_eff <- S1_black + E1_black + I1_black + R1_black + S2_black + E2_black + I2_black + R2_black 
  P_latino_eff <- S1_latino + E1_latino + I1_latino + R1_latino + S2_latino + E2_latino + I2_latino + R2_latino 
  P_other_eff <- S1_other + E1_other + I1_other + R1_other + S2_other + E2_other + I2_other + R2_other 
  P_white_eff <- S1_white + E1_white + I1_white + R1_white + S2_white + E2_white + I2_white + R2_white 
  
  P_asian_total <- P_asian_eff + S1_asian_L + E1_asian_L + I1_asian_L + R1_asian_L + S2_asian_L + E2_asian_L + I2_asian_L + R2_asian_L
  P_black_total <- P_black_eff + S1_black_L + E1_black_L + I1_black_L + R1_black_L + S2_black_L + E2_black_L + I2_black_L + R2_black_L  
  P_latino_total <- P_latino_eff + S1_latino_L + E1_latino_L + I1_latino_L + R1_latino_L + S2_latino_L + E2_latino_L + I2_latino_L + R2_latino_L  
  P_other_total <- P_other_eff + S1_other_L + E1_other_L + I1_other_L + R1_other_L + S2_other_L + E2_other_L + I2_other_L + R2_other_L  
  P_white_total <- P_white_eff + S1_white_L + E1_white_L + I1_white_L + R1_white_L + S2_white_L + E2_white_L + I2_white_L + R2_white_L
  
  ## CHANGING HERE TO P total
  prop_asian_eff <- I_asian_eff/P_asian_total
  prop_black_eff <- I_black_eff/P_black_total
  prop_latino_eff <- I_latino_eff/P_latino_total
  prop_other_eff <- I_other_eff/P_other_total
  prop_white_eff <- I_white_eff/P_white_total
  
  ## contacts 
  asian_c_asian <- aca(t)
  asian_c_black <- acb(t)
  asian_c_latino <- acl(t)
  asian_c_other <- aco(t)
  asian_c_white <- acw(t)
  
  black_c_asian <- bca(t)
  black_c_black <- bcb(t)
  black_c_latino <- bcl(t)
  black_c_other <- bco(t)
  black_c_white <- bcw(t)
  
  latino_c_asian <- lca(t)
  latino_c_black <- lcb(t)
  latino_c_latino <- lcl(t)
  latino_c_other <- lco(t)
  latino_c_white <- lcw(t)
  
  other_c_asian <- oca(t)
  other_c_black <- ocb(t)
  other_c_latino <- ocl(t)
  other_c_other <- oco(t)
  other_c_white <- ocw(t)
  
  white_c_asian <- wca(t)
  white_c_black <- wcb(t)
  white_c_latino <- wcl(t)
  white_c_other <- wco(t)
  white_c_white <- wcw(t)
  
  # The vaccination rates are based on data referenced to the whole population. In this model, only susceptible and recovered can get vaccinated (S1, R1, S2, R2). 
  # Thus, the pool to be vaccinated needs to be adjusted. If v is the percent daily that is observed in the data, then
  # v*P = v_adj * (S1 + R1 + S2 + R2)
  # so that v_adj = v*P/(S1 + R1 + S2 + R2)
  
  v_scalar_asian <- P_asian_total/(S1_asian + S1_asian_L + R1_asian + R1_asian_L + S2_asian + S2_asian_L + R2_asian  + R2_asian_L)
  v_scalar_black <- P_black_total/(S1_black + S1_black_L + R1_black + R1_black_L + S2_black + S2_black_L + R2_black  + R2_black_L)
  v_scalar_latino <- P_latino_total/(S1_latino + S1_latino_L + R1_latino + R1_latino_L + S2_latino + S2_latino_L + R2_latino  + R2_latino_L)
  v_scalar_other <- P_other_total/(S1_other + S1_other_L + R1_other + R1_other_L + S2_other + S2_other_L + R2_other  + R2_other_L)
  v_scalar_white <- P_white_total/(S1_white + S1_white_L + R1_white + R1_white_L + S2_white + S2_white_L + R2_white  + R2_white_L)
  
  nu_asian <- min(0.85*v_scalar_asian*vax_a(t),1)
  nu_black <- min(0.85*v_scalar_black*vax_b(t), 1)
  nu_latino <- min(0.85*v_scalar_latino*vax_l(t), 1)
  nu_other <- min(0.85*v_scalar_other*vax_o(t), 1)
  nu_white <- min(0.85*v_scalar_white*vax_w(t), 1)
  
  var <- variant(t)
  
  lambda1_asian <- var*mu*(aca(t)*prop_asian_eff + acb(t)*prop_black_eff + acl(t)*prop_latino_eff + aco(t)*prop_other_eff + acw(t)*prop_white_eff)
  lambda2_asian <- ir_2*lambda1_asian
  
  lambda1_black <- var*mu*(bca(t)*prop_asian_eff+ bcb(t)*prop_black_eff + bcl(t)*prop_latino_eff + bco(t)*prop_other_eff + bcw(t)*prop_white_eff)
  lambda2_black<- ir_2*lambda1_black
  
  lambda1_latino <- var*mu*(lca(t)*prop_asian_eff + lcb(t)*prop_black_eff + lcl(t)*prop_latino_eff + lco(t)*prop_other_eff + lcw(t)*prop_white_eff)
  lambda2_latino<- ir_2*lambda1_latino
  
  lambda1_other <- var*mu*(oca(t)*prop_asian_eff+ ocb(t)*prop_black_eff + ocl(t)*prop_latino_eff + oco(t)*prop_other_eff + ocw(t)*prop_white_eff)
  lambda2_other <- ir_2*lambda1_other
  
  lambda1_white <- var*mu*(wca(t)*prop_asian_eff+ wcb(t)*prop_black_eff + wcl(t)*prop_latino_eff + wco(t)*prop_other_eff + wcw(t)*prop_white_eff)
  lambda2_white <- ir_2*lambda1_white
  
  ## Updates
  dS1_asian <- - lambda1_asian*S1_asian + npi_out_asian*S1_asian_L - nu_asian*S1_asian
  dE1_asian <- + lambda1_asian*S1_asian - epsilon*E1_asian  + npi_out_asian*E1_asian_L 
  dI1_asian <- + epsilon*E1_asian - gamma*I1_asian + npi_out_asian*I1_asian_L
  dR1_asian <- + (1-sigma1_asian)*gamma*I1_asian - omega*R1_asian + npi_out_asian*R1_asian_L - nu_asian*R1_asian
  dS2_asian <- + omega*R1_asian - lambda2_asian*S2_asian + omega*R2_asian + npi_out_asian*S2_asian_L - nu_asian*S2_asian
  dE2_asian <- + lambda2_asian*S2_asian - epsilon*E2_asian + npi_out_asian*E2_asian_L
  dI2_asian <- + epsilon*E2_asian - gamma*I2_asian + npi_out_asian*I2_asian_L
  dR2_asian <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian - omega*R2_asian + npi_out_asian*R2_asian_L + nu_asian*(S1_asian + R1_asian + S2_asian)
  
  dS1_asian_L <-  - npi_out_asian*S1_asian_L - nu_asian*S1_asian_L - (1-npi_efficacy)*lambda1_asian*S1_asian_L
  dE1_asian_L <-  - epsilon*E1_asian_L  - npi_out_asian*E1_asian_L + (1-npi_efficacy)*lambda1_asian*S1_asian_L
  dI1_asian_L <- + epsilon*E1_asian_L - gamma*I1_asian_L - npi_out_asian*I1_asian_L
  dR1_asian_L <- + (1-sigma1_asian)*gamma*I1_asian_L - omega*R1_asian_L  - npi_out_asian*R1_asian_L - nu_asian*R1_asian_L
  dS2_asian_L <- + omega*R1_asian_L + omega*R2_asian_L - npi_out_asian*S2_asian_L - nu_asian*S2_asian_L - (1-npi_efficacy)*lambda2_asian*S2_asian_L
  dE2_asian_L <- - epsilon*E2_asian_L - npi_out_asian*E2_asian_L + (1-npi_efficacy)*lambda2_asian*S2_asian_L
  dI2_asian_L <- + epsilon*E2_asian_L - gamma*I2_asian_L  - npi_out_asian*I2_asian_L
  dR2_asian_L <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian_L - omega*R2_asian_L - npi_out_asian*R2_asian_L + nu_asian*(S1_asian_L + R1_asian_L + S2_asian_L)
  
  dD_asian <- sigma1_asian*gamma*(I1_asian + I1_asian_L) + ifr_2*sigma1_asian*gamma*(I2_asian + I2_asian_L)
  dInc_asian <- + epsilon*(E1_asian + E2_asian + E1_asian_L + E2_asian_L)
  dRepCases_asian <- rr_asian*epsilon*(E1_asian+ E2_asian + E1_asian_L + E2_asian_L)
  
  dS1_black <- - lambda1_black*S1_black + npi_out_black*S1_black_L - nu_black*S1_black
  dE1_black <- + lambda1_black*S1_black - epsilon*E1_black  + npi_out_black*E1_black_L 
  dI1_black <- + epsilon*E1_black - gamma*I1_black + npi_out_black*I1_black_L
  dR1_black <- + (1-sigma1_black)*gamma*I1_black - omega*R1_black + npi_out_black*R1_black_L - nu_black*R1_black
  dS2_black <- + omega*R1_black - lambda2_black*S2_black + omega*R2_black + npi_out_black*S2_black_L - nu_black*S2_black
  dE2_black <- + lambda2_black*S2_black - epsilon*E2_black + npi_out_black*E2_black_L
  dI2_black <- + epsilon*E2_black - gamma*I2_black + npi_out_black*I2_black_L
  dR2_black <- + (1-ifr_2*sigma1_black)*gamma*I2_black - omega*R2_black + npi_out_black*R2_black_L + nu_black*(S1_black + R1_black + S2_black)
  
  dS1_black_L <-  - npi_out_black*S1_black_L - nu_black*S1_black_L - (1-npi_efficacy)*lambda1_black*S1_black_L
  dE1_black_L <-  - epsilon*E1_black_L  - npi_out_black*E1_black_L + (1-npi_efficacy)*lambda1_black*S1_black_L
  dI1_black_L <- + epsilon*E1_black_L - gamma*I1_black_L - npi_out_black*I1_black_L
  dR1_black_L <- + (1-sigma1_black)*gamma*I1_black_L - omega*R1_black_L  - npi_out_black*R1_black_L - nu_black*R1_black_L
  dS2_black_L <- + omega*R1_black_L + omega*R2_black_L - npi_out_black*S2_black_L - nu_black*S2_black_L - (1-npi_efficacy)*lambda2_black*S2_black_L
  dE2_black_L <- - epsilon*E2_black_L - npi_out_black*E2_black_L + (1-npi_efficacy)*lambda2_black*S2_black_L
  dI2_black_L <- + epsilon*E2_black_L - gamma*I2_black_L  - npi_out_black*I2_black_L
  dR2_black_L <- + (1-ifr_2*sigma1_black)*gamma*I2_black_L - omega*R2_black_L - npi_out_black*R2_black_L + nu_black*(S1_black_L + R1_black_L + S2_black_L)
  
  dD_black <- sigma1_black*gamma*(I1_black + I1_black_L) + ifr_2*sigma1_black*gamma*(I2_black + I2_black_L)
  dInc_black <- + epsilon*(E1_black + E2_black + E1_black_L + E2_black_L)
  dRepCases_black <- rr_black*epsilon*(E1_black+ E2_black + E1_black_L + E2_black_L)
  
  dS1_latino <- - lambda1_latino*S1_latino + npi_out_latino*S1_latino_L - nu_latino*S1_latino
  dE1_latino <- + lambda1_latino*S1_latino - epsilon*E1_latino  + npi_out_latino*E1_latino_L 
  dI1_latino <- + epsilon*E1_latino - gamma*I1_latino + npi_out_latino*I1_latino_L
  dR1_latino <- + (1-sigma1_latino)*gamma*I1_latino - omega*R1_latino + npi_out_latino*R1_latino_L - nu_latino*R1_latino
  dS2_latino <- + omega*R1_latino - lambda2_latino*S2_latino + omega*R2_latino + npi_out_latino*S2_latino_L - nu_latino*S2_latino
  dE2_latino <- + lambda2_latino*S2_latino - epsilon*E2_latino + npi_out_latino*E2_latino_L
  dI2_latino <- + epsilon*E2_latino - gamma*I2_latino + npi_out_latino*I2_latino_L
  dR2_latino <- + (1-ifr_2*sigma1_latino)*gamma*I2_latino - omega*R2_latino + npi_out_latino*R2_latino_L + nu_latino*(S1_latino + R1_latino + S2_latino)
  
  dS1_latino_L <-  - npi_out_latino*S1_latino_L - nu_latino*S1_latino_L - (1-npi_efficacy)*lambda1_latino*S1_latino_L
  dE1_latino_L <-  - epsilon*E1_latino_L  - npi_out_latino*E1_latino_L + (1-npi_efficacy)*lambda1_latino*S1_latino_L
  dI1_latino_L <- + epsilon*E1_latino_L - gamma*I1_latino_L - npi_out_latino*I1_latino_L
  dR1_latino_L <- + (1-sigma1_latino)*gamma*I1_latino_L - omega*R1_latino_L  - npi_out_latino*R1_latino_L - nu_latino*R1_latino_L
  dS2_latino_L <- + omega*R1_latino_L + omega*R2_latino_L - npi_out_latino*S2_latino_L - nu_latino*S2_latino_L - (1-npi_efficacy)*lambda2_latino*S2_latino_L
  dE2_latino_L <- - epsilon*E2_latino_L - npi_out_latino*E2_latino_L + (1-npi_efficacy)*lambda2_latino*S2_latino_L
  dI2_latino_L <- + epsilon*E2_latino_L - gamma*I2_latino_L  - npi_out_latino*I2_latino_L
  dR2_latino_L <- + (1-ifr_2*sigma1_latino)*gamma*I2_latino_L - omega*R2_latino_L - npi_out_latino*R2_latino_L + nu_latino*(S1_latino_L + R1_latino_L + S2_latino_L)
  
  dD_latino <- sigma1_latino*gamma*(I1_latino + I1_latino_L) + ifr_2*sigma1_latino*gamma*(I2_latino + I2_latino_L)
  dInc_latino <- + epsilon*(E1_latino + E2_latino + E1_latino_L + E2_latino_L)
  dRepCases_latino <- rr_latino*epsilon*(E1_latino+ E2_latino + E1_latino_L + E2_latino_L)
  
  dS1_other <- - lambda1_other*S1_other + npi_out_other*S1_other_L - nu_other*S1_other
  dE1_other <- + lambda1_other*S1_other - epsilon*E1_other  + npi_out_other*E1_other_L 
  dI1_other <- + epsilon*E1_other - gamma*I1_other + npi_out_other*I1_other_L
  dR1_other <- + (1-sigma1_other)*gamma*I1_other - omega*R1_other + npi_out_other*R1_other_L - nu_other*R1_other
  dS2_other <- + omega*R1_other - lambda2_other*S2_other + omega*R2_other + npi_out_other*S2_other_L - nu_other*S2_other
  dE2_other <- + lambda2_other*S2_other - epsilon*E2_other + npi_out_other*E2_other_L
  dI2_other <- + epsilon*E2_other - gamma*I2_other + npi_out_other*I2_other_L
  dR2_other <- + (1-ifr_2*sigma1_other)*gamma*I2_other - omega*R2_other + npi_out_other*R2_other_L + nu_other*(S1_other + R1_other + S2_other)
  
  dS1_other_L <-  - npi_out_other*S1_other_L - nu_other*S1_other_L - (1-npi_efficacy)*lambda1_other*S1_other_L
  dE1_other_L <-  - epsilon*E1_other_L  - npi_out_other*E1_other_L + (1-npi_efficacy)*lambda1_other*S1_other_L
  dI1_other_L <- + epsilon*E1_other_L - gamma*I1_other_L - npi_out_other*I1_other_L
  dR1_other_L <- + (1-sigma1_other)*gamma*I1_other_L - omega*R1_other_L  - npi_out_other*R1_other_L - nu_other*R1_other_L
  dS2_other_L <- + omega*R1_other_L + omega*R2_other_L - npi_out_other*S2_other_L - nu_other*S2_other_L - (1-npi_efficacy)*lambda2_other*S2_other_L
  dE2_other_L <- - epsilon*E2_other_L - npi_out_other*E2_other_L + (1-npi_efficacy)*lambda2_other*S2_other_L
  dI2_other_L <- + epsilon*E2_other_L - gamma*I2_other_L  - npi_out_other*I2_other_L
  dR2_other_L <- + (1-ifr_2*sigma1_other)*gamma*I2_other_L - omega*R2_other_L - npi_out_other*R2_other_L + nu_other*(S1_other_L + R1_other_L + S2_other_L)
  
  dD_other <- sigma1_other*gamma*(I1_other + I1_other_L) + ifr_2*sigma1_other*gamma*(I2_other + I2_other_L)
  dInc_other <- + epsilon*(E1_other + E2_other + E1_other_L + E2_other_L)
  dRepCases_other <- rr_other*epsilon*(E1_other+ E2_other + E1_other_L + E2_other_L)
  
  dS1_white <- - lambda1_white*S1_white + npi_out_white*S1_white_L - nu_white*S1_white
  dE1_white <- + lambda1_white*S1_white - epsilon*E1_white  + npi_out_white*E1_white_L 
  dI1_white <- + epsilon*E1_white - gamma*I1_white + npi_out_white*I1_white_L
  dR1_white <- + (1-sigma1_white)*gamma*I1_white - omega*R1_white + npi_out_white*R1_white_L - nu_white*R1_white
  dS2_white <- + omega*R1_white - lambda2_white*S2_white + omega*R2_white + npi_out_white*S2_white_L - nu_white*S2_white
  dE2_white <- + lambda2_white*S2_white - epsilon*E2_white + npi_out_white*E2_white_L
  dI2_white <- + epsilon*E2_white - gamma*I2_white + npi_out_white*I2_white_L
  dR2_white <- + (1-ifr_2*sigma1_white)*gamma*I2_white - omega*R2_white + npi_out_white*R2_white_L + nu_white*(S1_white + R1_white + S2_white)
  
  dS1_white_L <-  - npi_out_white*S1_white_L - nu_white*S1_white_L - (1-npi_efficacy)*lambda1_white*S1_white_L
  dE1_white_L <-  - epsilon*E1_white_L  - npi_out_white*E1_white_L + (1-npi_efficacy)*lambda1_white*S1_white_L
  dI1_white_L <- + epsilon*E1_white_L - gamma*I1_white_L - npi_out_white*I1_white_L
  dR1_white_L <- + (1-sigma1_white)*gamma*I1_white_L - omega*R1_white_L  - npi_out_white*R1_white_L - nu_white*R1_white_L
  dS2_white_L <- + omega*R1_white_L + omega*R2_white_L - npi_out_white*S2_white_L - nu_white*S2_white_L - (1-npi_efficacy)*lambda2_white*S2_white_L
  dE2_white_L <- - epsilon*E2_white_L - npi_out_white*E2_white_L + (1-npi_efficacy)*lambda2_white*S2_white_L
  dI2_white_L <- + epsilon*E2_white_L - gamma*I2_white_L  - npi_out_white*I2_white_L
  dR2_white_L <- + (1-ifr_2*sigma1_white)*gamma*I2_white_L - omega*R2_white_L - npi_out_white*R2_white_L + nu_white*(S1_white_L + R1_white_L + S2_white_L)
  
  dD_white <- sigma1_white*gamma*(I1_white + I1_white_L) + ifr_2*sigma1_white*gamma*(I2_white + I2_white_L)
  dInc_white <- + epsilon*(E1_white + E2_white + E1_white_L + E2_white_L)
  dRepCases_white <- rr_white*epsilon*(E1_white+ E2_white + E1_white_L + E2_white_L)
  
  # Return list of gradients
  list(c(dS1_asian,dE1_asian,dI1_asian,dR1_asian,
         dS2_asian,dE2_asian,dI2_asian,dR2_asian,
         
         dS1_asian_L,dE1_asian_L,dI1_asian_L,dR1_asian_L,
         dS2_asian_L,dE2_asian_L,dI2_asian_L,dR2_asian_L,
         
         dD_asian, dInc_asian, dRepCases_asian,
         
         dS1_black,dE1_black,dI1_black,dR1_black,
         dS2_black,dE2_black,dI2_black,dR2_black,
         
         dS1_black_L,dE1_black_L,dI1_black_L,dR1_black_L,
         dS2_black_L,dE2_black_L,dI2_black_L,dR2_black_L,
         
         dD_black, dInc_black, dRepCases_black,
         
         dS1_latino,dE1_latino,dI1_latino,dR1_latino,
         dS2_latino,dE2_latino,dI2_latino,dR2_latino,
         
         dS1_latino_L,dE1_latino_L,dI1_latino_L,dR1_latino_L,
         dS2_latino_L,dE2_latino_L,dI2_latino_L,dR2_latino_L,
         
         dD_latino, dInc_latino, dRepCases_latino,
         
         dS1_other,dE1_other,dI1_other,dR1_other,
         dS2_other,dE2_other,dI2_other,dR2_other,
         
         dS1_other_L,dE1_other_L,dI1_other_L,dR1_other_L,
         dS2_other_L,dE2_other_L,dI2_other_L,dR2_other_L,
         
         dD_other, dInc_other, dRepCases_other,
         
         dS1_white,dE1_white,dI1_white,dR1_white,
         dS2_white,dE2_white,dI2_white,dR2_white,
         
         dS1_white_L,dE1_white_L,dI1_white_L,dR1_white_L,
         dS2_white_L,dE2_white_L,dI2_white_L,dR2_white_L,
         
         dD_white, dInc_white, dRepCases_white))
}


#### SEIR2 with vaccines ####
# for north carolina
seir2_nolatino_vax <- function(t, y, pars){
  # Parameters (fixed)
  epsilon <- pars["epsilon"] ## latent period
  omega <- pars["omega"] ## waning of immunity
  gamma <- pars["gamma"] ## infectious period
  ir_2 <- pars["ir_2"] ## scale infection risk for prior immunity
  ifr_2 <- pars["ifr_2"] ## scale infection fatality rate for prior immunity
  npi_efficacy <- pars["npi_efficacy"] ## npi efficacy (0-1)
  
  ## Parameters (fit)
  mu <- pars["mu"] ## risk of infection
  
  rr_asian <- pars["rr_asian"] # reporting rate
  rr_black <- pars["rr_black"]
  rr_other <- pars["rr_other"]
  rr_white <- pars["rr_white"]
  
  npi_out_asian <- pars["npi_out_asian"]
  npi_out_black <- pars["npi_out_black"]
  npi_out_other <- pars["npi_out_other"]
  npi_out_white <- pars["npi_out_white"]
  
  sigma1_asian <- pars["sigma1_asian"] ## infection fatality rate
  sigma1_black <- pars["sigma1_black"] 
  sigma1_other <- pars["sigma1_other"] 
  sigma1_white <- pars["sigma1_white"] 
  
  # State variables
  ## Asian
  S1_asian <- y[1]
  E1_asian <- y[2]
  I1_asian <- y[3]
  R1_asian <- y[4]
  S2_asian <- y[5]
  E2_asian <- y[6]
  I2_asian <- y[7]
  R2_asian <- y[8]
  
  S1_asian_L <- y[9]
  E1_asian_L <- y[10]
  I1_asian_L <- y[11]
  R1_asian_L <- y[12]
  S2_asian_L <- y[13]
  E2_asian_L <- y[14]
  I2_asian_L <- y[15]
  R2_asian_L <- y[16]
  
  D_asian <- y[17]
  Inc_asian <- y[18]
  RepCases_asian <- y[19]
  
  ## Black
  S1_black <- y[20]
  E1_black <- y[21]
  I1_black <- y[22]
  R1_black <- y[23]
  S2_black <- y[24]
  E2_black <- y[25]
  I2_black <- y[26]
  R2_black <- y[27]
  
  S1_black_L <- y[28]
  E1_black_L <- y[29]
  I1_black_L <- y[30]
  R1_black_L <- y[31]
  S2_black_L <- y[32]
  E2_black_L <- y[33]
  I2_black_L <- y[34]
  R2_black_L <- y[35]
  
  D_black <- y[36]
  Inc_black <- y[37]
  RepCases_black <- y[38]
  
  ## Other
  S1_other <- y[39]
  E1_other <- y[40]
  I1_other <- y[41]
  R1_other <- y[42]
  S2_other <- y[43]
  E2_other <- y[44]
  I2_other <- y[45]
  R2_other <- y[46]
  
  S1_other_L <- y[47]
  E1_other_L <- y[48]
  I1_other_L <- y[49]
  R1_other_L <- y[50]
  S2_other_L <- y[51]
  E2_other_L <- y[52]
  I2_other_L <- y[53]
  R2_other_L <- y[54]
  
  D_other <- y[55]
  Inc_other <- y[56]
  RepCases_other <- y[57]
  
  ## White
  S1_white <- y[58]
  E1_white <- y[59]
  I1_white <- y[60]
  R1_white <- y[61]
  S2_white <- y[62]
  E2_white <- y[63]
  I2_white <- y[64]
  R2_white <- y[65]
  
  S1_white_L <- y[66]
  E1_white_L <- y[67]
  I1_white_L <- y[68]
  R1_white_L <- y[69]
  S2_white_L <- y[70]
  E2_white_L <- y[71]
  I2_white_L <- y[72]
  R2_white_L <- y[73]
  
  D_white <- y[74]
  Inc_white <- y[75]
  RepCases_white <- y[76]
  
  ## calculations
  I_asian_eff <- I1_asian + I2_asian + (1-npi_efficacy)*(I1_asian_L + I2_asian_L)
  I_black_eff <- I1_black + I2_black + (1-npi_efficacy)*(I1_black_L + I2_black_L)
  I_other_eff <- I1_other + I2_other + (1-npi_efficacy)*(I1_other_L + I2_other_L)
  I_white_eff <- I1_white + I2_white + (1-npi_efficacy)*(I1_white_L + I2_white_L)
  
  P_asian_eff <- S1_asian + E1_asian + I1_asian + R1_asian + S2_asian + E2_asian + I2_asian + R2_asian  
  P_black_eff <- S1_black + E1_black + I1_black + R1_black + S2_black + E2_black + I2_black + R2_black 
  P_other_eff <- S1_other + E1_other + I1_other + R1_other + S2_other + E2_other + I2_other + R2_other 
  P_white_eff <- S1_white + E1_white + I1_white + R1_white + S2_white + E2_white + I2_white + R2_white 
  
  P_asian_total <- P_asian_eff + S1_asian_L + E1_asian_L + I1_asian_L + R1_asian_L + S2_asian_L + E2_asian_L + I2_asian_L + R2_asian_L
  P_black_total <- P_black_eff + S1_black_L + E1_black_L + I1_black_L + R1_black_L + S2_black_L + E2_black_L + I2_black_L + R2_black_L  
  P_other_total <- P_other_eff + S1_other_L + E1_other_L + I1_other_L + R1_other_L + S2_other_L + E2_other_L + I2_other_L + R2_other_L  
  P_white_total <- P_white_eff + S1_white_L + E1_white_L + I1_white_L + R1_white_L + S2_white_L + E2_white_L + I2_white_L + R2_white_L
  
  ## CHANGING HERE TO P total
  prop_asian_eff <- I_asian_eff/P_asian_total
  prop_black_eff <- I_black_eff/P_black_total
  prop_other_eff <- I_other_eff/P_other_total
  prop_white_eff <- I_white_eff/P_white_total
  
  asian_c_asian <- aca(t)
  asian_c_black <- acb(t)
  asian_c_other <- aco(t)
  asian_c_white <- acw(t)
  
  black_c_asian <- bca(t)
  black_c_black <- bcb(t)
  black_c_other <- bco(t)
  black_c_white <- bcw(t)
  
  other_c_asian <- oca(t)
  other_c_black <- ocb(t)
  other_c_other <- oco(t)
  other_c_white <- ocw(t)
  
  white_c_asian <- wca(t)
  white_c_black <- wcb(t)
  white_c_other <- wco(t)
  white_c_white <- wcw(t)
  
  # The vaccination rates are based on data referenced to the whole population. In this model, only susceptible and recovered can get vaccinated (S1, R1, S2, R2). 
  # Thus, the pool to be vaccinated needs to be adjusted. If v is the percent daily that is observed in the data, then
  # v_true*P = v_adj * (S1 + R1 + S2 + R2)
  # so that v_adj = v_true*P/(S1 + R1 + S2 + R2)
  
  v_scalar_asian <- P_asian_total/(S1_asian + S1_asian_L + R1_asian + R1_asian_L + S2_asian + S2_asian_L + R2_asian  + R2_asian_L)
  v_scalar_black <- P_black_total/(S1_black + S1_black_L + R1_black + R1_black_L + S2_black + S2_black_L + R2_black  + R2_black_L)
  v_scalar_other <- P_other_total/(S1_other + S1_other_L + R1_other + R1_other_L + S2_other + S2_other_L + R2_other  + R2_other_L)
  v_scalar_white <- P_white_total/(S1_white + S1_white_L + R1_white + R1_white_L + S2_white + S2_white_L + R2_white  + R2_white_L)
  
  nu_asian <- min(0.85*v_scalar_asian*vax_a(t),1)
  nu_black <- min(0.85*v_scalar_black*vax_b(t), 1)
  nu_other <- min(0.85*v_scalar_other*vax_o(t), 1)
  nu_white <- min(0.85*v_scalar_white*vax_w(t), 1)
  
  var <- variant(t)
  
  lambda1_asian <- var*mu*(aca(t)*prop_asian_eff + acb(t)*prop_black_eff  + aco(t)*prop_other_eff + acw(t)*prop_white_eff)
  lambda2_asian <- ir_2*lambda1_asian
  
  lambda1_black <- var*mu*(bca(t)*prop_asian_eff+ bcb(t)*prop_black_eff + bco(t)*prop_other_eff + bcw(t)*prop_white_eff)
  lambda2_black<- ir_2*lambda1_black
  
  lambda1_other <- var*mu*(oca(t)*prop_asian_eff+ ocb(t)*prop_black_eff + oco(t)*prop_other_eff + ocw(t)*prop_white_eff)
  lambda2_other <- ir_2*lambda1_other
  
  lambda1_white <- var*mu*(wca(t)*prop_asian_eff+ wcb(t)*prop_black_eff + wco(t)*prop_other_eff + wcw(t)*prop_white_eff)
  lambda2_white <- ir_2*lambda1_white
  
  ## Updates
  dS1_asian <- - lambda1_asian*S1_asian + npi_out_asian*S1_asian_L - nu_asian*S1_asian
  dE1_asian <- + lambda1_asian*S1_asian - epsilon*E1_asian + npi_out_asian*E1_asian_L
  dI1_asian <- + epsilon*E1_asian - gamma*I1_asian + npi_out_asian*I1_asian_L
  dR1_asian <- + (1-sigma1_asian)*gamma*I1_asian - omega*R1_asian + npi_out_asian*R1_asian_L - nu_asian*R1_asian
  dS2_asian <- + omega*R1_asian - lambda2_asian*S2_asian + omega*R2_asian + npi_out_asian*S2_asian_L - nu_asian*S2_asian
  dE2_asian <- + lambda2_asian*S2_asian - epsilon*E2_asian + npi_out_asian*E2_asian_L
  dI2_asian <- + epsilon*E2_asian - gamma*I2_asian + npi_out_asian*I2_asian_L
  dR2_asian <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian - omega*R2_asian + npi_out_asian*R2_asian_L + nu_asian*(S1_asian + R1_asian + S2_asian)
  
  dS1_asian_L <-  - npi_out_asian*S1_asian_L - nu_asian*S1_asian_L - (1-npi_efficacy)*lambda1_asian*S1_asian_L
  dE1_asian_L <-  - epsilon*E1_asian_L - npi_out_asian*E1_asian_L + (1-npi_efficacy)*lambda1_asian*S1_asian_L
  dI1_asian_L <- + epsilon*E1_asian_L - gamma*I1_asian_L - npi_out_asian*I1_asian_L
  dR1_asian_L <- + (1-sigma1_asian)*gamma*I1_asian_L - omega*R1_asian_L- npi_out_asian*R1_asian_L - nu_asian*R1_asian_L
  dS2_asian_L <- + omega*R1_asian_L + omega*R2_asian_L- npi_out_asian*S2_asian_L - nu_asian*S2_asian_L - (1-npi_efficacy)*lambda2_asian*S2_asian
  dE2_asian_L <- - epsilon*E2_asian_L - npi_out_asian*E2_asian_L + (1-npi_efficacy)*lambda2_asian*S2_asian
  dI2_asian_L <- + epsilon*E2_asian_L - gamma*I2_asian_L - npi_out_asian*I2_asian_L
  dR2_asian_L <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian_L - omega*R2_asian_L- npi_out_asian*R2_asian_L + nu_asian*(S1_asian_L + R1_asian_L + S2_asian_L)
  
  dD_asian <- sigma1_asian*gamma*(I1_asian + I1_asian_L) + ifr_2*sigma1_asian*gamma*(I2_asian + I2_asian_L)
  dInc_asian <- + epsilon*(E1_asian + E2_asian + E1_asian_L + E2_asian_L)
  dRepCases_asian <- rr_asian*epsilon*(E1_asian+ E2_asian + E1_asian_L + E2_asian_L)
  
  dS1_black <- - lambda1_black*S1_black + npi_out_black*S1_black_L - nu_black*S1_black
  dE1_black <- + lambda1_black*S1_black - epsilon*E1_black + npi_out_black*E1_black_L
  dI1_black <- + epsilon*E1_black - gamma*I1_black + npi_out_black*I1_black_L
  dR1_black <- + (1-sigma1_black)*gamma*I1_black - omega*R1_black + npi_out_black*R1_black_L - nu_black*R1_black
  dS2_black <- + omega*R1_black - lambda2_black*S2_black + omega*R2_black + npi_out_black*S2_black_L - nu_black*S2_black
  dE2_black <- + lambda2_black*S2_black - epsilon*E2_black + npi_out_black*E2_black_L
  dI2_black <- + epsilon*E2_black - gamma*I2_black + npi_out_black*I2_black_L
  dR2_black <- + (1-ifr_2*sigma1_black)*gamma*I2_black - omega*R2_black + npi_out_black*R2_black_L + nu_black*(S1_black + R1_black + S2_black)
  
  dS1_black_L <-  - npi_out_black*S1_black_L - nu_black*S1_black_L - (1-npi_efficacy)*lambda1_black*S1_black_L
  dE1_black_L <-  - epsilon*E1_black_L - npi_out_black*E1_black_L + (1-npi_efficacy)*lambda1_black*S1_black_L
  dI1_black_L <- + epsilon*E1_black_L - gamma*I1_black_L - npi_out_black*I1_black_L
  dR1_black_L <- + (1-sigma1_black)*gamma*I1_black_L - omega*R1_black_L- npi_out_black*R1_black_L - nu_black*R1_black_L
  dS2_black_L <- + omega*R1_black_L + omega*R2_black_L- npi_out_black*S2_black_L - nu_black*S2_black_L - (1-npi_efficacy)*lambda2_black*S2_black
  dE2_black_L <- - epsilon*E2_black_L - npi_out_black*E2_black_L + (1-npi_efficacy)*lambda2_black*S2_black
  dI2_black_L <- + epsilon*E2_black_L - gamma*I2_black_L - npi_out_black*I2_black_L
  dR2_black_L <- + (1-ifr_2*sigma1_black)*gamma*I2_black_L - omega*R2_black_L- npi_out_black*R2_black_L + nu_black*(S1_black_L + R1_black_L + S2_black_L)
  
  dD_black <- sigma1_black*gamma*(I1_black + I1_black_L) + ifr_2*sigma1_black*gamma*(I2_black + I2_black_L)
  dInc_black <- + epsilon*(E1_black + E2_black + E1_black_L + E2_black_L)
  dRepCases_black <- rr_black*epsilon*(E1_black+ E2_black + E1_black_L + E2_black_L)
  
  dS1_other <- - lambda1_other*S1_other + npi_out_other*S1_other_L - nu_other*S1_other
  dE1_other <- + lambda1_other*S1_other - epsilon*E1_other + npi_out_other*E1_other_L
  dI1_other <- + epsilon*E1_other - gamma*I1_other + npi_out_other*I1_other_L
  dR1_other <- + (1-sigma1_other)*gamma*I1_other - omega*R1_other + npi_out_other*R1_other_L - nu_other*R1_other
  dS2_other <- + omega*R1_other - lambda2_other*S2_other + omega*R2_other + npi_out_other*S2_other_L - nu_other*S2_other
  dE2_other <- + lambda2_other*S2_other - epsilon*E2_other + npi_out_other*E2_other_L
  dI2_other <- + epsilon*E2_other - gamma*I2_other + npi_out_other*I2_other_L
  dR2_other <- + (1-ifr_2*sigma1_other)*gamma*I2_other - omega*R2_other + npi_out_other*R2_other_L + nu_other*(S1_other + R1_other + S2_other)
  
  dS1_other_L <-  - npi_out_other*S1_other_L - nu_other*S1_other_L - (1-npi_efficacy)*lambda1_other*S1_other_L
  dE1_other_L <-  - epsilon*E1_other_L - npi_out_other*E1_other_L + (1-npi_efficacy)*lambda1_other*S1_other_L
  dI1_other_L <- + epsilon*E1_other_L - gamma*I1_other_L - npi_out_other*I1_other_L
  dR1_other_L <- + (1-sigma1_other)*gamma*I1_other_L - omega*R1_other_L- npi_out_other*R1_other_L - nu_other*R1_other_L
  dS2_other_L <- + omega*R1_other_L + omega*R2_other_L- npi_out_other*S2_other_L - nu_other*S2_other_L - (1-npi_efficacy)*lambda2_other*S2_other
  dE2_other_L <- - epsilon*E2_other_L - npi_out_other*E2_other_L + (1-npi_efficacy)*lambda2_other*S2_other
  dI2_other_L <- + epsilon*E2_other_L - gamma*I2_other_L - npi_out_other*I2_other_L
  dR2_other_L <- + (1-ifr_2*sigma1_other)*gamma*I2_other_L - omega*R2_other_L- npi_out_other*R2_other_L + nu_other*(S1_other_L + R1_other_L + S2_other_L)
  
  dD_other <- sigma1_other*gamma*(I1_other + I1_other_L) + ifr_2*sigma1_other*gamma*(I2_other + I2_other_L)
  dInc_other <- + epsilon*(E1_other + E2_other + E1_other_L + E2_other_L)
  dRepCases_other <- rr_other*epsilon*(E1_other+ E2_other + E1_other_L + E2_other_L)
  
  dS1_white <- - lambda1_white*S1_white + npi_out_white*S1_white_L - nu_white*S1_white
  dE1_white <- + lambda1_white*S1_white - epsilon*E1_white + npi_out_white*E1_white_L
  dI1_white <- + epsilon*E1_white - gamma*I1_white + npi_out_white*I1_white_L
  dR1_white <- + (1-sigma1_white)*gamma*I1_white - omega*R1_white + npi_out_white*R1_white_L - nu_white*R1_white
  dS2_white <- + omega*R1_white - lambda2_white*S2_white + omega*R2_white + npi_out_white*S2_white_L - nu_white*S2_white
  dE2_white <- + lambda2_white*S2_white - epsilon*E2_white + npi_out_white*E2_white_L
  dI2_white <- + epsilon*E2_white - gamma*I2_white + npi_out_white*I2_white_L
  dR2_white <- + (1-ifr_2*sigma1_white)*gamma*I2_white - omega*R2_white + npi_out_white*R2_white_L + nu_white*(S1_white + R1_white + S2_white)
  
  dS1_white_L <-  - npi_out_white*S1_white_L - nu_white*S1_white_L - (1-npi_efficacy)*lambda1_white*S1_white_L
  dE1_white_L <-  - epsilon*E1_white_L - npi_out_white*E1_white_L + (1-npi_efficacy)*lambda1_white*S1_white_L
  dI1_white_L <- + epsilon*E1_white_L - gamma*I1_white_L - npi_out_white*I1_white_L
  dR1_white_L <- + (1-sigma1_white)*gamma*I1_white_L - omega*R1_white_L- npi_out_white*R1_white_L - nu_white*R1_white_L
  dS2_white_L <- + omega*R1_white_L + omega*R2_white_L- npi_out_white*S2_white_L - nu_white*S2_white_L - (1-npi_efficacy)*lambda2_white*S2_white
  dE2_white_L <- - epsilon*E2_white_L - npi_out_white*E2_white_L + (1-npi_efficacy)*lambda2_white*S2_white
  dI2_white_L <- + epsilon*E2_white_L - gamma*I2_white_L - npi_out_white*I2_white_L
  dR2_white_L <- + (1-ifr_2*sigma1_white)*gamma*I2_white_L - omega*R2_white_L- npi_out_white*R2_white_L + nu_white*(S1_white_L + R1_white_L + S2_white_L)
  
  dD_white <- sigma1_white*gamma*(I1_white + I1_white_L) + ifr_2*sigma1_white*gamma*(I2_white + I2_white_L)
  dInc_white <- + epsilon*(E1_white + E2_white + E1_white_L + E2_white_L)
  dRepCases_white <- rr_white*epsilon*(E1_white+ E2_white + E1_white_L + E2_white_L)
  
  # Return list of gradients
  list(c(dS1_asian,dE1_asian,dI1_asian,dR1_asian,
         dS2_asian,dE2_asian,dI2_asian,dR2_asian,
         
         dS1_asian_L,dE1_asian_L,dI1_asian_L,dR1_asian_L,
         dS2_asian_L,dE2_asian_L,dI2_asian_L,dR2_asian_L,
         
         dD_asian, dInc_asian, dRepCases_asian,
         
         dS1_black,dE1_black,dI1_black,dR1_black,
         dS2_black,dE2_black,dI2_black,dR2_black,
         
         dS1_black_L,dE1_black_L,dI1_black_L,dR1_black_L,
         dS2_black_L,dE2_black_L,dI2_black_L,dR2_black_L,
         
         dD_black, dInc_black, dRepCases_black,
         
         dS1_other,dE1_other,dI1_other,dR1_other,
         dS2_other,dE2_other,dI2_other,dR2_other,
         
         dS1_other_L,dE1_other_L,dI1_other_L,dR1_other_L,
         dS2_other_L,dE2_other_L,dI2_other_L,dR2_other_L,
         
         dD_other, dInc_other, dRepCases_other,
         
         dS1_white,dE1_white,dI1_white,dR1_white,
         dS2_white,dE2_white,dI2_white,dR2_white,
         
         dS1_white_L,dE1_white_L,dI1_white_L,dR1_white_L,
         dS2_white_L,dE2_white_L,dI2_white_L,dR2_white_L,
         
         dD_white, dInc_white, dRepCases_white))
}



#### SEIR2 with vaccines ####
# For california
seir2_vax2 <- function(t, y, pars){
  # Parameters (fixed)
  epsilon <- pars["epsilon"] ## latent period
  omega <- pars["omega"] ## waning of immunity
  gamma <- pars["gamma"] ## infectious period
  ir_2 <- pars["ir_2"] ## scale infection risk for prior immunity
  ifr_2 <- pars["ifr_2"] ## scale infection fatality rate for prior immunity
  npi_efficacy <- pars["npi_efficacy"]
  
  ## Parameters (fit)
  mu <- pars["mu"] ## risk of infection
  
  rr_asian <- pars["rr_asian"] # reporting rate
  rr_black <- pars["rr_black"]
  rr_latino <- pars["rr_latino"]
  rr_other <- pars["rr_other"]
  rr_white <- pars["rr_white"]
  
  npi_out_asian <- pars["npi_out_asian"]
  npi_out_black <- pars["npi_out_black"]
  npi_out_latino <- pars["npi_out_latino"]
  npi_out_other <- pars["npi_out_other"]
  npi_out_white <- pars["npi_out_white"]
  
  sigma1_asian <- pars["sigma1_asian"] ## infection fatality rate
  sigma1_black <- pars["sigma1_black"] 
  sigma1_latino <- pars["sigma1_latino"] 
  sigma1_other <- pars["sigma1_other"] 
  sigma1_white <- pars["sigma1_white"]
  
  ## add for npi function
  npi_asian <- pars["npi_asian"]
  npi_black <- pars["npi_black"]
  npi_latino <- pars["npi_latino"]
  npi_other <- pars["npi_other"]
  npi_black <- pars["npi_black"]
  
  npi_secondwave <- pars["npi_secondwave"]
  npi_halfway <- pars["npi_halfway"]
  
  # State variables
  ## Asian
  S1_asian <- y[1]
  E1_asian <- y[2]
  I1_asian <- y[3]
  R1_asian <- y[4]
  S2_asian <- y[5]
  E2_asian <- y[6]
  I2_asian <- y[7]
  R2_asian <- y[8]
  
  S1_asian_L <- y[9]
  E1_asian_L <- y[10]
  I1_asian_L <- y[11]
  R1_asian_L <- y[12]
  S2_asian_L <- y[13]
  E2_asian_L <- y[14]
  I2_asian_L <- y[15]
  R2_asian_L <- y[16]
  
  D_asian <- y[17]
  Inc_asian <- y[18]
  RepCases_asian <- y[19]
  
  ## Black
  S1_black <- y[20]
  E1_black <- y[21]
  I1_black <- y[22]
  R1_black <- y[23]
  S2_black <- y[24]
  E2_black <- y[25]
  I2_black <- y[26]
  R2_black <- y[27]
  
  S1_black_L <- y[28]
  E1_black_L <- y[29]
  I1_black_L <- y[30]
  R1_black_L <- y[31]
  S2_black_L <- y[32]
  E2_black_L <- y[33]
  I2_black_L <- y[34]
  R2_black_L <- y[35]
  
  D_black <- y[36]
  Inc_black <- y[37]
  RepCases_black <- y[38]
  
  ## Latino
  S1_latino <- y[39]
  E1_latino <- y[40]
  I1_latino <- y[41]
  R1_latino <- y[42]
  S2_latino <- y[43]
  E2_latino <- y[44]
  I2_latino <- y[45]
  R2_latino <- y[46]
  
  S1_latino_L <- y[47]
  E1_latino_L <- y[48]
  I1_latino_L <- y[49]
  R1_latino_L <- y[50]
  S2_latino_L <- y[51]
  E2_latino_L <- y[52]
  I2_latino_L <- y[53]
  R2_latino_L <- y[54]
  
  D_latino <- y[55]
  Inc_latino <- y[56]
  RepCases_latino <- y[57]
  
  ## Other
  S1_other <- y[58]
  E1_other <- y[59]
  I1_other <- y[60]
  R1_other <- y[61]
  S2_other <- y[62]
  E2_other <- y[63]
  I2_other <- y[64]
  R2_other <- y[65]
  
  S1_other_L <- y[66]
  E1_other_L <- y[67]
  I1_other_L <- y[68]
  R1_other_L <- y[69]
  S2_other_L <- y[70]
  E2_other_L <- y[71]
  I2_other_L <- y[72]
  R2_other_L <- y[73]
  
  D_other <- y[74]
  Inc_other <- y[75]
  RepCases_other <- y[76]
  
  ## White
  S1_white <- y[77]
  E1_white <- y[78]
  I1_white <- y[79]
  R1_white <- y[80]
  S2_white <- y[81]
  E2_white <- y[82]
  I2_white <- y[83]
  R2_white <- y[84]
  
  S1_white_L <- y[85]
  E1_white_L <- y[86]
  I1_white_L <- y[87]
  R1_white_L <- y[88]
  S2_white_L <- y[89]
  E2_white_L <- y[90]
  I2_white_L <- y[91]
  R2_white_L <- y[92]
  
  D_white <- y[93]
  Inc_white <- y[94]
  RepCases_white <- y[95]
  
  ## calculations
  I_asian_eff <- I1_asian + I2_asian + (1-npi_efficacy)*(I1_asian_L + I2_asian_L)
  I_black_eff <- I1_black + I2_black + (1-npi_efficacy)*(I1_black_L + I2_black_L)
  I_latino_eff <- I1_latino + I1_latino + (1-npi_efficacy)*(I1_latino_L + I2_latino_L)
  I_other_eff <- I1_other + I2_other + (1-npi_efficacy)*(I1_other_L + I2_other_L)
  I_white_eff <- I1_white + I2_white + (1-npi_efficacy)*(I1_white_L + I2_white_L)
  
  P_asian_eff <- S1_asian + E1_asian + I1_asian + R1_asian + S2_asian + E2_asian + I2_asian + R2_asian  
  P_black_eff <- S1_black + E1_black + I1_black + R1_black + S2_black + E2_black + I2_black + R2_black 
  P_latino_eff <- S1_latino + E1_latino + I1_latino + R1_latino + S2_latino + E2_latino + I2_latino + R2_latino 
  P_other_eff <- S1_other + E1_other + I1_other + R1_other + S2_other + E2_other + I2_other + R2_other 
  P_white_eff <- S1_white + E1_white + I1_white + R1_white + S2_white + E2_white + I2_white + R2_white 
  
  P_asian_total <- P_asian_eff + S1_asian_L + E1_asian_L + I1_asian_L + R1_asian_L + S2_asian_L + E2_asian_L + I2_asian_L + R2_asian_L
  P_black_total <- P_black_eff + S1_black_L + E1_black_L + I1_black_L + R1_black_L + S2_black_L + E2_black_L + I2_black_L + R2_black_L  
  P_latino_total <- P_latino_eff + S1_latino_L + E1_latino_L + I1_latino_L + R1_latino_L + S2_latino_L + E2_latino_L + I2_latino_L + R2_latino_L  
  P_other_total <- P_other_eff + S1_other_L + E1_other_L + I1_other_L + R1_other_L + S2_other_L + E2_other_L + I2_other_L + R2_other_L  
  P_white_total <- P_white_eff + S1_white_L + E1_white_L + I1_white_L + R1_white_L + S2_white_L + E2_white_L + I2_white_L + R2_white_L
  
  ## CHANGING HERE TO P total
  prop_asian_eff <- I_asian_eff/P_asian_total
  prop_black_eff <- I_black_eff/P_black_total
  prop_latino_eff <- I_latino_eff/P_latino_total
  prop_other_eff <- I_other_eff/P_other_total
  prop_white_eff <- I_white_eff/P_white_total
  
  ## contacts 
  asian_c_asian <- aca(t)
  asian_c_black <- acb(t)
  asian_c_latino <- acl(t)
  asian_c_other <- aco(t)
  asian_c_white <- acw(t)
  
  black_c_asian <- bca(t)
  black_c_black <- bcb(t)
  black_c_latino <- bcl(t)
  black_c_other <- bco(t)
  black_c_white <- bcw(t)
  
  latino_c_asian <- lca(t)
  latino_c_black <- lcb(t)
  latino_c_latino <- lcl(t)
  latino_c_other <- lco(t)
  latino_c_white <- lcw(t)
  
  other_c_asian <- oca(t)
  other_c_black <- ocb(t)
  other_c_latino <- ocl(t)
  other_c_other <- oco(t)
  other_c_white <- ocw(t)
  
  white_c_asian <- wca(t)
  white_c_black <- wcb(t)
  white_c_latino <- wcl(t)
  white_c_other <- wco(t)
  white_c_white <- wcw(t)
  
  ## npi wave 2
  ## The max npi compliance is referenced to the whole population. Given that only non-L move to L, 
  ## if n is the rate moving into npis, then 
  # n_true * P = n_adj * (S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2)
  # n_adj = n_true * P/P_eff
  # maxing at 1
  
  npi_scalar_asian <- P_asian_total/P_asian_eff
  npi_scalar_black <- P_black_total/P_black_eff
  npi_scalar_latino <- P_latino_total/P_latino_eff
  npi_scalar_other <- P_other_total/P_other_eff
  npi_scalar_white <- P_white_total/P_white_eff
  
  npi_a <- npi_weekly(npi_secondwave*npi_asian, npi_halfway, 2)
  npi_b <- npi_weekly(npi_secondwave*npi_black, npi_halfway, 2)
  npi_l <- npi_weekly(npi_secondwave*npi_latino, npi_halfway, 2)
  npi_o <- npi_weekly(npi_secondwave*npi_other, npi_halfway, 2)
  npi_w <- npi_weekly(npi_secondwave*npi_white, npi_halfway, 2)
  
  npi_asian_rate <- min(npi_scalar_asian*npi_a(t-315), 1)
  npi_black_rate <- min(npi_scalar_black*npi_b(t-315), 1)
  npi_latino_rate <- min(npi_scalar_latino*npi_l(t-315), 1)
  npi_other_rate <- min(npi_scalar_other*npi_o(t-315), 1)
  npi_white_rate <- min(npi_scalar_white*npi_w(t-315), 1)
  
  # The vaccination rates are based on data referenced to the whole population. In this model, only susceptible and recovered can get vaccinated (S1, R1, S2, R2). 
  # Thus, the pool to be vaccinated needs to be adjusted. If v is the percent daily that is observed in the data, then
  # v*P = v_adj * (S1 + R1 + S2 + R2)
  # so that v_adj = v*P/(S1 + R1 + S2 + R2)
  
  v_scalar_asian <- P_asian_total/(S1_asian + S1_asian_L + R1_asian + R1_asian_L + S2_asian + S2_asian_L + R2_asian  + R2_asian_L)
  v_scalar_black <- P_black_total/(S1_black + S1_black_L + R1_black + R1_black_L + S2_black + S2_black_L + R2_black  + R2_black_L)
  v_scalar_latino <- P_latino_total/(S1_latino + S1_latino_L + R1_latino + R1_latino_L + S2_latino + S2_latino_L + R2_latino  + R2_latino_L)
  v_scalar_other <- P_other_total/(S1_other + S1_other_L + R1_other + R1_other_L + S2_other + S2_other_L + R2_other  + R2_other_L)
  v_scalar_white <- P_white_total/(S1_white + S1_white_L + R1_white + R1_white_L + S2_white + S2_white_L + R2_white  + R2_white_L)
  
  nu_asian <- min(0.85*v_scalar_asian*vax_a(t),1)
  nu_black <- min(0.85*v_scalar_black*vax_b(t), 1)
  nu_latino <- min(0.85*v_scalar_latino*vax_l(t), 1)
  nu_other <- min(0.85*v_scalar_other*vax_o(t), 1)
  nu_white <- min(0.85*v_scalar_white*vax_w(t), 1)
  
  var <- variant(t)
  
  lambda1_asian <- var*mu*(aca(t)*prop_asian_eff + acb(t)*prop_black_eff + acl(t)*prop_latino_eff + aco(t)*prop_other_eff + acw(t)*prop_white_eff)
  lambda2_asian <- ir_2*lambda1_asian
  
  lambda1_black <- var*mu*(bca(t)*prop_asian_eff+ bcb(t)*prop_black_eff + bcl(t)*prop_latino_eff + bco(t)*prop_other_eff + bcw(t)*prop_white_eff)
  lambda2_black<- ir_2*lambda1_black
  
  lambda1_latino <- var*mu*(lca(t)*prop_asian_eff + lcb(t)*prop_black_eff + lcl(t)*prop_latino_eff + lco(t)*prop_other_eff + lcw(t)*prop_white_eff)
  lambda2_latino<- ir_2*lambda1_latino
  
  lambda1_other <- var*mu*(oca(t)*prop_asian_eff+ ocb(t)*prop_black_eff + ocl(t)*prop_latino_eff + oco(t)*prop_other_eff + ocw(t)*prop_white_eff)
  lambda2_other <- ir_2*lambda1_other
  
  lambda1_white <- var*mu*(wca(t)*prop_asian_eff+ wcb(t)*prop_black_eff + wcl(t)*prop_latino_eff + wco(t)*prop_other_eff + wcw(t)*prop_white_eff)
  lambda2_white <- ir_2*lambda1_white
  
  ## Updates
  dS1_asian <- - lambda1_asian*S1_asian + npi_out_asian*S1_asian_L - nu_asian*S1_asian - npi_asian_rate*S1_asian
  dE1_asian <- + lambda1_asian*S1_asian - epsilon*E1_asian  + npi_out_asian*E1_asian_L - npi_asian_rate*E1_asian
  dI1_asian <- + epsilon*E1_asian - gamma*I1_asian + npi_out_asian*I1_asian_L - npi_asian_rate*I1_asian
  dR1_asian <- + (1-sigma1_asian)*gamma*I1_asian - omega*R1_asian + npi_out_asian*R1_asian_L - nu_asian*R1_asian - npi_asian_rate*R1_asian
  dS2_asian <- + omega*R1_asian - lambda2_asian*S2_asian + omega*R2_asian + npi_out_asian*S2_asian_L - nu_asian*S2_asian - npi_asian_rate*S2_asian
  dE2_asian <- + lambda2_asian*S2_asian - epsilon*E2_asian + npi_out_asian*E2_asian_L - npi_asian_rate*E2_asian
  dI2_asian <- + epsilon*E2_asian - gamma*I2_asian + npi_out_asian*I2_asian_L - npi_asian_rate*I2_asian
  dR2_asian <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian - omega*R2_asian + npi_out_asian*R2_asian_L + nu_asian*(S1_asian + R1_asian + S2_asian) - npi_asian_rate*R2_asian
  
  dS1_asian_L <-  - npi_out_asian*S1_asian_L - nu_asian*S1_asian_L - (1-npi_efficacy)*lambda1_asian*S1_asian_L + npi_asian_rate*S1_asian
  dE1_asian_L <-  - epsilon*E1_asian_L  - npi_out_asian*E1_asian_L + (1-npi_efficacy)*lambda1_asian*S1_asian_L + npi_asian_rate*E1_asian
  dI1_asian_L <- + epsilon*E1_asian_L - gamma*I1_asian_L - npi_out_asian*I1_asian_L + npi_asian_rate*I1_asian
  dR1_asian_L <- + (1-sigma1_asian)*gamma*I1_asian_L - omega*R1_asian_L  - npi_out_asian*R1_asian_L - nu_asian*R1_asian_L + npi_asian_rate*R1_asian
  dS2_asian_L <- + omega*R1_asian_L + omega*R2_asian_L - npi_out_asian*S2_asian_L - nu_asian*S2_asian_L - (1-npi_efficacy)*lambda2_asian*S2_asian_L + npi_asian_rate*S2_asian
  dE2_asian_L <- - epsilon*E2_asian_L - npi_out_asian*E2_asian_L + (1-npi_efficacy)*lambda2_asian*S2_asian_L + npi_asian_rate*E2_asian
  dI2_asian_L <- + epsilon*E2_asian_L - gamma*I2_asian_L  - npi_out_asian*I2_asian_L + npi_asian_rate*I2_asian
  dR2_asian_L <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian_L - omega*R2_asian_L - npi_out_asian*R2_asian_L + nu_asian*(S1_asian_L + R1_asian_L + S2_asian_L) + npi_asian_rate*R2_asian
  
  dD_asian <- sigma1_asian*gamma*(I1_asian + I1_asian_L) + ifr_2*sigma1_asian*gamma*(I2_asian + I2_asian_L)
  dInc_asian <- + epsilon*(E1_asian + E2_asian + E1_asian_L + E2_asian_L)
  dRepCases_asian <- rr_asian*epsilon*(E1_asian+ E2_asian + E1_asian_L + E2_asian_L)
  
  dS1_black <- - lambda1_black*S1_black + npi_out_black*S1_black_L - nu_black*S1_black - npi_black_rate*S1_black
  dE1_black <- + lambda1_black*S1_black - epsilon*E1_black  + npi_out_black*E1_black_L - npi_black_rate*E1_black
  dI1_black <- + epsilon*E1_black - gamma*I1_black + npi_out_black*I1_black_L - npi_black_rate*I1_black
  dR1_black <- + (1-sigma1_black)*gamma*I1_black - omega*R1_black + npi_out_black*R1_black_L - nu_black*R1_black - npi_black_rate*R1_black
  dS2_black <- + omega*R1_black - lambda2_black*S2_black + omega*R2_black + npi_out_black*S2_black_L - nu_black*S2_black - npi_black_rate*S2_black
  dE2_black <- + lambda2_black*S2_black - epsilon*E2_black + npi_out_black*E2_black_L - npi_black_rate*E2_black
  dI2_black <- + epsilon*E2_black - gamma*I2_black + npi_out_black*I2_black_L - npi_black_rate*I2_black
  dR2_black <- + (1-ifr_2*sigma1_black)*gamma*I2_black - omega*R2_black + npi_out_black*R2_black_L + nu_black*(S1_black + R1_black + S2_black) - npi_black_rate*R2_black
  
  dS1_black_L <-  - npi_out_black*S1_black_L - nu_black*S1_black_L - (1-npi_efficacy)*lambda1_black*S1_black_L + npi_black_rate*S1_black
  dE1_black_L <-  - epsilon*E1_black_L  - npi_out_black*E1_black_L + (1-npi_efficacy)*lambda1_black*S1_black_L + npi_black_rate*E1_black
  dI1_black_L <- + epsilon*E1_black_L - gamma*I1_black_L - npi_out_black*I1_black_L + npi_black_rate*I1_black
  dR1_black_L <- + (1-sigma1_black)*gamma*I1_black_L - omega*R1_black_L  - npi_out_black*R1_black_L - nu_black*R1_black_L + npi_black_rate*R1_black
  dS2_black_L <- + omega*R1_black_L + omega*R2_black_L - npi_out_black*S2_black_L - nu_black*S2_black_L - (1-npi_efficacy)*lambda2_black*S2_black_L + npi_black_rate*S2_black
  dE2_black_L <- - epsilon*E2_black_L - npi_out_black*E2_black_L + (1-npi_efficacy)*lambda2_black*S2_black_L + npi_black_rate*E2_black
  dI2_black_L <- + epsilon*E2_black_L - gamma*I2_black_L  - npi_out_black*I2_black_L + npi_black_rate*I2_black
  dR2_black_L <- + (1-ifr_2*sigma1_black)*gamma*I2_black_L - omega*R2_black_L - npi_out_black*R2_black_L + nu_black*(S1_black_L + R1_black_L + S2_black_L) + npi_black_rate*R2_black
  
  dD_black <- sigma1_black*gamma*(I1_black + I1_black_L) + ifr_2*sigma1_black*gamma*(I2_black + I2_black_L)
  dInc_black <- + epsilon*(E1_black + E2_black + E1_black_L + E2_black_L)
  dRepCases_black <- rr_black*epsilon*(E1_black+ E2_black + E1_black_L + E2_black_L)
  
  dS1_latino <- - lambda1_latino*S1_latino + npi_out_latino*S1_latino_L - nu_latino*S1_latino - npi_latino_rate*S1_latino
  dE1_latino <- + lambda1_latino*S1_latino - epsilon*E1_latino  + npi_out_latino*E1_latino_L - npi_latino_rate*E1_latino
  dI1_latino <- + epsilon*E1_latino - gamma*I1_latino + npi_out_latino*I1_latino_L - npi_latino_rate*I1_latino
  dR1_latino <- + (1-sigma1_latino)*gamma*I1_latino - omega*R1_latino + npi_out_latino*R1_latino_L - nu_latino*R1_latino - npi_latino_rate*R1_latino
  dS2_latino <- + omega*R1_latino - lambda2_latino*S2_latino + omega*R2_latino + npi_out_latino*S2_latino_L - nu_latino*S2_latino - npi_latino_rate*S2_latino
  dE2_latino <- + lambda2_latino*S2_latino - epsilon*E2_latino + npi_out_latino*E2_latino_L - npi_latino_rate*E2_latino
  dI2_latino <- + epsilon*E2_latino - gamma*I2_latino + npi_out_latino*I2_latino_L - npi_latino_rate*I2_latino
  dR2_latino <- + (1-ifr_2*sigma1_latino)*gamma*I2_latino - omega*R2_latino + npi_out_latino*R2_latino_L + nu_latino*(S1_latino + R1_latino + S2_latino) - npi_latino_rate*R2_latino
  
  dS1_latino_L <-  - npi_out_latino*S1_latino_L - nu_latino*S1_latino_L - (1-npi_efficacy)*lambda1_latino*S1_latino_L + npi_latino_rate*S1_latino
  dE1_latino_L <-  - epsilon*E1_latino_L  - npi_out_latino*E1_latino_L + (1-npi_efficacy)*lambda1_latino*S1_latino_L + npi_latino_rate*E1_latino
  dI1_latino_L <- + epsilon*E1_latino_L - gamma*I1_latino_L - npi_out_latino*I1_latino_L + npi_latino_rate*I1_latino
  dR1_latino_L <- + (1-sigma1_latino)*gamma*I1_latino_L - omega*R1_latino_L  - npi_out_latino*R1_latino_L - nu_latino*R1_latino_L + npi_latino_rate*R1_latino
  dS2_latino_L <- + omega*R1_latino_L + omega*R2_latino_L - npi_out_latino*S2_latino_L - nu_latino*S2_latino_L - (1-npi_efficacy)*lambda2_latino*S2_latino_L + npi_latino_rate*S2_latino
  dE2_latino_L <- - epsilon*E2_latino_L - npi_out_latino*E2_latino_L + (1-npi_efficacy)*lambda2_latino*S2_latino_L + npi_latino_rate*E2_latino
  dI2_latino_L <- + epsilon*E2_latino_L - gamma*I2_latino_L  - npi_out_latino*I2_latino_L + npi_latino_rate*I2_latino
  dR2_latino_L <- + (1-ifr_2*sigma1_latino)*gamma*I2_latino_L - omega*R2_latino_L - npi_out_latino*R2_latino_L + nu_latino*(S1_latino_L + R1_latino_L + S2_latino_L) + npi_latino_rate*R2_latino
  
  dD_latino <- sigma1_latino*gamma*(I1_latino + I1_latino_L) + ifr_2*sigma1_latino*gamma*(I2_latino + I2_latino_L)
  dInc_latino <- + epsilon*(E1_latino + E2_latino + E1_latino_L + E2_latino_L)
  dRepCases_latino <- rr_latino*epsilon*(E1_latino+ E2_latino + E1_latino_L + E2_latino_L)
  
  dS1_other <- - lambda1_other*S1_other + npi_out_other*S1_other_L - nu_other*S1_other - npi_other_rate*S1_other
  dE1_other <- + lambda1_other*S1_other - epsilon*E1_other  + npi_out_other*E1_other_L - npi_other_rate*E1_other
  dI1_other <- + epsilon*E1_other - gamma*I1_other + npi_out_other*I1_other_L - npi_other_rate*I1_other
  dR1_other <- + (1-sigma1_other)*gamma*I1_other - omega*R1_other + npi_out_other*R1_other_L - nu_other*R1_other - npi_other_rate*R1_other
  dS2_other <- + omega*R1_other - lambda2_other*S2_other + omega*R2_other + npi_out_other*S2_other_L - nu_other*S2_other - npi_other_rate*S2_other
  dE2_other <- + lambda2_other*S2_other - epsilon*E2_other + npi_out_other*E2_other_L - npi_other_rate*E2_other
  dI2_other <- + epsilon*E2_other - gamma*I2_other + npi_out_other*I2_other_L - npi_other_rate*I2_other
  dR2_other <- + (1-ifr_2*sigma1_other)*gamma*I2_other - omega*R2_other + npi_out_other*R2_other_L + nu_other*(S1_other + R1_other + S2_other) - npi_other_rate*R2_other
  
  dS1_other_L <-  - npi_out_other*S1_other_L - nu_other*S1_other_L - (1-npi_efficacy)*lambda1_other*S1_other_L + npi_other_rate*S1_other
  dE1_other_L <-  - epsilon*E1_other_L  - npi_out_other*E1_other_L + (1-npi_efficacy)*lambda1_other*S1_other_L + npi_other_rate*E1_other
  dI1_other_L <- + epsilon*E1_other_L - gamma*I1_other_L - npi_out_other*I1_other_L + npi_other_rate*I1_other
  dR1_other_L <- + (1-sigma1_other)*gamma*I1_other_L - omega*R1_other_L  - npi_out_other*R1_other_L - nu_other*R1_other_L + npi_other_rate*R1_other
  dS2_other_L <- + omega*R1_other_L + omega*R2_other_L - npi_out_other*S2_other_L - nu_other*S2_other_L - (1-npi_efficacy)*lambda2_other*S2_other_L + npi_other_rate*S2_other
  dE2_other_L <- - epsilon*E2_other_L - npi_out_other*E2_other_L + (1-npi_efficacy)*lambda2_other*S2_other_L + npi_other_rate*E2_other
  dI2_other_L <- + epsilon*E2_other_L - gamma*I2_other_L  - npi_out_other*I2_other_L + npi_other_rate*I2_other
  dR2_other_L <- + (1-ifr_2*sigma1_other)*gamma*I2_other_L - omega*R2_other_L - npi_out_other*R2_other_L + nu_other*(S1_other_L + R1_other_L + S2_other_L) + npi_other_rate*R2_other
  
  dD_other <- sigma1_other*gamma*(I1_other + I1_other_L) + ifr_2*sigma1_other*gamma*(I2_other + I2_other_L)
  dInc_other <- + epsilon*(E1_other + E2_other + E1_other_L + E2_other_L)
  dRepCases_other <- rr_other*epsilon*(E1_other+ E2_other + E1_other_L + E2_other_L)
  
  dS1_white <- - lambda1_white*S1_white + npi_out_white*S1_white_L - nu_white*S1_white - npi_white_rate*S1_white
  dE1_white <- + lambda1_white*S1_white - epsilon*E1_white  + npi_out_white*E1_white_L - npi_white_rate*E1_white
  dI1_white <- + epsilon*E1_white - gamma*I1_white + npi_out_white*I1_white_L - npi_white_rate*I1_white
  dR1_white <- + (1-sigma1_white)*gamma*I1_white - omega*R1_white + npi_out_white*R1_white_L - nu_white*R1_white - npi_white_rate*R1_white
  dS2_white <- + omega*R1_white - lambda2_white*S2_white + omega*R2_white + npi_out_white*S2_white_L - nu_white*S2_white - npi_white_rate*S2_white
  dE2_white <- + lambda2_white*S2_white - epsilon*E2_white + npi_out_white*E2_white_L - npi_white_rate*E2_white
  dI2_white <- + epsilon*E2_white - gamma*I2_white + npi_out_white*I2_white_L - npi_white_rate*I2_white
  dR2_white <- + (1-ifr_2*sigma1_white)*gamma*I2_white - omega*R2_white + npi_out_white*R2_white_L + nu_white*(S1_white + R1_white + S2_white) - npi_white_rate*R2_white
  
  dS1_white_L <-  - npi_out_white*S1_white_L - nu_white*S1_white_L - (1-npi_efficacy)*lambda1_white*S1_white_L + npi_white_rate*S1_white
  dE1_white_L <-  - epsilon*E1_white_L  - npi_out_white*E1_white_L + (1-npi_efficacy)*lambda1_white*S1_white_L + npi_white_rate*E1_white
  dI1_white_L <- + epsilon*E1_white_L - gamma*I1_white_L - npi_out_white*I1_white_L + npi_white_rate*I1_white
  dR1_white_L <- + (1-sigma1_white)*gamma*I1_white_L - omega*R1_white_L  - npi_out_white*R1_white_L - nu_white*R1_white_L + npi_white_rate*R1_white
  dS2_white_L <- + omega*R1_white_L + omega*R2_white_L - npi_out_white*S2_white_L - nu_white*S2_white_L - (1-npi_efficacy)*lambda2_white*S2_white_L + npi_white_rate*S2_white
  dE2_white_L <- - epsilon*E2_white_L - npi_out_white*E2_white_L + (1-npi_efficacy)*lambda2_white*S2_white_L + npi_white_rate*E2_white
  dI2_white_L <- + epsilon*E2_white_L - gamma*I2_white_L  - npi_out_white*I2_white_L + npi_white_rate*I2_white
  dR2_white_L <- + (1-ifr_2*sigma1_white)*gamma*I2_white_L - omega*R2_white_L - npi_out_white*R2_white_L + nu_white*(S1_white_L + R1_white_L + S2_white_L) + npi_white_rate*R2_white
  
  dD_white <- sigma1_white*gamma*(I1_white + I1_white_L) + ifr_2*sigma1_white*gamma*(I2_white + I2_white_L)
  dInc_white <- + epsilon*(E1_white + E2_white + E1_white_L + E2_white_L)
  dRepCases_white <- rr_white*epsilon*(E1_white+ E2_white + E1_white_L + E2_white_L)
  
  # Return list of gradients
  list(c(dS1_asian,dE1_asian,dI1_asian,dR1_asian,
         dS2_asian,dE2_asian,dI2_asian,dR2_asian,
         
         dS1_asian_L,dE1_asian_L,dI1_asian_L,dR1_asian_L,
         dS2_asian_L,dE2_asian_L,dI2_asian_L,dR2_asian_L,
         
         dD_asian, dInc_asian, dRepCases_asian,
         
         dS1_black,dE1_black,dI1_black,dR1_black,
         dS2_black,dE2_black,dI2_black,dR2_black,
         
         dS1_black_L,dE1_black_L,dI1_black_L,dR1_black_L,
         dS2_black_L,dE2_black_L,dI2_black_L,dR2_black_L,
         
         dD_black, dInc_black, dRepCases_black,
         
         dS1_latino,dE1_latino,dI1_latino,dR1_latino,
         dS2_latino,dE2_latino,dI2_latino,dR2_latino,
         
         dS1_latino_L,dE1_latino_L,dI1_latino_L,dR1_latino_L,
         dS2_latino_L,dE2_latino_L,dI2_latino_L,dR2_latino_L,
         
         dD_latino, dInc_latino, dRepCases_latino,
         
         dS1_other,dE1_other,dI1_other,dR1_other,
         dS2_other,dE2_other,dI2_other,dR2_other,
         
         dS1_other_L,dE1_other_L,dI1_other_L,dR1_other_L,
         dS2_other_L,dE2_other_L,dI2_other_L,dR2_other_L,
         
         dD_other, dInc_other, dRepCases_other,
         
         dS1_white,dE1_white,dI1_white,dR1_white,
         dS2_white,dE2_white,dI2_white,dR2_white,
         
         dS1_white_L,dE1_white_L,dI1_white_L,dR1_white_L,
         dS2_white_L,dE2_white_L,dI2_white_L,dR2_white_L,
         
         dD_white, dInc_white, dRepCases_white))
}


#### SEIR2 with vaccines ####
# for north carolina
seir2_nolatino_vax2 <- function(t, y, pars){
  # Parameters (fixed)
  epsilon <- pars["epsilon"] ## latent period
  omega <- pars["omega"] ## waning of immunity
  gamma <- pars["gamma"] ## infectious period
  ir_2 <- pars["ir_2"] ## scale infection risk for prior immunity
  ifr_2 <- pars["ifr_2"] ## scale infection fatality rate for prior immunity
  npi_efficacy <- pars["npi_efficacy"] ## npi efficacy (0-1)
  
  ## Parameters (fit)
  mu <- pars["mu"] ## risk of infection
  
  rr_asian <- pars["rr_asian"] # reporting rate
  rr_black <- pars["rr_black"]
  rr_other <- pars["rr_other"]
  rr_white <- pars["rr_white"]
  
  npi_out_asian <- pars["npi_out_asian"]
  npi_out_black <- pars["npi_out_black"]
  npi_out_other <- pars["npi_out_other"]
  npi_out_white <- pars["npi_out_white"]
  
  sigma1_asian <- pars["sigma1_asian"] ## infection fatality rate
  sigma1_black <- pars["sigma1_black"] 
  sigma1_other <- pars["sigma1_other"] 
  sigma1_white <- pars["sigma1_white"] 
  
  ## add for npi function
  npi_asian <- pars["npi_asian"]
  npi_black <- pars["npi_black"]
  npi_latino <- pars["npi_latino"]
  npi_other <- pars["npi_other"]
  npi_black <- pars["npi_black"]
  
  npi_secondwave <- pars["npi_secondwave"]
  npi_halfway <- pars["npi_halfway"]
  
  # State variables
  ## Asian
  S1_asian <- y[1]
  E1_asian <- y[2]
  I1_asian <- y[3]
  R1_asian <- y[4]
  S2_asian <- y[5]
  E2_asian <- y[6]
  I2_asian <- y[7]
  R2_asian <- y[8]
  
  S1_asian_L <- y[9]
  E1_asian_L <- y[10]
  I1_asian_L <- y[11]
  R1_asian_L <- y[12]
  S2_asian_L <- y[13]
  E2_asian_L <- y[14]
  I2_asian_L <- y[15]
  R2_asian_L <- y[16]
  
  D_asian <- y[17]
  Inc_asian <- y[18]
  RepCases_asian <- y[19]
  
  ## Black
  S1_black <- y[20]
  E1_black <- y[21]
  I1_black <- y[22]
  R1_black <- y[23]
  S2_black <- y[24]
  E2_black <- y[25]
  I2_black <- y[26]
  R2_black <- y[27]
  
  S1_black_L <- y[28]
  E1_black_L <- y[29]
  I1_black_L <- y[30]
  R1_black_L <- y[31]
  S2_black_L <- y[32]
  E2_black_L <- y[33]
  I2_black_L <- y[34]
  R2_black_L <- y[35]
  
  D_black <- y[36]
  Inc_black <- y[37]
  RepCases_black <- y[38]
  
  ## Other
  S1_other <- y[39]
  E1_other <- y[40]
  I1_other <- y[41]
  R1_other <- y[42]
  S2_other <- y[43]
  E2_other <- y[44]
  I2_other <- y[45]
  R2_other <- y[46]
  
  S1_other_L <- y[47]
  E1_other_L <- y[48]
  I1_other_L <- y[49]
  R1_other_L <- y[50]
  S2_other_L <- y[51]
  E2_other_L <- y[52]
  I2_other_L <- y[53]
  R2_other_L <- y[54]
  
  D_other <- y[55]
  Inc_other <- y[56]
  RepCases_other <- y[57]
  
  ## White
  S1_white <- y[58]
  E1_white <- y[59]
  I1_white <- y[60]
  R1_white <- y[61]
  S2_white <- y[62]
  E2_white <- y[63]
  I2_white <- y[64]
  R2_white <- y[65]
  
  S1_white_L <- y[66]
  E1_white_L <- y[67]
  I1_white_L <- y[68]
  R1_white_L <- y[69]
  S2_white_L <- y[70]
  E2_white_L <- y[71]
  I2_white_L <- y[72]
  R2_white_L <- y[73]
  
  D_white <- y[74]
  Inc_white <- y[75]
  RepCases_white <- y[76]
  
  ## calculations
  I_asian_eff <- I1_asian + I2_asian + (1-npi_efficacy)*(I1_asian_L + I2_asian_L)
  I_black_eff <- I1_black + I2_black + (1-npi_efficacy)*(I1_black_L + I2_black_L)
  I_other_eff <- I1_other + I2_other + (1-npi_efficacy)*(I1_other_L + I2_other_L)
  I_white_eff <- I1_white + I2_white + (1-npi_efficacy)*(I1_white_L + I2_white_L)
  
  P_asian_eff <- S1_asian + E1_asian + I1_asian + R1_asian + S2_asian + E2_asian + I2_asian + R2_asian  
  P_black_eff <- S1_black + E1_black + I1_black + R1_black + S2_black + E2_black + I2_black + R2_black 
  P_other_eff <- S1_other + E1_other + I1_other + R1_other + S2_other + E2_other + I2_other + R2_other 
  P_white_eff <- S1_white + E1_white + I1_white + R1_white + S2_white + E2_white + I2_white + R2_white 
  
  P_asian_total <- P_asian_eff + S1_asian_L + E1_asian_L + I1_asian_L + R1_asian_L + S2_asian_L + E2_asian_L + I2_asian_L + R2_asian_L
  P_black_total <- P_black_eff + S1_black_L + E1_black_L + I1_black_L + R1_black_L + S2_black_L + E2_black_L + I2_black_L + R2_black_L  
  P_other_total <- P_other_eff + S1_other_L + E1_other_L + I1_other_L + R1_other_L + S2_other_L + E2_other_L + I2_other_L + R2_other_L  
  P_white_total <- P_white_eff + S1_white_L + E1_white_L + I1_white_L + R1_white_L + S2_white_L + E2_white_L + I2_white_L + R2_white_L
  
  ## CHANGING HERE TO P total
  prop_asian_eff <- I_asian_eff/P_asian_total
  prop_black_eff <- I_black_eff/P_black_total
  prop_other_eff <- I_other_eff/P_other_total
  prop_white_eff <- I_white_eff/P_white_total
  
  asian_c_asian <- aca(t)
  asian_c_black <- acb(t)
  asian_c_other <- aco(t)
  asian_c_white <- acw(t)
  
  black_c_asian <- bca(t)
  black_c_black <- bcb(t)
  black_c_other <- bco(t)
  black_c_white <- bcw(t)
  
  other_c_asian <- oca(t)
  other_c_black <- ocb(t)
  other_c_other <- oco(t)
  other_c_white <- ocw(t)
  
  white_c_asian <- wca(t)
  white_c_black <- wcb(t)
  white_c_other <- wco(t)
  white_c_white <- wcw(t)
  
  ## npi wave 2
  ## The max npi compliance is referenced to the whole population. Given that only non-L move to L, 
  ## if n is the rate moving into npis, then 
  # n_true * P = n_adj * (S1 + E1 + I1 + R1 + S2 + E2 + I2 + R2)
  # n_adj = n_true * P/P_eff
  # maxing at 1
  
  npi_scalar_asian <- P_asian_total/P_asian_eff
  npi_scalar_black <- P_black_total/P_black_eff
  npi_scalar_other <- P_other_total/P_other_eff
  npi_scalar_white <- P_white_total/P_white_eff
  
  npi_a <- npi_weekly(npi_secondwave*npi_asian, npi_halfway, 2)
  npi_b <- npi_weekly(npi_secondwave*npi_black, npi_halfway, 2)
  npi_o <- npi_weekly(npi_secondwave*npi_other, npi_halfway, 2)
  npi_w <- npi_weekly(npi_secondwave*npi_white, npi_halfway, 2)
  
  npi_asian_rate <- min(npi_scalar_asian*npi_a(t-315), 1)
  npi_black_rate <- min(npi_scalar_black*npi_b(t-315), 1)
  npi_other_rate <- min(npi_scalar_other*npi_o(t-315), 1)
  npi_white_rate <- min(npi_scalar_white*npi_w(t-315), 1)
  
  # The vaccination rates are based on data referenced to the whole population. In this model, only susceptible and recovered can get vaccinated (S1, R1, S2, R2). 
  # Thus, the pool to be vaccinated needs to be adjusted. If v is the percent daily that is observed in the data, then
  # v_true*P = v_adj * (S1 + R1 + S2 + R2)
  # so that v_adj = v_true*P/(S1 + R1 + S2 + R2)
  
  v_scalar_asian <- P_asian_total/(S1_asian + S1_asian_L + R1_asian + R1_asian_L + S2_asian + S2_asian_L + R2_asian  + R2_asian_L)
  v_scalar_black <- P_black_total/(S1_black + S1_black_L + R1_black + R1_black_L + S2_black + S2_black_L + R2_black  + R2_black_L)
  v_scalar_other <- P_other_total/(S1_other + S1_other_L + R1_other + R1_other_L + S2_other + S2_other_L + R2_other  + R2_other_L)
  v_scalar_white <- P_white_total/(S1_white + S1_white_L + R1_white + R1_white_L + S2_white + S2_white_L + R2_white  + R2_white_L)
  
  nu_asian <- min(0.85*v_scalar_asian*vax_a(t),1)
  nu_black <- min(0.85*v_scalar_black*vax_b(t), 1)
  nu_other <- min(0.85*v_scalar_other*vax_o(t), 1)
  nu_white <- min(0.85*v_scalar_white*vax_w(t), 1)
  
  var <- variant(t)
  
  lambda1_asian <- var*mu*(aca(t)*prop_asian_eff + acb(t)*prop_black_eff  + aco(t)*prop_other_eff + acw(t)*prop_white_eff)
  lambda2_asian <- ir_2*lambda1_asian
  
  lambda1_black <- var*mu*(bca(t)*prop_asian_eff+ bcb(t)*prop_black_eff + bco(t)*prop_other_eff + bcw(t)*prop_white_eff)
  lambda2_black<- ir_2*lambda1_black
  
  lambda1_other <- var*mu*(oca(t)*prop_asian_eff+ ocb(t)*prop_black_eff + oco(t)*prop_other_eff + ocw(t)*prop_white_eff)
  lambda2_other <- ir_2*lambda1_other
  
  lambda1_white <- var*mu*(wca(t)*prop_asian_eff+ wcb(t)*prop_black_eff + wco(t)*prop_other_eff + wcw(t)*prop_white_eff)
  lambda2_white <- ir_2*lambda1_white
  
  ## Updates
  dS1_asian <- - lambda1_asian*S1_asian + npi_out_asian*S1_asian_L - nu_asian*S1_asian - npi_asian_rate*S1_asian
  dE1_asian <- + lambda1_asian*S1_asian - epsilon*E1_asian  + npi_out_asian*E1_asian_L - npi_asian_rate*E1_asian
  dI1_asian <- + epsilon*E1_asian - gamma*I1_asian + npi_out_asian*I1_asian_L - npi_asian_rate*I1_asian
  dR1_asian <- + (1-sigma1_asian)*gamma*I1_asian - omega*R1_asian + npi_out_asian*R1_asian_L - nu_asian*R1_asian - npi_asian_rate*R1_asian
  dS2_asian <- + omega*R1_asian - lambda2_asian*S2_asian + omega*R2_asian + npi_out_asian*S2_asian_L - nu_asian*S2_asian - npi_asian_rate*S2_asian
  dE2_asian <- + lambda2_asian*S2_asian - epsilon*E2_asian + npi_out_asian*E2_asian_L - npi_asian_rate*E2_asian
  dI2_asian <- + epsilon*E2_asian - gamma*I2_asian + npi_out_asian*I2_asian_L - npi_asian_rate*I2_asian
  dR2_asian <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian - omega*R2_asian + npi_out_asian*R2_asian_L + nu_asian*(S1_asian + R1_asian + S2_asian) - npi_asian_rate*R2_asian
  
  dS1_asian_L <-  - npi_out_asian*S1_asian_L - nu_asian*S1_asian_L - (1-npi_efficacy)*lambda1_asian*S1_asian_L + npi_asian_rate*S1_asian
  dE1_asian_L <-  - epsilon*E1_asian_L  - npi_out_asian*E1_asian_L + (1-npi_efficacy)*lambda1_asian*S1_asian_L + npi_asian_rate*E1_asian
  dI1_asian_L <- + epsilon*E1_asian_L - gamma*I1_asian_L - npi_out_asian*I1_asian_L + npi_asian_rate*I1_asian
  dR1_asian_L <- + (1-sigma1_asian)*gamma*I1_asian_L - omega*R1_asian_L  - npi_out_asian*R1_asian_L - nu_asian*R1_asian_L + npi_asian_rate*R1_asian
  dS2_asian_L <- + omega*R1_asian_L + omega*R2_asian_L - npi_out_asian*S2_asian_L - nu_asian*S2_asian_L - (1-npi_efficacy)*lambda2_asian*S2_asian_L + npi_asian_rate*S2_asian
  dE2_asian_L <- - epsilon*E2_asian_L - npi_out_asian*E2_asian_L + (1-npi_efficacy)*lambda2_asian*S2_asian_L + npi_asian_rate*E2_asian
  dI2_asian_L <- + epsilon*E2_asian_L - gamma*I2_asian_L  - npi_out_asian*I2_asian_L + npi_asian_rate*I2_asian
  dR2_asian_L <- + (1-ifr_2*sigma1_asian)*gamma*I2_asian_L - omega*R2_asian_L - npi_out_asian*R2_asian_L + nu_asian*(S1_asian_L + R1_asian_L + S2_asian_L) + npi_asian_rate*R2_asian
  
  dD_asian <- sigma1_asian*gamma*(I1_asian + I1_asian_L) + ifr_2*sigma1_asian*gamma*(I2_asian + I2_asian_L)
  dInc_asian <- + epsilon*(E1_asian + E2_asian + E1_asian_L + E2_asian_L)
  dRepCases_asian <- rr_asian*epsilon*(E1_asian+ E2_asian + E1_asian_L + E2_asian_L)
  
  dS1_black <- - lambda1_black*S1_black + npi_out_black*S1_black_L - nu_black*S1_black - npi_black_rate*S1_black
  dE1_black <- + lambda1_black*S1_black - epsilon*E1_black  + npi_out_black*E1_black_L - npi_black_rate*E1_black
  dI1_black <- + epsilon*E1_black - gamma*I1_black + npi_out_black*I1_black_L - npi_black_rate*I1_black
  dR1_black <- + (1-sigma1_black)*gamma*I1_black - omega*R1_black + npi_out_black*R1_black_L - nu_black*R1_black - npi_black_rate*R1_black
  dS2_black <- + omega*R1_black - lambda2_black*S2_black + omega*R2_black + npi_out_black*S2_black_L - nu_black*S2_black - npi_black_rate*S2_black
  dE2_black <- + lambda2_black*S2_black - epsilon*E2_black + npi_out_black*E2_black_L - npi_black_rate*E2_black
  dI2_black <- + epsilon*E2_black - gamma*I2_black + npi_out_black*I2_black_L - npi_black_rate*I2_black
  dR2_black <- + (1-ifr_2*sigma1_black)*gamma*I2_black - omega*R2_black + npi_out_black*R2_black_L + nu_black*(S1_black + R1_black + S2_black) - npi_black_rate*R2_black
  
  dS1_black_L <-  - npi_out_black*S1_black_L - nu_black*S1_black_L - (1-npi_efficacy)*lambda1_black*S1_black_L + npi_black_rate*S1_black
  dE1_black_L <-  - epsilon*E1_black_L  - npi_out_black*E1_black_L + (1-npi_efficacy)*lambda1_black*S1_black_L + npi_black_rate*E1_black
  dI1_black_L <- + epsilon*E1_black_L - gamma*I1_black_L - npi_out_black*I1_black_L + npi_black_rate*I1_black
  dR1_black_L <- + (1-sigma1_black)*gamma*I1_black_L - omega*R1_black_L  - npi_out_black*R1_black_L - nu_black*R1_black_L + npi_black_rate*R1_black
  dS2_black_L <- + omega*R1_black_L + omega*R2_black_L - npi_out_black*S2_black_L - nu_black*S2_black_L - (1-npi_efficacy)*lambda2_black*S2_black_L + npi_black_rate*S2_black
  dE2_black_L <- - epsilon*E2_black_L - npi_out_black*E2_black_L + (1-npi_efficacy)*lambda2_black*S2_black_L + npi_black_rate*E2_black
  dI2_black_L <- + epsilon*E2_black_L - gamma*I2_black_L  - npi_out_black*I2_black_L + npi_black_rate*I2_black
  dR2_black_L <- + (1-ifr_2*sigma1_black)*gamma*I2_black_L - omega*R2_black_L - npi_out_black*R2_black_L + nu_black*(S1_black_L + R1_black_L + S2_black_L) + npi_black_rate*R2_black
  
  dD_black <- sigma1_black*gamma*(I1_black + I1_black_L) + ifr_2*sigma1_black*gamma*(I2_black + I2_black_L)
  dInc_black <- + epsilon*(E1_black + E2_black + E1_black_L + E2_black_L)
  dRepCases_black <- rr_black*epsilon*(E1_black+ E2_black + E1_black_L + E2_black_L)
  
  dS1_other <- - lambda1_other*S1_other + npi_out_other*S1_other_L - nu_other*S1_other - npi_other_rate*S1_other
  dE1_other <- + lambda1_other*S1_other - epsilon*E1_other  + npi_out_other*E1_other_L - npi_other_rate*E1_other
  dI1_other <- + epsilon*E1_other - gamma*I1_other + npi_out_other*I1_other_L - npi_other_rate*I1_other
  dR1_other <- + (1-sigma1_other)*gamma*I1_other - omega*R1_other + npi_out_other*R1_other_L - nu_other*R1_other - npi_other_rate*R1_other
  dS2_other <- + omega*R1_other - lambda2_other*S2_other + omega*R2_other + npi_out_other*S2_other_L - nu_other*S2_other - npi_other_rate*S2_other
  dE2_other <- + lambda2_other*S2_other - epsilon*E2_other + npi_out_other*E2_other_L - npi_other_rate*E2_other
  dI2_other <- + epsilon*E2_other - gamma*I2_other + npi_out_other*I2_other_L - npi_other_rate*I2_other
  dR2_other <- + (1-ifr_2*sigma1_other)*gamma*I2_other - omega*R2_other + npi_out_other*R2_other_L + nu_other*(S1_other + R1_other + S2_other) - npi_other_rate*R2_other
  
  dS1_other_L <-  - npi_out_other*S1_other_L - nu_other*S1_other_L - (1-npi_efficacy)*lambda1_other*S1_other_L + npi_other_rate*S1_other
  dE1_other_L <-  - epsilon*E1_other_L  - npi_out_other*E1_other_L + (1-npi_efficacy)*lambda1_other*S1_other_L + npi_other_rate*E1_other
  dI1_other_L <- + epsilon*E1_other_L - gamma*I1_other_L - npi_out_other*I1_other_L + npi_other_rate*I1_other
  dR1_other_L <- + (1-sigma1_other)*gamma*I1_other_L - omega*R1_other_L  - npi_out_other*R1_other_L - nu_other*R1_other_L + npi_other_rate*R1_other
  dS2_other_L <- + omega*R1_other_L + omega*R2_other_L - npi_out_other*S2_other_L - nu_other*S2_other_L - (1-npi_efficacy)*lambda2_other*S2_other_L + npi_other_rate*S2_other
  dE2_other_L <- - epsilon*E2_other_L - npi_out_other*E2_other_L + (1-npi_efficacy)*lambda2_other*S2_other_L + npi_other_rate*E2_other
  dI2_other_L <- + epsilon*E2_other_L - gamma*I2_other_L  - npi_out_other*I2_other_L + npi_other_rate*I2_other
  dR2_other_L <- + (1-ifr_2*sigma1_other)*gamma*I2_other_L - omega*R2_other_L - npi_out_other*R2_other_L + nu_other*(S1_other_L + R1_other_L + S2_other_L) + npi_other_rate*R2_other
  
  dD_other <- sigma1_other*gamma*(I1_other + I1_other_L) + ifr_2*sigma1_other*gamma*(I2_other + I2_other_L)
  dInc_other <- + epsilon*(E1_other + E2_other + E1_other_L + E2_other_L)
  dRepCases_other <- rr_other*epsilon*(E1_other+ E2_other + E1_other_L + E2_other_L)
  
  dS1_white <- - lambda1_white*S1_white + npi_out_white*S1_white_L - nu_white*S1_white - npi_white_rate*S1_white
  dE1_white <- + lambda1_white*S1_white - epsilon*E1_white  + npi_out_white*E1_white_L - npi_white_rate*E1_white
  dI1_white <- + epsilon*E1_white - gamma*I1_white + npi_out_white*I1_white_L - npi_white_rate*I1_white
  dR1_white <- + (1-sigma1_white)*gamma*I1_white - omega*R1_white + npi_out_white*R1_white_L - nu_white*R1_white - npi_white_rate*R1_white
  dS2_white <- + omega*R1_white - lambda2_white*S2_white + omega*R2_white + npi_out_white*S2_white_L - nu_white*S2_white - npi_white_rate*S2_white
  dE2_white <- + lambda2_white*S2_white - epsilon*E2_white + npi_out_white*E2_white_L - npi_white_rate*E2_white
  dI2_white <- + epsilon*E2_white - gamma*I2_white + npi_out_white*I2_white_L - npi_white_rate*I2_white
  dR2_white <- + (1-ifr_2*sigma1_white)*gamma*I2_white - omega*R2_white + npi_out_white*R2_white_L + nu_white*(S1_white + R1_white + S2_white) - npi_white_rate*R2_white
  
  dS1_white_L <-  - npi_out_white*S1_white_L - nu_white*S1_white_L - (1-npi_efficacy)*lambda1_white*S1_white_L + npi_white_rate*S1_white
  dE1_white_L <-  - epsilon*E1_white_L  - npi_out_white*E1_white_L + (1-npi_efficacy)*lambda1_white*S1_white_L + npi_white_rate*E1_white
  dI1_white_L <- + epsilon*E1_white_L - gamma*I1_white_L - npi_out_white*I1_white_L + npi_white_rate*I1_white
  dR1_white_L <- + (1-sigma1_white)*gamma*I1_white_L - omega*R1_white_L  - npi_out_white*R1_white_L - nu_white*R1_white_L + npi_white_rate*R1_white
  dS2_white_L <- + omega*R1_white_L + omega*R2_white_L - npi_out_white*S2_white_L - nu_white*S2_white_L - (1-npi_efficacy)*lambda2_white*S2_white_L + npi_white_rate*S2_white
  dE2_white_L <- - epsilon*E2_white_L - npi_out_white*E2_white_L + (1-npi_efficacy)*lambda2_white*S2_white_L + npi_white_rate*E2_white
  dI2_white_L <- + epsilon*E2_white_L - gamma*I2_white_L  - npi_out_white*I2_white_L + npi_white_rate*I2_white
  dR2_white_L <- + (1-ifr_2*sigma1_white)*gamma*I2_white_L - omega*R2_white_L - npi_out_white*R2_white_L + nu_white*(S1_white_L + R1_white_L + S2_white_L) + npi_white_rate*R2_white
  
  dD_white <- sigma1_white*gamma*(I1_white + I1_white_L) + ifr_2*sigma1_white*gamma*(I2_white + I2_white_L)
  dInc_white <- + epsilon*(E1_white + E2_white + E1_white_L + E2_white_L)
  dRepCases_white <- rr_white*epsilon*(E1_white+ E2_white + E1_white_L + E2_white_L)
  
  # Return list of gradients
  list(c(dS1_asian,dE1_asian,dI1_asian,dR1_asian,
         dS2_asian,dE2_asian,dI2_asian,dR2_asian,
         
         dS1_asian_L,dE1_asian_L,dI1_asian_L,dR1_asian_L,
         dS2_asian_L,dE2_asian_L,dI2_asian_L,dR2_asian_L,
         
         dD_asian, dInc_asian, dRepCases_asian,
         
         dS1_black,dE1_black,dI1_black,dR1_black,
         dS2_black,dE2_black,dI2_black,dR2_black,
         
         dS1_black_L,dE1_black_L,dI1_black_L,dR1_black_L,
         dS2_black_L,dE2_black_L,dI2_black_L,dR2_black_L,
         
         dD_black, dInc_black, dRepCases_black,
         
         dS1_other,dE1_other,dI1_other,dR1_other,
         dS2_other,dE2_other,dI2_other,dR2_other,
         
         dS1_other_L,dE1_other_L,dI1_other_L,dR1_other_L,
         dS2_other_L,dE2_other_L,dI2_other_L,dR2_other_L,
         
         dD_other, dInc_other, dRepCases_other,
         
         dS1_white,dE1_white,dI1_white,dR1_white,
         dS2_white,dE2_white,dI2_white,dR2_white,
         
         dS1_white_L,dE1_white_L,dI1_white_L,dR1_white_L,
         dS2_white_L,dE2_white_L,dI2_white_L,dR2_white_L,
         
         dD_white, dInc_white, dRepCases_white))
}








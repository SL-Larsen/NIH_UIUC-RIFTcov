library(deSolve)
library(tidyverse)

#### General parameters ####

## CA
Pop_T_CA <- 39346023
Pop_Asian_CA <- 5743983
Pop_Black_CA <- 2142371
Pop_Latino_CA <- 15380929
Pop_White_CA <- 14365145
Pop_Other_CA <- 1713595

## NC
Pop_T_NC <- 10698973
Pop_Asian_NC <- 341052
Pop_Black_NC <- 2155650
Pop_White_NC <- 6497519
Pop_Other_NC <- 1704752


names_var <- c("S1_asian","E1_asian","I1_asian","R1_asian",
               "S2_asian","E2_asian","I2_asian","R2_asian",
               "S1_asian_L","E1_asian_L","I1_asian_L","R1_asian_L",
               "S2_asian_L","E2_asian_L","I2_asian_L","R2_asian_L",
               "D_asian", "Inc_asian", "RepCases_asian",
               
               "S1_black","E1_black","I1_black","R1_black",
               "S2_black","E2_black","I2_black","R2_black",
               "S1_black_L","E1_black_L","I1_black_L","R1_black_L",
               "S2_black_L","E2_black_L","I2_black_L","R2_black_L",              
               "D_black", "Inc_black","RepCases_black",
               
               "S1_latino","E1_latino","I1_latino","R1_latino",
               "S2_latino","E2_latino","I2_latino","R2_latino",
               "S1_latino_L","E1_latino_L","I1_latino_L","R1_latino_L",
               "S2_latino_L","E2_latino_L","I2_latino_L","R2_latino_L",              
               "D_latino", "Inc_latino","RepCases_latino",
               
               "S1_other","E1_other","I1_other","R1_other",
               "S2_other","E2_other","I2_other","R2_other",
               "S1_other_L","E1_other_L","I1_other_L","R1_other_L",
               "S2_other_L","E2_other_L","I2_other_L","R2_other_L",
               "D_other", "Inc_other", "RepCases_other",
              
               "S1_white","E1_white","I1_white","R1_white",
               "S2_white","E2_white","I2_white","R2_white",
               "S1_white_L","E1_white_L","I1_white_L","R1_white_L",
               "S2_white_L","E2_white_L","I2_white_L","R2_white_L",
               "D_white", "Inc_white", "RepCases_white")

names_var_nolatino <- c("S1_asian","E1_asian","I1_asian","R1_asian",
               "S2_asian","E2_asian","I2_asian","R2_asian",
               "S1_asian_L","E1_asian_L","I1_asian_L","R1_asian_L",
               "S2_asian_L","E2_asian_L","I2_asian_L","R2_asian_L",
               "D_asian", "Inc_asian", "RepCases_asian",
               
               "S1_black","E1_black","I1_black","R1_black",
               "S2_black","E2_black","I2_black","R2_black",
               "S1_black_L","E1_black_L","I1_black_L","R1_black_L",
               "S2_black_L","E2_black_L","I2_black_L","R2_black_L",              
               "D_black", "Inc_black","RepCases_black",
               
               
               "S1_other","E1_other","I1_other","R1_other",
               "S2_other","E2_other","I2_other","R2_other",
               "S1_other_L","E1_other_L","I1_other_L","R1_other_L",
               "S2_other_L","E2_other_L","I2_other_L","R2_other_L",
               "D_other", "Inc_other", "RepCases_other",
               
               "S1_white","E1_white","I1_white","R1_white",
               "S2_white","E2_white","I2_white","R2_white",
               "S1_white_L","E1_white_L","I1_white_L","R1_white_L",
               "S2_white_L","E2_white_L","I2_white_L","R2_white_L",
               "D_white", "Inc_white", "RepCases_white")

# sigma are set to the case fatality rate for each group during calibration for phase 1
CA_parms <- c(epsilon = 1/3, omega = 1/180, gamma = 1/5,
              ir_2 = 0.5, ifr_2 = 0.25,
              npi_efficacy = 1,
              
              mu = 0.03, 
              rr_asian = 0.3, rr_black = 0.3, rr_latino = 0.3, 
              rr_other = 0.3, rr_white = 0.3,
              
              npi_out_asian = 1/7, npi_out_black = 1/7, npi_out_latino = 1/7,
              npi_out_other = 1/7, npi_out_white = 1/7, 
              
              sigma1_asian = 	0.049739437, sigma1_black = 0.040823002, sigma1_latino = 0.019558332, 
              sigma1_other = 0.001036306, sigma1_white = 0.038037242,
              
              npi_asian = 0, npi_black = 0, npi_latino = 0, npi_other = 0, npi_white = 0,
              npi_secondwave = 0, npi_halfway = 10
              )


# sigma are set to the case fatality rate for each group during calibration for phase 1
NC_parms <- c(epsilon = 1/3, omega = 1/180, gamma = 1/5,
              ir_2 = 0.5, ifr_2 = 0.25,
              npi_efficacy = 1,

              mu = 0.03, 
              rr_asian = 0.3, rr_black = 0.3,
              rr_other = 0.3, rr_white = 0.3,
              
              npi_out_asian = 1/93.2, npi_out_black = 1/121,
              npi_out_other = 1/83, npi_out_white = 1/220,
              
              sigma1_asian = 0.006012413, sigma1_black = 0.027064695, 
              sigma1_other = 0.007673184, sigma1_white = 0.018500013,
              
              npi_asian = 0, npi_black = 0, npi_other = 0, npi_white = 0,
              npi_secondwave = 0, npi_halfway = 10
              )

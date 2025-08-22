#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 8/22/2025 #
#         Modified on 8/22/2025          #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")

# standardize wetland data (invert and birds already standardized)
wetland <- wetland %>%
  mutate(across(where(is.numeric), scale))

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

birds <- birds %>% 
  mutate(PlasmaDetection = case_when(
    PlasmaDetection == "Y" ~ 1,
    PlasmaDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  InvertPesticideDetection = case_when(
    InvertPesticideDetection == "Y" ~ 1,
    InvertPesticideDetection == "N" ~ 0,
    TRUE ~ NA_real_), 
  AnyDetection = case_when(
    AnyDetection == "Y" ~ 1,
    AnyDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ 1,
    Permanence == "Semipermanent" ~ 2,
    Permanence == "Permanent" ~ 3,
    TRUE ~ NA_real_))


invert <- invert %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Buffered = case_when(
      Buffered == "Y" ~ 1,
      Buffered == "N" ~ 0,
      TRUE ~ NA_real_))

  
wetland <- wetland %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  InvertPesticideDetection = case_when(
    InvertPesticideDetection == "Y" ~ 1,
    InvertPesticideDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ 1,
    Permanence == "Semipermanent" ~ 2,
    Permanence == "Permanent" ~ 3,
    TRUE ~ NA_real_),
  EnvDetection = case_when(
    EnvDetection == "Y" ~ 1,
    EnvDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  Buffered = case_when(
    Buffered == "Y" ~ 1,
    Buffered == "N" ~ 0,
    TRUE ~ NA_real_))
  

#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

# ...plasma neonic detection (pd) ----
pd.percentag <- glm(PlasmaDetection ~ PercentAg, family = "binomial",
                        data = birds)

pd.


# ...water neonic detection (wd) ----


# ...environmental pesticide detection (ed) ----

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 



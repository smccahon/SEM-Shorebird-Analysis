#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 8/22/2025 #
#         Modified on 8/25/2025          #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)

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
  EnvDetection = case_when(
    EnvDetection == "Y" ~ 1,
    EnvDetection == "N" ~ 0,
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
#                         variable correlations                             ----                        
#------------------------------------------------------------------------------# 

# drought and precipitation
cor(birds$SPEI, birds$AnnualSnowfall_in) # 0.30
cor(birds$SPEI, birds$DaysSinceLastPrecipitation_5mm) # -0.51
cor(birds$SPEI, birds$PrecipitationAmount_7days) # 0.52

# julian and season
birds$Event_num <- as.numeric(as.factor(birds$Event))
cor(birds$Event_num, birds$Julian) # -0.78

# local vegetation and buffer presence
cor(invert$PercentLocalVeg_50m, invert$Buffered) # 0.24

# buffer presence and % cropland cover
cor(invert$Buffered, invert$PercentAg) #-0.75

# local vegetation and % cropland cover
cor(invert$PercentLocalVeg_50m, invert$PercentAg) # -0.27

# biomass and diversity
cor(birds$Biomass, birds$Diversity, use = "complete.obs") # 0.38

# water quality and pH
cor(invert$WaterQuality, invert$pH_probe) # 0.27

  
#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

# check about scaling!! scaled predictors, unscaled responses? do I scale before?

m1 <- glm(EnvDetection ~ PercentAg + SPEI + AnnualSnowfall_in + 
            DaysSinceLastPrecipitation_5mm + Julian, 
            family = "binomial", data = wetland)

m2 <- glm(PlasmaDetection ~ EnvDetection + seconds_since_midnight + PercentAg +
                SPEI + Julian + Species, family = "binomial", data = birds)

m3 <- lm(Biomass ~ PercentAg + Julian + EnvDetection +
           PercentLocalVeg_50m + WaterQuality + pH_probe, data = invert)

m4 <- lm(BodyCondition ~ Species + seconds_since_midnight + Biomass +
           Diversity + Julian, 
         data = birds)

m5 <- lm(PercentLocalVeg_50m ~ PercentAg + WaterQuality + Julian + pH_probe +
           AnnualSnowfall_in, 
         data = invert)

# won't run as is due to scaling issue
m6 <- glm(Diversity ~ PercentAg + Julian + EnvDetection +
           PercentLocalVeg_50m + WaterQuality + pH_probe, data = invert,
          family = "poisson")


#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

model <- psem(m1, m2, m3)
summary(model, conserve = TRUE)

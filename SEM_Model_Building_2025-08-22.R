#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 8/22/2025 #
#         Modified on 8/26/2025          #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)
library(multcompView)
library(cplm)
library(DHARMa)

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")

# filter to only include 2023 data
birds <- birds %>% 
  filter(Event %in% c("Fall 2023", "Spring 2023")) # 126 birds

wetland <- wetland %>% 
  filter(Year == "2023") # 79 wetland surveys surveys

# unresolved issues:
# Julian or season? Julian seems to cause more problems but its continuous
# Dominant crop or % cropland cover? Could run both models and interchange
# variables and use AIC to help decide
# permanence as factor or ordinal variable?
# might need to handle Diversity differently (could be zero-inflated)

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
  Season = case_when(
    Season == "Spring" ~ 1,
    Season == "Fall" ~ 0),
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
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0),
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
  EnvDetection = case_when(
    EnvDetection == "Y" ~ 1,
    EnvDetection == "N" ~ 0,
    TRUE ~ NA_real_),
  Season = case_when(
    Season == "Spring" ~ 1,
    Season == "Fall" ~ 0),
  Buffered = case_when(
    Buffered == "Y" ~ 1,
    Buffered == "N" ~ 0,
    TRUE ~ NA_real_))

#------------------------------------------------------------------------------#
#                         variable correlations                             ----                        
#------------------------------------------------------------------------------# 

# julian and season
cor(birds$Season, birds$Julian) # -0.95

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

# dominant crop type and % cropland cover
summary(aov(data = wetland, PercentAg ~ DominantCrop)) # p < 0.001

  
#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

wetland.1 <- wetland %>% 
  mutate(across(c(AnnualSnowfall_in, PercentAg, DaysSinceLastPrecipitation_5mm), 
                scale))

m1 <- glm(EnvDetection ~ PercentAg + AnnualSnowfall_in + 
            DaysSinceLastPrecipitation_5mm + Season + Permanence, 
            family = "binomial", data = wetland.1)

#---

birds.2 <- birds %>% 
  mutate(across(c(seconds_since_midnight, PercentAg, Julian), 
                scale))

m2 <- glmmTMB(PlasmaDetection ~ EnvDetection + seconds_since_midnight + 
                PercentAg + Julian + (1 | Site), family = "binomial", 
              data = birds.2)

#---

invert.3 <- invert %>% 
  mutate(across(c(PercentLocalVeg_50m, PercentAg, Julian, pH_probe), 
                scale))

m3 <- cpglm(Biomass ~ PercentAg + Julian + EnvDetection +
           PercentLocalVeg_50m + WaterQuality + pH_probe, data = invert.3,
         link = "log")

#---

birds.4 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, Julian), 
                scale))

m4 <- lm(BodyCondition ~ seconds_since_midnight + Biomass +
           Diversity + Julian, 
         data = birds.4)

#---

invert.5 <- invert %>% 
  mutate(across(c(PercentAg, pH_probe, Julian, AnnualSnowfall_in), 
                scale))

m5 <- lm(PercentLocalVeg_50m ~ PercentAg + WaterQuality + Julian + pH_probe +
           AnnualSnowfall_in, 
         data = invert.5)

#---

invert.6 <- invert %>% 
  mutate(across(c(PercentAg, PercentLocalVeg_50m, Julian, pH_probe), 
                scale))

m6 <- glm(Diversity ~ PercentAg + Julian + EnvDetection +
           PercentLocalVeg_50m + WaterQuality + pH_probe, data = invert.6,
          family = "poisson")

#---

birds.7 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, Julian), 
                scale))

m7 <- lm(FatteningIndex ~ seconds_since_midnight + Biomass +
           Diversity + Julian, 
         data = birds.7)

# continue building the rest of the models and then check model diagnostics

#---

model <- psem(m1, m2, m3, m4, m5, m6, m7)
summary(model)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# ...environmental detection model (GOOD) ----
overdispersion <- sum(resid(m1, type = "pearson")^2) / df.residual(m1)
overdispersion  # 1.06

DHARM <- simulateResiduals(fittedModel = m1)
plot(DHARM)
testDispersion(DHARM) # overdispersion: 1.02
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation

# ...plasma detection model (GOOD) ----
DHARM <- simulateResiduals(fittedModel = m2)
plot(DHARM)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation


# ...biomass model (GOOD) ----
plot(residuals(m3), main = "Deviance Residuals", ylab = "Residuals", 
     xlab = "Index") # no clear pattern

#...body condition model ----
DHARM <- simulateResiduals(fittedModel = m4)
plot(DHARM) # combined adjusted quantile test significant
plotResiduals(DHARM, form = birds.4$Julian)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation






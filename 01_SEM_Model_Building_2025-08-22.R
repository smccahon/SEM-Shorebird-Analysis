#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 8/22/2025 #
#         Modified on 8/29/2025          #
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
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
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
    TRUE ~ NA_real_),
  Permanence = case_when(
    Permanence %in% c("Temporary", "Seasonal") ~ 1,
    Permanence == "Semipermanent" ~ 2,
    Permanence == "Permanent" ~ 3,
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
  mutate(across(c(AnnualSnowfall_in, PercentAg, DaysSinceLastPrecipitation_5mm,
                  SPEI), 
               ~ as.numeric(scale(.))))

m1 <- glm(EnvDetection ~ PercentAg + AnnualSnowfall_in + 
            DaysSinceLastPrecipitation_5mm + Season + Permanence + SPEI, 
            family = "binomial", data = wetland.1)

#---

birds.2 <- birds %>% 
  mutate(across(c(seconds_since_midnight, PercentAg, Julian, SPEI), 
                ~ as.numeric(scale(.))))

m2 <- glmmTMB(PlasmaDetection ~ EnvDetection + seconds_since_midnight + 
                PercentAg + Julian + SPEI + (1 | Site), family = "binomial", 
              data = birds.2)

#---

invert.3 <- invert %>% 
  mutate(across(c(PercentLocalVeg_50m, PercentAg, Julian), 
                ~ as.numeric(scale(.))))

m3 <- glm(Biomass ~ PercentAg + Julian + EnvDetection +
           PercentLocalVeg_50m + WaterQuality + Permanence, data = invert.3)

#---

birds.4 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, Julian, SPEI), 
                ~ as.numeric(scale(.))))

m4 <- lmer(BodyCondition ~ seconds_since_midnight + Biomass + Julian + SPEI +
           PlasmaDetection + (1 | Site), 
         data = birds.4)

#---

invert.5 <- invert %>% 
  mutate(across(c(PercentAg, Julian, AnnualSnowfall_in), 
                ~ as.numeric(scale(.))))

m5 <- lm(PercentLocalVeg_50m ~ PercentAg + Julian + Permanence +
           AnnualSnowfall_in, 
         data = invert.5)

#---

birds.6 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, SPEI, Julian),
                ~ as.numeric(scale(.))))


m6 <- lmer(FatteningIndex ~ seconds_since_midnight + Biomass + PlasmaDetection +
           SPEI + Julian + (1 | Site), data = birds.6)

#---

wetland.7 <- wetland %>% 
  mutate(across(c(AnnualSnowfall_in, Julian),
                ~ as.numeric(scale(.))))

m7 <- lm(SPEI ~ AnnualSnowfall_in + Julian, data = wetland.7)
  

#---

wetland.8 <- wetland %>% 
  mutate(across(c(AnnualSnowfall_in, Julian, SPEI),
                ~ as.numeric(scale(.))))

m8 <- glm(Permanence ~ AnnualSnowfall_in + Julian + SPEI, data = wetland.8,
          family = "poisson")

#---

invert.9 <- invert %>% 
  mutate(across(c(PercentAg, PercentLocalVeg_50m),
                ~ as.numeric(scale(.))))

m9 <- lm(WaterQuality ~ PercentLocalVeg_50m + PercentAg, data = invert.9)

#---

model <- psem(m1, m2, m3, m4, m5, m6, m7, m8, m9)
summary(model)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# ...environmental detection model (GOOD) ----
overdispersion <- sum(resid(m1, type = "pearson")^2) / df.residual(m1)
overdispersion  # 1.08

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






#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/16/2025          #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)
library(multcompView)
library(DHARMa)
library(car)
library(statmod)

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-08-29.csv")

# combine species into bill length groupings
birds <- birds %>%
  mutate(Group = case_when(
    Species %in% c("Marbled Godwit", "American Avocet", "Shortbilled Dowitcher",
                   "Longbilled Dowitcher", "Greater Yellowlegs",
                   "Willet") ~ "Long",
    Species %in% c("Lesser Yellowlegs", "Pectoral Sandpiper", 
                   "Wilsons Phalarope") ~ "Medium",
    Species %in% c("Least Sandpiper", "Killdeer", 
                   "Semipalmated Sandpiper") ~ "Short")) %>%
  mutate(Group = factor(Group, levels = c("Short", "Medium", "Long")))

# filter to only include 2023 data
birds <- birds %>% 
  filter(Event %in% c("Fall 2023", "Spring 2023")) # 126 birds

wetland <- wetland %>% 
  filter(Year == "2023") # 79 wetland surveys surveys



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

full <- full %>% 
  mutate(PlasmaDetection = case_when(
    PlasmaDetection == "Y" ~ 1,
    PlasmaDetection == "N" ~ 0,
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

#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

# biomass model

invert$log_PercentAg <- log(invert$PercentAg)

# transforming 0's to small numbers to allow the use of a Gamma distribution
invert$Biomass_adj <- ifelse(invert$Biomass == 0, 0.0001, 
                             invert$Biomass)

birds$Biomass_adj <- ifelse(birds$Biomass == 0, 0.0001, 
                            birds$Biomass)

m1 <- glm(Biomass_adj ~ PercentAg,
              family = Gamma(link = "log"), data = invert)

m2 <- lm(BodyCondition ~ Biomass_adj, data = birds)


# psem expects a random effect in glmmTMB so I created a dummy random effect
invert$dummy <- 1

# tweedie model (psem does not support)
  # m1 <- glmmTMB(Biomass ~ log_PercentAg + Julian + (1 | dummy), 
  #               data = invert, 
  #               family = glmmTMB::tweedie(link = "log"))
  


# zero-inflated Gamma model (psem does not support)
 # m1 <- glmmTMB(Biomass ~ log_PercentAg + Julian + (1 | dummy),
 #                 ziformula = ~ log_PercentAg + Julian + (1 | dummy),
 #                 family = ziGamma(link = "log"),
 #                 data = invert)


# plasma detection model
m1 <- lm(BodyCondition ~ PlasmaDetection, data = birds)

m2 <- glmmTMB(PlasmaDetection ~ PercentAg + Biomass + (1|Site), data = birds, 
              family = binomial(link = "logit"))


# body condition model
m3 <- glmmTMB(BodyCondition ~ PlasmaDetection + (1|Site),
              data = birds)


# gamma model runs...but need to calculate standardized estimates manually
model <- psem(m1, m2, m3)
print(model)
summary(model, conserve = TRUE)

#------------------------------------------------------------------------------#
#                     model diagnostics with DHARMa                         ----                        
#------------------------------------------------------------------------------# 


# m1 --- 
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testZeroInflation(m1) # good
testUniformity(simulationOutput) # good,
testOutliers(simulationOutput) # good
testQuantiles(simulationOutput) # good, but step failure warning

plotResiduals(simulationOutput, form = invert$Julian) # issue
plotResiduals(simulationOutput, form = invert$PercentAg) #  good

# m2 --- 
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput) # good,
testOutliers(simulationOutput) # good
testQuantiles(simulationOutput) # good

plotResiduals(simulationOutput, form = birds.cs$Species) # good
plotResiduals(simulationOutput, form = birds$PercentAg) # mostly good






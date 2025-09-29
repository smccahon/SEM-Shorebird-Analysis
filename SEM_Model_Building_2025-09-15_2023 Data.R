#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/29/2025          #
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
library(AICcmodavg)
library(cplm)
library(statmod)

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")

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


#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# model wetlands with biomass > 0 only, given tweedie is not supported (n = 66)
invert.pos <- subset(invert, Biomass > 0)

# saturated model
m1 <- glm(Biomass ~ PercentAg + EnvDetection + WaterQuality + 
            PercentLocalVeg_50m + Season, data = invert.pos,
          na.action = na.omit,
          family = Gamma(link = "log"))

# extract standardized coefficients manually
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert.pos$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std

#---

# ...plasma detection model ----

var <- c("PlasmaDetection", "PercentAg", "EnvDetection")

# Subset data to complete cases for those vars
birds_complete_plasma <- birds[complete.cases(birds[, var]), ]

# saturated
m2 <- glm(PlasmaDetection ~ PercentAg + EnvDetection,
            data = birds_complete_plasma,
            family = binomial(link = "logit"))

#---

# ...body condition model ----

# saturated model
m3 <- lm(BCI ~ Biomass + PercentAg + SPEI + PlasmaDetection +
           time_hours + EnvDetection,
         data = birds,
         na.action = na.omit)

# ...fattening index model ----

# saturated model
m4 <- lm(FatteningIndex ~ Biomass)






model <- psem(m1, m2, m3)
print(model)
summary(model, conserve = TRUE)


#------------------------------------------------------------------------------#
#                     model diagnostics with DHARMa                         ----                        
#------------------------------------------------------------------------------# 


# m1 --- GOOD
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testZeroInflation(m1)
testUniformity(simulationOutput) 
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = invert.pos$Season)
plotResiduals(simulationOutput, form = invert.pos$PercentAg) # not perfect, but okay

# m2 --- GOOD
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = birds_complete_plasma$PercentAg)
plotResiduals(simulationOutput, form = birds_complete_plasma$EnvDetection)

# m3 ---
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

birds_complete <- na.omit(birds[, c("FatteningIndex", "Season", "Group", 
                                    "PercentAg", "Species")])

plotResiduals(simulationOutput, form = birds_complete$Group)





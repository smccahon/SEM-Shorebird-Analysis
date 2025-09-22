#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/22/2025          #
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

# reduced model
# m1 <- glm(Biomass ~ PercentAg , data = invert.pos,
#           na.action = na.omit,
#           family = Gamma(link = "log"))

# extract standardized coefficients manually
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert.pos$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std

#---

# ...plasma detection model----

# var <- c("PlasmaDetection", "PercentAg", "EnvDetection", "Season")
# 
# Subset data to complete cases for those vars
# birds_complete <- birds[complete.cases(birds[, var]), ]

# saturated
m2 <- glm(PlasmaDetection ~ PercentAg + EnvDetection,
            data = birds,
            family = binomial(link = "logit"))

#---

# body condition model

# saturate
m3 <- lm(SMI ~ Biomass + PercentAg + SPEI + PlasmaDetection +
           time_hours + Julian + EnvDetection,
         data = birds,
         na.action = na.omit)

ggplot(birds, aes(x = Biomass, y = SMI)) + geom_point()



# model_names <- paste0("m", 3:4)
# 
# models <- mget(model_names)
# 
# aictab(models, modnames = model_names)


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

# m2 --- not good
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = birds_complete$PercentAg)

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





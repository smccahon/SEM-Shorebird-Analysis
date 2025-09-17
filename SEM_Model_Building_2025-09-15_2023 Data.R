#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/17/2025          #
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

# fix numeric instability issue
birds$time_hours <- birds$seconds_since_midnight / 3600 

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

# ...biomass model ----

# model wetlands with biomass > 0 only (n = 66)
invert.pos <- subset(invert, Biomass > 0)

# saturated

m1 <- glm(Biomass ~ PercentAg + EnvDetection + WaterQuality + 
            PercentLocalVeg_50m + Season, data = invert.pos,
          na.action = na.omit,
          family = Gamma(link = "log"))

#---

# ...plasma detection model----

var <- c("PlasmaDetection", "PercentAg", "EnvDetection", "Season")

# 2. Subset data to complete cases for those vars
birds_complete <- birds[complete.cases(birds[, var]), ]

# saturated
m2 <- glm(PlasmaDetection ~ PercentAg + EnvDetection,
            data = birds_complete,
            family = binomial(link = "logit"))

#---

# body condition model

# saturated (without Species/Group)

m3 <- lm(FatteningIndex ~ Biomass + PercentAg + SPEI + PlasmaDetection +
           time_hours + Species,
         data = birds,
         na.action = na.omit)

m3 <- lm(FatteningIndex ~ Group,
         data = birds,
         na.action = na.omit)


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


# m1 --- GOOD, % ag not perfect within saturated model though
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testZeroInflation(m1)
testUniformity(simulationOutput) 
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = invert.pos$Season)
plotResiduals(simulationOutput, form = invert.pos$PercentAg)

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

birds_complete <- na.omit(birds[, c("Mass", "Season", "Group", 
                                    "PercentAg", "Species")])

plotResiduals(simulationOutput, form = birds_complete$Group)

#------------------------------------------------------------------------------#
#             manually calculate standardized coefficients                  ----                        
#------------------------------------------------------------------------------# 

# --- 1. Biomass ~ PercentAg + Season (Gamma GLM with log link) ---
m1 <- glm(Biomass ~ PercentAg + Season, data = invert.pos,
          family = Gamma(link = "log"))

coef_est1 <- coef(summary(m1))[, "Estimate"]
sd_PercentAg_1 <- sd(invert.pos$PercentAg, na.rm = TRUE)
sd_Season_1 <- sd(invert.pos$Season, na.rm = TRUE)

eta1 <- predict(m1, type = "link")
sd_eta1 <- sd(eta1)

std_beta_PercentAg_1 <- coef_est1["PercentAg"] * sd_PercentAg_1 / sd_eta1
std_beta_Season_1 <- coef_est1["Season"] * sd_Season_1 / sd_eta1

# --- 2. PlasmaDetection ~ PercentAg + Biomass (Binomial GLMM, logit link) ---
m2 <- glmmTMB(PlasmaDetection ~ PercentAg + Biomass + (1|Site), data = birds,
              family = binomial(link = "logit"))

coef_est2 <- fixef(m2)$cond
sd_PercentAg_2 <- sd(birds$PercentAg, na.rm = TRUE)
sd_Biomass_2 <- sd(birds$Biomass, na.rm = TRUE)

# Use fixed latent SD for logistic link
sd_latent_logit <- pi / sqrt(3)  # ~1.814

std_beta_PercentAg_2 <- coef_est2["PercentAg"] * sd_PercentAg_2 / sd_latent_logit
std_beta_Biomass_2 <- coef_est2["Biomass"] * sd_Biomass_2 / sd_latent_logit

# --- 3. BodyCondition ~ PlasmaDetection + Biomass (Gaussian mixed model) ---
m3 <- glmmTMB(BodyCondition ~ PlasmaDetection + Biomass + (1|Site), data = birds)

coef_est3 <- fixef(m3)$cond
sd_PlasmaDetection_3 <- sd(birds$PlasmaDetection, na.rm = TRUE)
sd_Biomass_3 <- sd(birds$Biomass, na.rm = TRUE)
sd_BodyCondition_3 <- sd(birds$BodyCondition, na.rm = TRUE)

std_beta_PlasmaDetection_3 <- coef_est3["PlasmaDetection"] * sd_PlasmaDetection_3 / sd_BodyCondition_3
std_beta_Biomass_3 <- coef_est3["Biomass"] * sd_Biomass_3 / sd_BodyCondition_3

# --- Combine all standardized coefficients ---
std_coefs <- list(
  Biomass = c(PercentAg = std_beta_PercentAg_1, Season = std_beta_Season_1),
  PlasmaDetection = c(PercentAg = std_beta_PercentAg_2, Biomass = std_beta_Biomass_2),
  BodyCondition = c(PlasmaDetection = std_beta_PlasmaDetection_3, Biomass = std_beta_Biomass_3)
)

print(std_coefs)




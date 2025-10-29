#----------------------------------------#
#           SEM Model Building           #
#          Full Dataset Overlap          #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 10/27/2025         #
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

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

# birds: each row is an individual bird (2021-2023)
# invert: each row is a biomass estimate from each wetland (2023)
# wetland: each row is a wetland (2021-2023)
# full: each row is an individual bird from 2023 (same dataset as birds subsetted to 2023)
# birds and full kept separate because FI and BCI were calculated respectively
# for each dataset
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("cleaned_data/wetland_data_cleaned_2025-09-30.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")

# filter to only include 2023 data
birds.2023 <- birds %>% 
  filter(Event %in% c("Fall 2023", "Spring 2023")) # 122 birds

wetland.2023 <- wetland %>% 
  filter(Year == "2023") # 79 wetland surveys surveys

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))


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
    MigStatus = case_when(
      MigStatus == "Migratory" ~ 1,
      MigStatus == "Resident" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

birds.2023 <- birds.2023 %>% 
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
    MigStatus = case_when(
      MigStatus == "Migratory" ~ 1,
      MigStatus == "Resident" ~ 0,
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
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

wetland.2023 <- wetland.2023 %>% 
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
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

full <- full %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    PlasmaDetection = case_when(
      PlasmaDetection == "Y" ~ 1,
      PlasmaDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0),
    MigStatus = case_when(
      MigStatus == "Migratory" ~ 1,
      MigStatus == "Resident" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# model wetlands with biomass > 0 only, given tweedie is not supported (n = 66)
# invert.pos <- subset(invert, Biomass > 0)

# can we use full invert dataset for this path? answer is no (very biased estimate)
# proceed with full dataset from "full" (all biomass estimates > 0 by default)
# m1 <- glm(Biomass ~ PercentAg, data = invert.pos,
#           family = Gamma(link = "log")) # B = -2.48
# 
# m2 <- glm(Biomass ~ PercentAg, data = site_data,
#           family = Gamma(link = "log")) # B = -0.04

# extract one row per site to avoid pseudoreplication in analysis (n = 11 wetlands)
site_data <- full %>%
  distinct(Site, Biomass, PercentAg, Season)

# gamma distribution or log-transformation? We know with the full dataset that
# gamma is a better fit so use this for consistency with other analysis

#### ...final biomass model----
m1 <- glm(Biomass ~ PercentAg, data = site_data,
          family = Gamma(link = "log"))

# we know season is influential, should it be included despite overfitting?
# answer is no, model without season is a much better model due to low sample size
# across seasons (3 in spring, 8 in fall)

# m2 <- glm(Biomass ~ PercentAg + Season, data = site_data,
#           family = Gamma(link = "log"))
# 
# ...model comparison
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

# view individual relationships
# ggplot(site_data, aes(x = PercentAg, y = Biomass)) + 
#  geom_point() + my_theme

# extract standardized coefficients manually from gamma distribution
  beta <- coef(m1)["PercentAg"]
  sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
                 trigamma(1 / summary(m1)$dispersion)) # observation-level variance
  sd_x <- sd(site_data$PercentAg)
  beta_std <- beta * (sd_x / sd_y)
  beta_std # -0.541 is standardized estimate


#---

# ...plasma detection model ----

# model with complete dataset -- SPEI removed because only 1 year of data
  #### final plasma detection model ----
m2 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
             MigStatus,
             data = full,
             na.action = na.omit,
             family = binomial(link = "logit"))
  
  
# season or julian? season is a better fit...Julian may make more sense but
# use season for consistency in SEM. either way, effect is not significant.
# m1 <- glm(PlasmaDetection ~ PercentAg + Julian + EnvDetection + time_hours + 
#             MigStatus,
#             data = full,
#             na.action = na.omit,
#             family = binomial(link = "logit"))
# 
# m2 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
#             MigStatus,
#           data = full,
#           na.action = na.omit,
#           family = binomial(link = "logit"))

# ...model comparison
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)
  
# can we use full bird dataset for this path? yes...but some bias
# m1 <- glm(PlasmaDetection ~ PercentAg + Julian + EnvDetection + time_hours + 
#               MigStatus + SPEI,
#             data = full,
#             na.action = na.omit,
#             family = binomial(link = "logit"))
#   
# m2 <- glm(PlasmaDetection ~ PercentAg + Julian + EnvDetection + 
#                 time_hours + MigStatus + SPEI,
#             data = birds,
#             na.action = na.omit,
#             family = binomial(link = "logit"))

# summary(m1)
# summary(m2)
   
# is random effect of site needed or helpful? yes for this model. simpler 
# model without random effect preferred though
# m1 <- glmmTMB(PlasmaDetection ~ PercentAg + Julian + EnvDetection + time_hours + 
#                MigStatus,
#              data = full,
#              na.action = na.omit,
#              family = binomial(link = "logit")) # preferred by AIC
#    
# m2 <- glmmTMB(PlasmaDetection ~ PercentAg + Julian + EnvDetection + time_hours + 
#                    MigStatus + (1|Site),
#                  data = full,
#                  na.action = na.omit,
#                  family = binomial(link = "logit")) # preferred by AIC
# 
# summary(m2) # variance is 0.2674, sd = 0.5171
#    
 # model_names <- paste0("m", 1:2)
 # models <- mget(model_names)
 # aictab(models, modnames = model_names)
  
# view individual relationships
# ggplot(birds, aes(y = time_hours, x = as.factor(PlasmaDetection))) + 
#   geom_boxplot() + my_theme

#---

# ...body condition model ----
# BCI = body condition index corrected for species and season
# values > 0 = individual is in better condition relative to its species captured
# during that season
# values < 0 = individual is in worse condition relative to its species captured
# during that season
# don't need to account for season as a predictor because it's already accounted for in BCI

#### ...final body condition index model ----
m3 <- lm(BCI ~ Biomass + PercentAg + PlasmaDetection +
           time_hours,
           data = full,
           na.action = na.omit)

# view individual relationships
# no clear relationships with BCI
 # ggplot(birds, aes(x = PercentAg, y = BCI)) + geom_point() +
 #   geom_hline(yintercept = 0, linetype = "dashed", color = "red",
 #              size = 1) + my_theme
 #  

# is random effect of site needed or helpful? no
# m1 <- glmmTMB(BCI ~ Biomass + PercentAg + PlasmaDetection +
#             time_hours + (1|Site),
#           family = gaussian(),
#           data = full,
#           na.action = na.omit)
# 
# m2 <- glmmTMB(BCI ~ Biomass + PercentAg + PlasmaDetection +
#             time_hours,
#           family = gaussian(),
#           data = full,
#           na.action = na.omit) # simpler model preferred by AIC
# 
# summary(m1) # variability ~ 0
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


#---

# ...fattening index (FI) model ----

#### ...final fattening index model ----
m4 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
          PlasmaDetection + time_hours + BCI + EnvDetection,
          family = gaussian(),
          data = full,
          na.action = na.omit)

  
# season or julian as time variable? they're about the same...season is better
# and makes more sense conceptually anyways
# m1 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Julian +
#             PlasmaDetection + time_hours + BCI + EnvDetection,
#             family = gaussian(),
#             data = full,
#             na.action = na.omit)
#   
# m2 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#             PlasmaDetection + time_hours + BCI + EnvDetection,
#             family = gaussian(),
#             data = full,
#             na.action = na.omit)
#   
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)
# 
# # plot relationships
# ggplot(data = full, aes(y = FatteningIndex, x = Julian)) + geom_point()
  

# is random effect of site needed or helpful? NO! Causes issues
# m1 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Julian +
#             PlasmaDetection + time_hours + BCI + EnvDetection,
#           family = gaussian(),
#           data = full,
#           na.action = na.omit)
# 
# m2 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Julian +
#             PlasmaDetection + time_hours + BCI + EnvDetection +
#               (1|Site),
#           family = gaussian(),
#           data = full,
#           na.action = na.omit)
# 
# summary(m2) # variance ~ 0; sd ~ 0
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

# view individual relationships
# ggplot(birds.all, aes(x = SPEI, y = FatteningIndex)) + geom_point() +
#   geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
#    geom_hline(yintercept = 0, linetype = "dashed", color = "red",
#               size = 1) + my_theme

#---

# ...environmental detection model ----
# extract one row per site to avoid pseudoreplication in analysis (n = 11 wetlands)
site_data_wetland <- full %>%
  distinct(Site, PercentAg, Season, AnnualSnowfall_in,
           SPEI, DaysSinceLastPrecipitation_5mm, EnvDetection,
           PesticideInvert_ng.g, WaterNeonicConc)
  
#  model with complete dataset with main variable of interest (PercentAg)
#### ...final environmental detection model ----
m5 <- glm(EnvDetection ~ PercentAg,
          family = binomial(link = "logit"),
         data = site_data_wetland,
         na.action = na.omit)

# can we use full wetland dataset? answer is no, not reliably across predictors
# will need to test each covariate individually to avoid overfitting
# ...AnnualSnowfall_in has no bias, relationship is not significant though anyways
# ...PercentAg has BIG bias, relationship is not significant though anyways
# ...Season has BIG bias, relationship is not significant though anyways
# ...SPEI has some bias, relationship is not significant though anyways
# ...DaysSinceLastPrecip has no bias, relationship is not significant though anyways

# m1 <- glm(EnvDetection ~ DaysSinceLastPrecipitation_5mm,
#           family = binomial(link = "logit"),
#           data = wetland,
#           na.action = na.omit)
# 
# m2 <- glm(EnvDetection ~ DaysSinceLastPrecipitation_5mm,
#           family = binomial(link = "logit"),
#           data = site_data_wetland,
#           na.action = na.omit)
# 
# summary(m1)
# summary(m2)

# view individual relationships
# no clear relationships with any covariates
# ggplot(wetland, aes(y = SPEI, x = as.factor(EnvDetection))) +
#    geom_boxplot(aes(group = EnvDetection)) + my_theme
  
  
# do concentrations or detections have more explanatory power? DETECTIONS DO
# water conc. aren't a good variable because there were no neonic detections
# m10 <- lm(PesticideInvert_ng.g ~ PercentAg,
#           data = site_data_wetland,
#           na.action = na.omit)
# 
# summary(m10)
# 
# ggplot(site_data_wetland, aes(y = PesticideInvert_ng.g, x = PercentAg)) +
#      my_theme + geom_point() + geom_smooth(method = "lm")
# 
# # worse fit
# model.conc <- psem(m1, m2, m3, m4, m10)
# summary(model.conc, conserve = TRUE)

  
#---
  
# ...pectoral muscle model ----
# m6 <- lm(Standardized.Pec ~ Biomass + MigStatus + PercentAg + Season +
#         PlasmaDetection + time_hours + EnvDetection,
#         na.action = na.omit,
#         data = full)
  
# view individual relationships
# no clear relationships
# ggplot(full, aes(x = as.factor(EnvDetection), y = Standardized.Pec)) + geom_point() +
#    geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "red",
#                size = 1) + my_theme
  
  
#---
  
# ...fat model ISSUES WITH MODEL DIAGNOSTICS AND RIGHT SKEWED DATA ----
# m7 <- lm(Fat ~ Biomass + MigStatus + PercentAg + Season +
#          PlasmaDetection + time_hours + EnvDetection,
#          na.action = na.omit,
#          data = full)
  
# view individual relationships
# higher fat in migrants
# lower fat in spring
# lower fat in birds with plasma detections (but there's a correlation with season)
# ggplot(full, aes(x = time_hours,  y = Fat)) + geom_point() +
#     geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) + my_theme
  
#---
  
# ...uric acid model  ----
# m7 <- lm(Uric ~ Biomass + MigStatus + PercentAg + Season +
#           PlasmaDetection + time_hours + BCI,
#           na.action = na.omit,
#           data = full)
#   
# removed because it forces new missing paths that don't make any sense

#------------------------------------------------------------------------------#
#                            run piecewise SEMs                             ----                        
#------------------------------------------------------------------------------# 

# FINAL SEM --------------------------------------------------------------------
# BCI as primary body condition index (BEST & SIMPLEST MODEL)                  #  
model.BCI <- psem(m1, m2, m3, m4, m5)
summary(model.BCI, conserve = TRUE)
# Fisher's C = 19.926 with P=value = 0.588 and on 22 degrees of freedom
# AIC 264.024


# Create fake data
set.seed(1)

data <- data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Create SEM using `psem`
modelList <- psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data),
  data
)

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))

# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)

summary(modelList2)

basisSet(model.BCI)
#------------------------------------------------------------------------------#

# Standardized pectoral muscle size as primary body condition index
# m4 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#             PlasmaDetection + time_hours + Standardized.Pec + EnvDetection,
#           family = gaussian(),
#           data = full,
#           na.action = na.omit)
# 
# model.Pec <- psem(m1, m2, m4, m5, m6)
# summary(model.Pec, conserve = TRUE)
# # Fisher's C = 218.88 with P-value = 0 and on 16 degrees of freedom
# # AIC 313.507
# # missing paths that don't make sense (except biomass ~ envdetection but 
# # we know apriori that relationship is not significant)
# 
  
  
# # BCI and Standardized Pectoral Muscle in same model
# m3 <- lm(BCI ~ Biomass + PercentAg + PlasmaDetection +
#            time_hours + Standardized.Pec,
#          data = full,
#          na.action = na.omit)
# 
# m4 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#             PlasmaDetection + time_hours + BCI + Standardized.Pec +
#             EnvDetection,
#           family = gaussian(),
#           data = full,
#           na.action = na.omit)
# 
# model.BCIandPec <- psem(m1, m2, m3, m4, m5, m6)
# summary(model.BCIandPec, conserve = TRUE)
# Fisher's C = 225.564 with P-value = 0 and on 22 degrees of freedom
# AIC 229.842
# lots of missing paths that don't make sense
  
  
# Model that also includes uric acid 
# model.BCIandUric <- psem(m1, m2, m3, m4, m5, m7)
# summary(model.BCIandUric, conserve = TRUE)
# lots of missing paths that don't make sense, brought in new correlations
# extremely high AIC (1306.766)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 

# biomass model (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# m2 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m2)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Julian) # good

# m3 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m3)$Biomass) # good


# m4 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m4)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$BCI) # good

# env detection model (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m5, type = "deviance") ~ fitted(m5))
abline(h = 0, col = "red") # reasonable fit

# m6 --- GOOD, no severe violations overall
# simulationOutput <- simulateResiduals(fittedModel = m6) 
# plot(simulationOutput)
# testDispersion(m6) 
# testUniformity(simulationOutput)
# testOutliers(simulationOutput) 
# testQuantiles(simulationOutput) 
# 
# plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg) # good
# plotResiduals(simulationOutput, form = model.frame(m6)$Biomass) # good
# plotResiduals(simulationOutput, form = model.frame(m6)$MigStatus) # good
# plotResiduals(simulationOutput, form = model.frame(m6)$time_hours) # good
# plotResiduals(simulationOutput, form = model.frame(m6)$Season) # some pattern
# plotResiduals(simulationOutput, form = model.frame(m6)$PlasmaDetection) # good
# plotResiduals(simulationOutput, form = model.frame(m6)$EnvDetection) # good
# 
# 
# # m7 --- GOOD, no violations
# simulationOutput <- simulateResiduals(fittedModel = m7) 
# plot(simulationOutput)
# testDispersion(m7) 
# testUniformity(simulationOutput)
# testOutliers(simulationOutput) 
# testQuantiles(simulationOutput) 
# 
# plotResiduals(simulationOutput, form = model.frame(m7)$PercentAg) # good
# plotResiduals(simulationOutput, form = model.frame(m7)$Biomass) # good
# plotResiduals(simulationOutput, form = model.frame(m7)$MigStatus) # good
# plotResiduals(simulationOutput, form = model.frame(m7)$time_hours) # good
# plotResiduals(simulationOutput, form = model.frame(m7)$Season) # some pattern
# plotResiduals(simulationOutput, form = model.frame(m7)$PlasmaDetection) # good
# plotResiduals(simulationOutput, form = model.frame(m7)$EnvDetection) # good

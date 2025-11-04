#----------------------------------------#
#           SEM Model Building           #
#          Full Dataset Overlap          #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 10/30/2025         #
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

options(tibble.print_max = Inf)

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


# use complete cases of linked variables (n = 114)
full <- full %>%
 filter(complete.cases(BCI.NoEvent, Standardized.Pec.NoEvent,
                       PlasmaDetection))

# Only include species with at least three individuals (n = 114)
full <- full %>% 
   group_by(Species) %>% 
    filter(n() >= 3) %>% 
    ungroup()

# I've decided to model response body condition metrics as correlated errors, 
# and not specify direction...otherwise I get spurious results (birds with 
# smaller pectoral muscles are refueling faster?? but birds with more mass
# are also refueling faster?)

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
  distinct(Site, Biomass, PercentAg, Season, EnvDetection)

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

# what about for envdetection? no, much worse
# m1 <- glm(Biomass ~ PercentAg, data = site_data,
#           family = Gamma(link = "log"))
# 
# m2 <- glm(Biomass ~ PercentAg + EnvDetection, data = site_data,
#           family = Gamma(link = "log"))


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

# WITH RANDOM EFFECT OF SPECIES -- causing singular boundary fit error
# variance ~ 0 
  
# glmmTMB (Failed to converge in SEM)
 # m2 <- glmmTMB(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
 #              MigStatus + (1|Species),
 #              data = full.3,
 #              na.action = na.omit,
 #              family = binomial(link = "logit"))
  
# no random effect
m2 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
                  MigStatus,
                data = full,
                na.action = na.omit,
                family = binomial(link = "logit"))
  
# glmer caused boundary (singular) fit warning

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
# model without random effect preferred though (AICc wt = 0.67)
# m1 <- glmmTMB(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
#                 MigStatus,
#               data = full,
#               na.action = na.omit,
#               family = binomial(link = "logit")) # preferred by AIC
#     
# m2 <- glmmTMB(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
#                    MigStatus + (1|Site),
#                   data = full,
#                   na.action = na.omit,
#                   family = binomial(link = "logit")) # preferred by AIC
#  
# summary(m2) # variance is 0.2165, sd = 0.4653
# #    
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
m3 <- lm(BCI.NoEvent ~ Biomass + PercentAg + PlasmaDetection +
         time_hours + Season + Standardized.Pec.NoEvent,
         data = full,
         na.action = na.omit)
  
# does including fattening index here improve model fit?
# is the relationship BCI -> FI or FI -> BCI?
# Including fattening index in this model was significantly WORSE!
# likely not FI -> BCI
# m1 <- lm(BCI.NoEvent ~ Biomass + PercentAg + PlasmaDetection +
#            time_hours + Season + Standardized.Pec.NoEvent +
#            FatteningIndex,
#            data = full,
#            na.action = na.omit)
#   
# m2 <- lm(BCI.NoEvent ~ Biomass + PercentAg + PlasmaDetection +
#            time_hours + Season + Standardized.Pec.NoEvent,
#            data = full,
#            na.action = na.omit)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

  

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
  
# WITH RANDOM EFFECT

# variance = 0.1227
m4 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
                 PlasmaDetection + time_hours + EnvDetection + (1|Species),
                data = full,
                family = gaussian(),
                 na.action = na.omit)
  
# does including pectoral muscle index improve model fit? 
# yes! but in a negative way?? 
  
# ggplot(full, aes(x = Standardized.Pec.NoEvent, y = FatteningIndex)) +
#   geom_point() + geom_hline(yintercept = 0) + geom_smooth(method = "lm")
# 
# m <- glmmTMB(FatteningIndex ~ Standardized.Pec.NoEvent + BCI.NoEvent +
#           (1|Species), data = full)

# m1 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                  PlasmaDetection + time_hours + EnvDetection + BCI.NoEvent +
#                 (1|Species),
#                 data = full,
#                 family = gaussian(),
#                  na.action = na.omit)
#   
# m2 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                  PlasmaDetection + time_hours + EnvDetection + BCI.NoEvent +
#                 (1|Species) + Standardized.Pec.NoEvent,
#                 data = full,
#                 family = gaussian(),
#                  na.action = na.omit)
# 
# m2 <- lm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                 PlasmaDetection + time_hours + EnvDetection + BCI.NoEvent +
#                 Standardized.Pec.NoEvent,
#               data = full,
#               na.action = na.omit)
  
# does including bci improve model fit? is it BCI -> FI or FI -> BCI
# definitely BCI -> FI
# AIC STRONGLY preferred model with FI ~ BCI
# m1 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                 PlasmaDetection + time_hours + EnvDetection + (1|Species) +
#                 BCI.NoEvent,
#                 data = full,
#                 family = gaussian(),
#                 na.action = na.omit)
  
# m2 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                 PlasmaDetection + time_hours + EnvDetection + (1|Species),
#               data = full,
#               family = gaussian(),
#               na.action = na.omit)
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)



# no random effect
# m4 <- lm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                 PlasmaDetection + time_hours,
#               data = full,
#               na.action = na.omit)

  
# m4 <- lmer(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                 PlasmaDetection + time_hours + (1|Species),
#                 data = full.c, REML = FALSE,
#                 na.action = na.omit)
#   
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
#  models <- mget(model_names)
#  aictab(models, modnames = model_names)
# 
# # plot relationships
# ggplot(data = full, aes(y = FatteningIndex, x = Julian)) + geom_point()
  

# is random effect of species needed or helpful? no!
# m1 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                PlasmaDetection + time_hours,
#              family = gaussian(),
#              data = full,
#              na.action = na.omit)
#   
# m2 <- glmmTMB(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
#                PlasmaDetection + time_hours +
#                  (1|Species),
#             family = gaussian(),
#              data = full,
#              na.action = na.omit)
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
m6 <- lm(Standardized.Pec.NoEvent ~ Biomass + MigStatus + PercentAg + Season +
         PlasmaDetection + time_hours + EnvDetection,
         na.action = na.omit,
         data = full)
  
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
 #          data = full.c)
 #  
# view individual relationships
# higher fat in migrants
# lower fat in spring
# lower fat in birds with plasma detections (but there's a correlation with season)
# ggplot(full, aes(x = time_hours,  y = Fat)) + geom_point() +
#     geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) + my_theme
  
#---
  
# ...uric acid model  ----
# WITH RANDOM EFFECT OF SPECIES
# removed BCI to avoid overfitting

# no random effect (causes whole model to fail to converge)
m7 <- lm(Uric ~ Biomass + MigStatus + PercentAg + Season +
                PlasmaDetection + time_hours,
              na.action = na.omit,
              data = full)

# random effect (Variance ~ 0)
# m7 <- glmmTMB(Uric ~ Biomass + MigStatus + PercentAg + Season +
#              PlasmaDetection + time_hours + (1|Species),
#             na.action = na.omit,
#             data = full.3)


#------------------------------------------------------------------------------#
#                            run piecewise SEMs                             ----                        
#------------------------------------------------------------------------------# 
# SEM
model <- psem(m1, m2, m3, m4, m5, m6, m7)
summary(model, conserve = TRUE)

# add correlated errors for nonsensical independence claims and correlated
# body condition metrics (double arrow)
model2 <- update(model, 
                 EnvDetection %~~% time_hours,
                 Biomass %~~% EnvDetection,
                 Standardized.Pec.NoEvent %~~% BCI.NoEvent,
                 FatteningIndex %~~% BCI.NoEvent)

summary(model2, conserve = TRUE)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 

# biomass model (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# m2 --- GOOD, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput) # significant
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # pattern but not significant

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # some pattern, n.s.
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m2)$time_hours) # some pattern
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # good

# m3 --- GOOD, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m3)$Biomass) # good


# m4 --- GOOD, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # significant, but all else good
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
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # some pattern but not significant

plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m6)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m6)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m6)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m6)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m6)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m6)$EnvDetection) # good

 
# # m7 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m7) 
plot(simulationOutput)
testDispersion(m7) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 
# 
plotResiduals(simulationOutput, form = model.frame(m7)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m7)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m7)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m7)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m7)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m7)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m7)$EnvDetection) # good


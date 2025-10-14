#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 10/14/2025          #
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
wetland <- read.csv("cleaned_data/wetland_data_cleaned_2025-09-30.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")

# filter to only include 2023 data
birds <- birds %>% 
  filter(Event %in% c("Fall 2023", "Spring 2023")) # 122 birds

wetland <- wetland %>% 
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


#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# model wetlands with biomass > 0 only, given tweedie is not supported (n = 66)
# must do the same for all datasets or error related to gamma is thrown
invert.pos <- subset(invert, Biomass > 0)
birds.pos <- subset(birds, Biomass > 0)
wetland.pos <- subset(wetland, Biomass > 0)
full.pos <- subset(full, Biomass > 0) # 11 wetlands total

full.om <- subset(full, Biomass < 3) # 10 wetlands total

# can we use full invert dataset? answer is no (biased estimate)
m1 <- glm(Biomass ~ PercentAg, data = invert.pos,
          family = Gamma(link = "log")) # B = -2.48

m2 <- glm(Biomass ~ PercentAg, data = full.pos,
          family = Gamma(link = "log")) # B = -0.04

# extract one row per site to avoid pseudoreplication (n = 11)
site_data <- full %>%
  distinct(Site, Biomass, PercentAg, Season)

#  log transformation needed?
m1 <- lm(Biomass ~ PercentAg, data = site_data)
m2 <- lm(Biomass ~ PercentAg, data = site_data)

summary(m1)
summary(m2)

model_names <- paste0("m", 1:2)

models <- mget(model_names)

aictab(models, modnames = model_names)




m2 <- glm(Biomass ~ PercentAg, data = site_data,
          family = Gamma(link = "log"))







# view individual relationships
# clear effects of season and % cropland cover on biomass
 ggplot(full, aes(x = PercentAg, y = Biomass)) + 
   geom_point() + my_theme

# extract standardized coefficients manually
  beta <- coef(m1)["PercentAg"]
  sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
                 trigamma(1 / summary(m1)$dispersion)) # observation-level variance
  sd_x <- sd(invert.pos$PercentAg)
  beta_std <- beta * (sd_x / sd_y)
  beta_std

#---

# ...plasma detection model ----

# saturated
m2 <- glm(PlasmaDetection ~ PercentAg + EnvDetection + SPEI + time_hours +
                MigStatus,
            data = birds.pos,
            family = binomial(link = "logit"))

# view individual relationships
# no clear relationships with plasma detection
# ggplot(birds.pos, aes(y = time_hours, x = as.factor(PlasmaDetection))) + 
#   geom_boxplot() + my_theme

#---

# ...body condition model ----

# saturated model
m3 <- lm(BCI ~ Biomass + PercentAg + SPEI + PlasmaDetection +
           time_hours,
         data = birds.pos,
         na.action = na.omit)

# view individual relationships
# no clear relationships with BCI
# ggplot(birds.pos, aes(x = PlasmaDetection, y = BCI)) + geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red",
#              size = 1) + my_theme

#---

# ...fattening index model ----

# saturated model
m4 <- lm(FatteningIndex ~ Biomass + MigStatus + PercentAg + SPEI +
           PlasmaDetection + time_hours + BCI + EnvDetection,
         data = birds.pos,
         na.action = na.omit)

# view individual relationships
# positive relationship between FI and BCI
# positive relationship between FI and time
# positive relationship between FI and SPEI
# ggplot(birds.pos, aes(x = SPEI, y = FatteningIndex)) + geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red",
#              size = 1) + my_theme

#---

# ...environmental detection ----

# saturated model
m5 <- glm(EnvDetection ~ AnnualSnowfall_in + PercentAg + 
           SPEI + DaysSinceLastPrecipitation_5mm,
          family = binomial(link = "logit"),
         data = wetland.pos,
         na.action = na.omit)

# view individual relationships
# no clear relationships
# ggplot(wetland.pos, aes(y = DaysSinceLastPrecipitation_5mm, x = EnvDetection)) +
#   geom_boxplot(aes(group = EnvDetection)) + my_theme

#---

# ...water quality and veg models ISSUES WITH MODEL FIT----

# saturated model
m6 <- lm(WaterQuality ~ PercentAg,
         data = wetland.pos,
         na.action = na.omit)

# view individual relationships
 # ggplot(wetland.pos, aes(x = PercentAg, y = WaterQuality)) +
 #   geom_point() + my_theme + geom_hline(yintercept = 0)


#---
 
# ...vegetation cover ----
m7 <- lm(PercentLocalVeg_50m ~ WaterQuality + Season + AnnualSnowfall_in, 
         data = wetland.pos)
 

model <- psem(m1, m2, m3, m4, m5, m6, m7)
summary(model, conserve = TRUE)
# print(model)

#------------------------------------------------------------------------------#
#            random invert code               ----                        
#------------------------------------------------------------------------------# 


# TEST (mixed datasets)
m1 <- lm(Biomass ~ PercentAg + Season, data = invert) # contains extra invert information

m2 <- lm(Mass ~ PercentAg + Season + Biomass, data = birds)

model <- psem(m1,m2)
summary(model)


# TEST (full dataset)

m3 <- lm(Biomass ~ PercentAg + Season, data = full)

m4 <- lm(Mass ~ PercentAg + Season + Biomass, data =  full)

m.full <- psem(m3,m4)
summary(m.full)






# m1 --- GOOD, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testZeroInflation(m1)
testUniformity(simulationOutput) 
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = model.frame(m1)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # slight pattern
plotResiduals(simulationOutput, form = model.frame(m1)$PercentLocalVeg_50m) # good
plotResiduals(simulationOutput, form = model.frame(m1)$WaterQuality) # slight pattern
plotResiduals(simulationOutput, form = model.frame(m1)$EnvDetection) # good

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


# m3 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # very slight pattern
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours) # good


# m4 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # very slight pattern
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours) # good

# m5 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m5)$DaysSinceLastPrecipitation_5mm) # good
plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m5)$time_hours) # good


# m6 --- ISSUES
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$SPEI)
plotResiduals(simulationOutput, form = model.frame(m6)$PercentLocalVeg_50m) 
plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg)

# m7 --- ISSUES
simulationOutput <- simulateResiduals(fittedModel = m7) 
plot(simulationOutput)
testDispersion(m7) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m7)$WaterQuality)



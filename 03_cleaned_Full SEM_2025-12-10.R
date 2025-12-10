#-----------------------------------------#
#           SEM Model Building            #
#       2023 Full Dataset with Inverts    #
# Created by Shelby McCahon on 12/10/2025 #
#         Modified on 12/10/2025          #
#-----------------------------------------#

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
#                        load data and organize dataset                     ----                        
#------------------------------------------------------------------------------# 

full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")

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

full <- full %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Fat.Binomial = case_when(
      Fat >= 2 ~ 1,
      Fat < 2 ~ 0,
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

# use complete cases of linked variables (n = 67)
full <- full %>%
  filter(complete.cases(BCI.NoEvent, Fat.Binomial,
                        PlasmaDetection,
                        FatteningIndex,
                        Uric))


# notes from previous code and analyses:
# originally, I reduced dataset down to species that had >3 individuals so I 
# could consider species as a random effect. However, the random effect was
# not explaining any variation and SEM would not compute marginal R2. so I
# ultimately decided to include the four birds back in and not consider
# species as a random effect.
# For all analyses, I also fit models to the full dataset (using birds from 
# all years) to maximize the use of our data, but it led to biased results. So
# ultimately I only used complete cases.
# I considered pectoral muscle in analysis originally, but I would have to 
# further reduce my dataset, and it's highly correlated with mass anyways. So 
# I did not consider it any further. I did not find any other interesting
# relationships either with the variable.

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# extract one row per site to avoid pseudoreplication (n = 9 wetlands)
site_data <- full %>%
  distinct(Site, Biomass, PercentAg, Season, EnvDetection,
           Diversity)

# model
m1 <- glm(Biomass ~ PercentAg, data = site_data,
          family = Gamma(link = "log"))

# extract standardized coefficients manually from gamma distribution
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(site_data$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.586 is standardized estimate

# notes from previous code and analyses:
# tweedie distribution would be optimal but piecewiseSEM does not support this.
# I then decided between gamma or log-transformation of biomass with a small
# constant added, but I know from the full dataset that gamma is a much better
# fit so I used gamma for consistency with other analyses.
# I wanted to consider season as well, but model was overfit and the model
# without season was a much better model.
# I also considered wetland pesticide detection but this also worsened
# model fit. Ultimately, I just included % cropland.

#---

# diversity model

m2 <- glm(Diversity ~ PercentAg, data = site_data, na.action = na.omit,
          family = poisson())

# ...plasma detection model ----

# notes from previous code and analyses:
# I did not consider SPEI because I only had one year of data. I also
# considered Julian day vs Season but season had a better AICc score. Both
# were highly correlated though anyways.

m3 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
            MigStatus,
          data = full,
          na.action = na.omit,
          family = binomial(link = "logit"))

#---

# ...body condition model ----

# notes from previous code and analyses:
# body condition is corrected for species. I did not consider season as a 
# predictor because it's accounted for in the index.
# values > 0 = individual is in better condition relative to its species captured
# during that season
# values < 0 = individual is in worse condition relative to its species captured
# during that season

m4 <- lm(BCI.NoEvent ~ Biomass + Diversity + PercentAg + PlasmaDetection +
           time_hours + Season,
         data = full)

#---

# ...fattening index (FI) model ----

m5 <- lm(FatteningIndex ~ Biomass + Diversity + MigStatus + PercentAg + Season +
           PlasmaDetection + time_hours + EnvDetection,
         data = full,
         na.action = na.omit)
#---

# ...environmental detection model ----

# notes from previous code and analyses:
# detections have more explanatory power than concentrations

m6 <- glm(EnvDetection ~ PercentAg,
          family = binomial(link = "logit"),
          data = site_data,
          na.action = na.omit)

#---

# ...uric acid model ----

m7 <- lm(Uric ~ time_hours + PercentAg + Season +
           PlasmaDetection + MigStatus + Biomass + Diversity, data = full)

# relationship looks real
# ggplot(full, aes(x = PercentAg, y = Uric)) + geom_point() + my_theme

#---

#...fat score model ----

# notes from previous code and analyses:
# I tried a lot of different transformations and distributions for fat, but the
# best model fit and diagnostics came from when I treated it is binary. 
# Splitting the variable into low fat (0/1) and high fat (2-5) made the most 
# sense and had good diagnostics (n = 42 low fat, 25 high fat)

m8 <- glm(Fat.Binomial ~ PercentAg + Season + PlasmaDetection + 
            MigStatus + Biomass + Diversity,
          data = full,
          family = binomial(link = "logit"))

#------------------------------------------------------------------------------#
#                            run piecewise SEMs                             ----                        
#------------------------------------------------------------------------------# 
# SEM
model <- psem(m1, m2, m3, m4, m5, m6, m7, m8)
summary(model, conserve = TRUE)

# add correlated errors
sem <- update(model,
              Fat.Binomial %~~% BCI.NoEvent,
              Biomass %~~% EnvDetection,
              Biomass %~~% Diversity,
              FatteningIndex %~~% BCI.NoEvent,
              Fat.Binomial %~~% FatteningIndex)

summary(sem, conserve  = TRUE)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 

# biomass model (DHARMa diagnostics unreliable with effect size of 9)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# diversity model (DHARMa diagnostics unreliable with effect size of 9)
plot(residuals(m2, type = "deviance") ~ fitted(m2))
abline(h = 0, col = "red") # reasonable fit

# m3 --- GOOD, no significant violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput) 
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # pattern but not significant

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # some pattern, n.s.
plotResiduals(simulationOutput, form = model.frame(m3)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours) # some pattern
plotResiduals(simulationOutput, form = model.frame(m3)$Season) # good

# m4 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m4)$Biomass) # good


# m5 --- GOOD, no severe violations
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # significant, but all else good
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m5)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m5)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m5)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m5)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m5)$BCI) # good

# env detection model (DHARMa diagnostics unreliable with effect size of 9)
plot(residuals(m6, type = "deviance") ~ fitted(m6))
abline(h = 0, col = "red") # reasonable fit

# m7 --- GOOD, no violations
simulationOutput <- simulateResiduals(fittedModel = m7) 
plot(simulationOutput)
testDispersion(m7) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = model.frame(m7)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m7)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m7)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m7)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m7)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m7)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m7)$EnvDetection) # good

# # m8 --- GOOD, no significant violations
simulationOutput <- simulateResiduals(fittedModel = m8) 
plot(simulationOutput)
testDispersion(m8) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m8)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m8)$Biomass) # good
plotResiduals(simulationOutput, form = model.frame(m8)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m8)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m8)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m8)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m8)$EnvDetection) # good

# diversity model (DHARMa diagnostics unreliable with effect size of 9)
plot(residuals(m8, type = "deviance") ~ fitted(m8))
abline(h = 0, col = "red") # reasonable fit


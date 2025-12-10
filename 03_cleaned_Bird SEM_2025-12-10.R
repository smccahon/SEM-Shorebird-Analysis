#-----------------------------------------#
#           SEM Model Building            #
#   Drivers of Shorebird Body Condition   #
# Created by Shelby McCahon on 12/10/2025 #
#         Modified on 12/10/2025          #
#-----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(nlme)
library(glmmTMB)
library(multcompView)
library(DHARMa)
library(car)
library(AICcmodavg)
library(statmod)
library(MuMIn)
library(MASS)

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

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
    Fat.Binomial = case_when(
      Fat >= 2 ~ 1,
      Fat < 2 ~ 0,
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
      Season == "Fall" ~ 0,
      TRUE ~ NA_real_),
    Sex = case_when(
      Sex == "F" ~ 1,
      Sex == "M" ~ 0,
      TRUE ~ NA_real_),
    MigStatus = case_when(
      MigStatus == "Migratory" ~ 1,
      MigStatus == "Resident" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

birds <- birds %>%
  filter(complete.cases(PlasmaDetection,
                        Fat, BCI.NoEvent,
                        FatteningIndex,
                        Uric))

# notes from previous code and analyses:
# I originally reduced the dataset down to species with >3 individuals so that
# I could consider species as a random effect. However, the only model that it
# would converge in is the plasma neonic detection model. The fixed effects
# on plasma detection were very similar and in the same direction, but 
# migratory status was significant only in the fixed effects model. However,
# I ultimately decided to remove species as a random effect for consistency 
# with the other bird and invert analysis in 2023 and so that I could include
# the additional 4 birds that I had to remove originally.
# I also ran all of these analyses with Julian and Season, and season performed
# better according to AIC.
# I also tested for correlations in previous scripts and found that the only
# correlated variables of interest were season and drought

# extract one row per site to avoid pseudoreplication in analysis (n = 24 wetlands)
site_data <- birds %>%
  distinct(Site, SPEI, PercentAg, Season, EnvDetection, AnnualSnowfall_in,
           DaysSinceLastPrecipitation_5mm, PrecipitationAmount_7days)

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...wetland pesticide detection model ----

# notes from previous analysis:
# I removed season to avoid overfitting the data. The model without season is 
# also preferred (AICcwt = 0.84).

m1 <- glm(EnvDetection ~ PercentAg + SPEI, 
          data = site_data, 
          family = "binomial")

# ...plasma detection model ----

# notes from previous analyses:
# julian is a better fit than season (AICc wt = 0.81), but I used season
# for consistency with other analyses. Julian has high
# correlation with Season anyways (r = 0.95). Random effect of species 
# worsened model fit (AICc wt in model without random effect = 0.59) so I 
# decided ultimately to remove it. I also considered precipitation variables
# but it did not have a significant or interesting effect.

m2 <- glm(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
                MigStatus + time_hours,
              data = birds, family = "binomial")

# ...body condition index model ----

m3 <- lm(BCI.NoEvent ~ time_hours + PlasmaDetection + SPEI + PercentAg + 
           Season,
         data = birds)


# ...fattening index model ----

m4 <- lm(FatteningIndex ~ MigStatus + Season + SPEI + 
           PercentAg + time_hours + PlasmaDetection + EnvDetection, 
         data = birds)

# ...fat model ----

# notes from previous analysis:
# I had a lot of issues at first with fat as a response variable. I tried
# a lot of different transformations and distributions and ultimately 
# decided to treat fat as binomial (low fat: 0/1 [n = 57] and 
# high fat: 2-5 [n = 27])

m5 <- glm(Fat.Binomial ~ PercentAg + Season + SPEI + PlasmaDetection + 
            MigStatus,
          data = birds,
          family = binomial(link = "logit"))

#...uric acid level model ----

m6 <- lm(Uric ~ time_hours + PercentAg + Season + SPEI +
           PlasmaDetection + MigStatus, data = birds)

#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 

# run piecewiseSEM
model <- psem(m1,m2,m3,m4,m5,m6)
summary(model, conserve = TRUE)
sem <- update(model,Season %~~% SPEI,
              FatteningIndex %~~% BCI.NoEvent,
              Fat.Binomial %~~% BCI.NoEvent,
              Fat.Binomial %~~% FatteningIndex)

summary(sem, conserve = TRUE)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# wetland pesticide detection model 
# (DHARMa diagnostics unreliable with effect size of 16)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# plasma detection model --- no severe violations
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # some pattern

# body condition index model --- all good, no severe patterns
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # significant but no apparent problems
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m3)$FatteningIndex) # good

# fattening index model --- no severe issues
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # significant
testQuantiles(simulationOutput) # great

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m4)$BCI.NoEvent) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # good

plot(residuals(m4, type = "deviance") ~ fitted(m4))
abline(h = 0, col = "red") # residuals look just fine though

# fat model --- no severe violations
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) #  some pattern
plotResiduals(simulationOutput, form = model.frame(m5)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m5)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) # some pattern

# pectoral muscle model --- all good, no patterns
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg)
plotResiduals(simulationOutput, form = model.frame(m6)$PlasmaDetection)
plotResiduals(simulationOutput, form = model.frame(m6)$time_hours)  
plotResiduals(simulationOutput, form = model.frame(m6)$SPEI)
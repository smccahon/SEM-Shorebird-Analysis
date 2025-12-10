#-----------------------------------------#
#           SEM Model Building            #
#   Drivers of Shorebird Body Condition   #
# Created by Shelby McCahon on 10/27/2025 #
#         Modified on 11/12/2025          #
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
      Fat > 0 ~ 1,
      Fat == 0 ~ 0,
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


# extract one row per site to avoid pseudoreplication in analysis (n = 24 wetlands)
site_data <- birds %>%
  distinct(Site, SPEI, PercentAg, Season, EnvDetection, AnnualSnowfall_in,
           DaysSinceLastPrecipitation_5mm, PrecipitationAmount_7days)

birds <- birds %>%
  filter(complete.cases(PlasmaDetection,
                        FatteningIndex,
                        Uric))

# Only include species with at least three individuals (n = 80)
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

#------------------------------------------------------------------------------#
#                           test for correlations                           ----                        
#------------------------------------------------------------------------------# 

# bird correlations
cor(birds$Season, birds$SPEI) # yes
cor(birds$PercentAg, birds$SPEI) # no
cor(birds$time_hours, birds$Season) # no
cor(birds$PlasmaDetection, birds$EnvDetection, use = "complete.obs") # no
cor(birds$PlasmaDetection, birds$PercentAg, use = "complete.obs") # no
cor(birds$MigStatus, birds$Season) # no
cor(birds$MigStatus, birds$PlasmaDetection, use = "complete.obs") # no
cor(birds$SPEI, birds$AnnualSnowfall_in) # no
cor(birds$Season, birds$AnnualSnowfall_in, use = "complete.obs") # no
cor(birds$Uric, birds$FatteningIndex) # no
cor(birds$Uric, birds$SPEI) # no
cor(birds$Uric, birds$EnvDetection) # no
cor(birds$Uric, birds$PercentAg) # no
cor(birds$Uric, birds$time_hours) # no
cor(birds$Uric, birds$PlasmaDetection) # no
cor(birds$Uric, birds$MigStatus) # no

# birds_complete correlations
cor(birds_complete$Season, birds_complete$SPEI) # yes
cor(birds_complete$BCI.NoEvent, birds_complete$FatteningIndex, use = "complete.obs") # no
cor(birds_complete$PercentAg, birds_complete$SPEI) # no
cor(birds_complete$time_hours, birds_complete$Season) # no
cor(birds_complete$PlasmaDetection, birds_complete$EnvDetection, use = "complete.obs") # no
cor(birds_complete$PlasmaDetection, birds_complete$PercentAg, use = "complete.obs") # no
cor(birds_complete$MigStatus, birds_complete$Season) # no
cor(birds_complete$MigStatus, birds_complete$PlasmaDetection, use = "complete.obs") # no
cor(birds_complete$SPEI, birds_complete$AnnualSnowfall_in) # no
cor(birds_complete$Season, birds_complete$AnnualSnowfall_in, use = "complete.obs") # no

# site data correlations (none of interest)
cor(site_data$PercentAg, site_data$SPEI) # no
cor(site_data$SPEI, site_data$AnnualSnowfall_in) # no
cor(site_data$PercentAg, site_data$AnnualSnowfall_in) # no

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

## ULTIMATELY DECIDED TO INCLUDE RANDOM EFFECTS of SPECIES IN PLASMA DETECTION,
## FATTENING INDEX, AND URIC ACID MODELS
## excluding random effects made some nonsensical independence claims sign.
## and the only way to exclude those is to use a completely reduced dataset
## with no missing data, which would bring sample size down to 79
## There are no random effects of species in body condition index or 
## pectoral muscle size index because size and species is already accounted
## for in the residual analysis (log(Mass) ~ log(Wing) + Species)

# ...wetland pesticide detection model ----

m1 <- glm(EnvDetection ~ PercentAg + SPEI, 
          data = site_data, 
          family = "binomial")

# ...plasma detection model ----

# WITH RANDOM EFFECT
m2 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
             MigStatus + time_hours + (1|Species),
           data = birds, family = "binomial")

# summary(m2) # random effect warranted

# ...uric acid model ----
# m3 <- glmmTMB(Uric ~ MigStatus + Season + SPEI + EnvDetection + 
#                 PercentAg + time_hours + PlasmaDetection + (1|Species), 
#          data = birds)

# outlier significant (in migrant group), remove outlier and rerun analysis
# when outlier is removed, birds captured in wetlands with neonic detections in
# environment had lower uric acid levels...?

# birds.or <- birds %>% 
#   filter(Uric < 2000)
# 
# m3 <- glmmTMB(Uric ~ MigStatus + Season + SPEI + 
#                 PercentAg + time_hours + PlasmaDetection + EnvDetection +
#                 (1|Species), 
#               data = birds.or)
# 
# summary(m3)

 m3 <- lm(Uric ~ MigStatus + Season + SPEI + 
                 PercentAg + time_hours + PlasmaDetection + EnvDetection, 
               data = birds)

# random effect explains no variation so it's removed

# ...fattening index model ----

# WITH RANDOM EFFECT
# m4 <- glmmTMB(FatteningIndex ~ MigStatus + Season + SPEI + 
#                 PercentAg + time_hours + PlasmaDetection + EnvDetection +
#                 (1|Species), 
#          data = birds)

m4 <- lm(FatteningIndex ~ MigStatus + Season + SPEI + 
               PercentAg + time_hours + PlasmaDetection + EnvDetection, 
              data = birds)
  
# plot(m4)

# random effect explains no variation so it's removed

#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 

# run piecewiseSEM
model <- psem(m1,m2,m3,m4)
summary(model, conserve = TRUE)
sem <- update(model,Season %~~% SPEI, # must include due to known correlation
              EnvDetection %~~% time_hours) # nonsensical claim
summary(sem, conserve  = TRUE)


#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# wetland pesticide detection model 
# (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# plasma detection model --- no severe issues 
# (season and SPEI are correlated, so that may cause patterns)
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput) # quantile test significant
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # sign.
plotResiduals(simulationOutput, form = model.frame(m2)$time_hours) # sign.
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # sign.

# uric acid model --- no severe problems 
# there is an outlier in a migrant, but its a real point
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # significant
testQuantiles(simulationOutput) # good

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m3)$EnvDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m3)$MigStatus) # good

# migrants have higher uric acid levels (makes sense)
ggplot(birds, aes(x = as.factor(MigStatus), y = Uric)) +
  geom_boxplot() + my_theme + 
  labs(x = "Migration Status",
       y = "Uric Acid Levels")

hist(birds$Uric)

ggplot(birds.or, aes(x = as.factor(EnvDetection), y = Uric)) +
  geom_boxplot() + my_theme + 
  labs(x = "Wetland Neonic Detection",
       y = "Uric Acid Levels")

# fattening index model --- no severe violations
# outlier test significant, some patterns but not significant
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) # outlier significant
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # cyclical but okay
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$BCI.NoEvent) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # some pattern

ggplot(birds, aes(x = as.factor(PlasmaDetection), y = FatteningIndex)) +
  geom_boxplot() + my_theme + 
  geom_hline(yintercept = 0,
             color = "red",
             linetype = "dashed",
             size = 1) +
  labs(x = "Plasma Neonic Detection",
       y = "Fattening Index")

plot(residuals(m4, type = "deviance") ~ fitted(m4))
abline(h = 0, col = "red") # residuals look just fine though



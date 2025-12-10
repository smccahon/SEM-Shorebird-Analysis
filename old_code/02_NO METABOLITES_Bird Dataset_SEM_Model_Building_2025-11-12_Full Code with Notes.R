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
      Fat >= 2 ~ 1,
      Fat <= 1 ~ 0,
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

# including pectoral muscle cuts out 17 birds
birds <- birds %>%
  filter(complete.cases(PlasmaDetection, Fat.Binomial, BCI.NoEvent))

# Only include species with at least three individuals (n = 168)
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

#------------------------------------------------------------------------------#
#                           test for correlations                           ----                        
#------------------------------------------------------------------------------# 

# 
# # EnvDetection (none)
# cor(birds$EnvDetection, birds$PlasmaDetection, use="complete.obs") 
# cor(birds$EnvDetection, birds$BCI.NoEvent, use="complete.obs")
# cor(birds$EnvDetection, birds$Fat.Binomial, use="complete.obs") 
# cor(birds$EnvDetection, birds$Standardized.Pec.NoEvent, use="complete.obs") 
# cor(birds$EnvDetection, birds$PercentAg, use="complete.obs") 
# cor(birds$EnvDetection, birds$SPEI, use="complete.obs") 
# cor(birds$EnvDetection, birds$Season, use="complete.obs") 
# cor(birds$EnvDetection, birds$time_hours, use="complete.obs") 
# cor(birds$EnvDetection, birds$MigStatus, use="complete.obs") 
# cor(birds$EnvDetection, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # PlasmaDetection (none)
# cor(birds$PlasmaDetection, birds$BCI.NoEvent, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$Fat.Binomial, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$Standardized.Pec.NoEvent, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$PercentAg, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$SPEI, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$Season, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$time_hours, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$MigStatus, use="complete.obs") 
# cor(birds$PlasmaDetection, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # BCI.NoEvent (none)
# cor(birds$BCI.NoEvent, birds$Fat.Binomial, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$Standardized.Pec.NoEvent, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$PercentAg, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$SPEI, use="complete.obs")
# cor(birds$BCI.NoEvent, birds$Season, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$time_hours, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$MigStatus, use="complete.obs") 
# cor(birds$BCI.NoEvent, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # Fat.Binomial (FAT AND SEASON!!)
# cor(birds$Fat.Binomial, birds$Standardized.Pec.NoEvent, use="complete.obs") 
# cor(birds$Fat.Binomial, birds$PercentAg, use="complete.obs")
# cor(birds$Fat.Binomial, birds$SPEI, use="complete.obs") 
# cor(birds$Fat.Binomial, birds$Season, use="complete.obs") 
# cor(birds$Fat.Binomial, birds$time_hours, use="complete.obs") 
# cor(birds$Fat.Binomial, birds$MigStatus, use="complete.obs") 
# cor(birds$Fat.Binomial, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # Standardized.Pec.NoEvent (none)
# cor(birds$Standardized.Pec.NoEvent, birds$PercentAg, use="complete.obs")
# cor(birds$Standardized.Pec.NoEvent, birds$SPEI, use="complete.obs") 
# cor(birds$Standardized.Pec.NoEvent, birds$Season, use="complete.obs") 
# cor(birds$Standardized.Pec.NoEvent, birds$time_hours, use="complete.obs") 
# cor(birds$Standardized.Pec.NoEvent, birds$MigStatus, use="complete.obs") 
# cor(birds$Standardized.Pec.NoEvent, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # PercentAg (none)
# cor(birds$PercentAg, birds$SPEI, use="complete.obs") 
# cor(birds$PercentAg, birds$Season, use="complete.obs") 
# cor(birds$PercentAg, birds$time_hours, use="complete.obs") 
# cor(birds$PercentAg, birds$MigStatus, use="complete.obs") 
# cor(birds$PercentAg, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # SPEI (season and SPEI)
# cor(birds$SPEI, birds$Season, use="complete.obs") 
# cor(birds$SPEI, birds$time_hours, use="complete.obs")
# cor(birds$SPEI, birds$MigStatus, use="complete.obs") 
# cor(birds$SPEI, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # Season (none)
# cor(birds$Season, birds$time_hours, use="complete.obs")
# cor(birds$Season, birds$MigStatus, use="complete.obs") 
# cor(birds$Season, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # time_hours (none)
# cor(birds$time_hours, birds$MigStatus, use="complete.obs") 
# cor(birds$time_hours, birds$AnnualSnowfall_in, use="complete.obs") 
# 
# # MigStatus (none)
# cor(birds$MigStatus, birds$AnnualSnowfall_in, use="complete.obs")
# 
# # site data correlations (none of interest)
# cor(site_data$PercentAg, site_data$SPEI) # no
# cor(site_data$SPEI, site_data$AnnualSnowfall_in) # no
# cor(site_data$PercentAg, site_data$AnnualSnowfall_in) # no


#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...wetland pesticide detection model ----

m1 <- glm(EnvDetection ~ PercentAg + SPEI, 
          data = site_data, 
          family = "binomial")


# ...plasma detection model ----

# WITH RANDOM EFFECT
m2 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
             MigStatus + time_hours + (1|Species),
           data = birds, family = "binomial")

# ...body condition index model ----

m3 <- lm(BCI.NoEvent ~ time_hours + PlasmaDetection + SPEI + PercentAg + 
           Season,
         data = birds)

# ...fat model ----
# heavily right skewed
# ordinal would be ideal, but not possible in piecewiseSEM
# I tried treating data as continuous with various transformations and variance
# structures but heteroscedastsicity is still highly prevalent 
# gamma does not work
# not appropriate to treat as count data
# best solution is to split fat into low fat (0-1) and high fat (2-5)
# RANDOM EFFECT NOT NEEDED

m4 <- glm(Fat.Binomial ~ PercentAg + Season + SPEI + PlasmaDetection + 
                 MigStatus + time_hours,
          data = birds,
          family = binomial(link = "logit"))


# ...pectoral muscle size model ----

# is removing the 17 birds worth it so I can include this?
# corrected for size and species -- random effect and mig. status not needed
# m5 <- lm(Standardized.Pec.NoEvent ~ time_hours + PercentAg + SPEI + 
#                Season + PlasmaDetection, 
#                data = birds)

# plot(m5) # good fit


#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 

# run piecewiseSEM
# not including fat due to model diagnostics issues
model <- psem(m1,m2,m3,m4)
summary(model, conserve = TRUE)
sem <- update(model,Season %~~% SPEI, # must include due to known correlation
              EnvDetection %~~% time_hours, # nonsensical claim
              Fat.Binomial %~~% BCI.NoEvent) # identified correlation
summary(sem, conserve  = TRUE)

# ggplot(birds, aes(x = as.factor(PlasmaDetection), y = PercentAg)) +
#   geom_boxplot() + my_theme

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# wetland pesticide detection model 
# (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# plasma detection model --- all great, no significant patterns
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$AnnualSnowfall_in) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m2)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # good

# body condition index model --- all good, no severe patterns
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # pattern
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # pattern
plotResiduals(simulationOutput, form = model.frame(m3)$FatteningIndex) # good

# fat model --- good, no violations
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # good


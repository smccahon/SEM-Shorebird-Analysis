#-----------------------------------------#
#       SEM Body Condition Comparisons    #
# Created by Shelby McCahon on 10/20/2025 #
#         Modified on 10/20/2025          #
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
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))

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
      TRUE ~ NA_real_))

full$Fat <- as.numeric(full$Fat)

#------------------------------------------------------------------------------#
#     full SEM with BCI as primary body condition index (no random effects) ----                        
#------------------------------------------------------------------------------# 

# extract one row per site to avoid pseudoreplication in analysis (n = 11 wetlands)
site_data <- full %>%
  distinct(Site, Biomass, PercentAg, Season)

#### ...final biomass model----
m1 <- glm(Biomass ~ PercentAg, data = site_data,
          family = Gamma(link = "log"))

#### final plasma detection model ----
m2 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
            MigStatus,
          data = full,
          na.action = na.omit,
          family = binomial(link = "logit"))

#### ...final body condition index model ----
m3 <- lm(BCI ~ Biomass + PercentAg + PlasmaDetection +
           time_hours,
         data = full,
         na.action = na.omit)

#### ...final fattening index model ----
m4 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
            PlasmaDetection + time_hours + BCI + EnvDetection,
          family = gaussian(),
          data = full,
          na.action = na.omit)

# extract one row per site to avoid pseudoreplication in analysis (n = 11 wetlands)
site_data_wetland <- full %>%
  distinct(Site, PercentAg, Season, AnnualSnowfall_in,
           SPEI, DaysSinceLastPrecipitation_5mm, EnvDetection)

#  model with complete dataset with main variable of interest (PercentAg)
#### ...final environmental detection model ----
m5 <- glm(EnvDetection ~ PercentAg,
          family = binomial(link = "logit"),
          data = site_data_wetland,
          na.action = na.omit)


model <- psem(m1, m2, m3, m4, m5)
summary(model, conserve = TRUE) # AIC = 263.765

#------------------------------------------------------------------------------#
#     full SEM with fat as primary body condition index (no random effects) ----                        
#------------------------------------------------------------------------------# 
#### ...final biomass model----
m1 <- glm(Biomass ~ PercentAg, data = site_data,
          family = Gamma(link = "log"))

#### final plasma detection model ----
m2 <- glm(PlasmaDetection ~ PercentAg + Season + EnvDetection + time_hours + 
            MigStatus,
          data = full,
          na.action = na.omit,
          family = binomial(link = "logit"))

#### ...final fat model ----
m3 <- lm(Fat ~ Biomass + PercentAg + PlasmaDetection +
           time_hours + Season + MigStatus,
         data = full,
         na.action = na.omit)

# plot individual relationships
ggplot(data = full, aes(x = Biomass, y = Fat)) + geom_point()

#### ...final fattening index model ----
m4 <- glm(FatteningIndex ~ Biomass + MigStatus + PercentAg + Season +
            PlasmaDetection + time_hours + BCI + EnvDetection,
          family = gaussian(),
          data = full,
          na.action = na.omit)

# extract one row per site to avoid pseudoreplication in analysis (n = 11 wetlands)
site_data_wetland <- full %>%
  distinct(Site, PercentAg, Season, AnnualSnowfall_in,
           SPEI, DaysSinceLastPrecipitation_5mm, EnvDetection)

#  model with complete dataset with main variable of interest (PercentAg)
#### ...final environmental detection model ----
m5 <- glm(EnvDetection ~ PercentAg,
          family = binomial(link = "logit"),
          data = site_data_wetland,
          na.action = na.omit)

model <- psem(m1, m2, m3, m4, m5)
summary(model, conserve = TRUE) # AIC = 750.351


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

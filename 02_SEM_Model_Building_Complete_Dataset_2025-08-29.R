#----------------------------------------#
#           SEM Model Building           #
# Created by Shelby McCahon on 8/29/2025 #
#         Modified on 8/29/2025          #
#----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)
library(multcompView)
library(cplm)
library(DHARMa)
library(corrplot)
library(car)

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-08-29.csv")

# unresolved issues or ideas to explore
# maybe do a multigroup analysis with season separated...season is correlated
# with a lot of variables (Julian and Permanence)
# what do I do with species?
# issues with models converging

data(shipley)

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
#                         variable correlations                             ----                        
#------------------------------------------------------------------------------# 

# Select only numeric columns
numeric_vars <- full %>% select(where(is.numeric))

# Calculate correlation matrix (Pearson by default)
cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")

# Set a correlation threshold
threshold <- 0.60

# Get upper triangle of correlation matrix without diagonal
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA

# Find pairs with absolute correlation >= threshold
strong_corr <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)

# Format output nicely
result <- data.frame(
  Var1 = rownames(cor_matrix)[strong_corr[, 1]],
  Var2 = colnames(cor_matrix)[strong_corr[, 2]],
  Correlation = cor_matrix[strong_corr]
)

print(result)


# dominant crop type and % cropland cover
summary(aov(data = full, PercentAg ~ DominantCrop)) # p < 0.001

# season and permanence
summary(aov(data = full, Season ~ Permanence)) # p < 0.001

# season and detection
summary(aov(data = full, Season ~ EnvDetection)) # not sign. diff.

# season and percentag
summary(aov(data = full, Season ~ PercentAg)) # p < 0.001


#------------------------------------------------------------------------------#
#              fit individual models (structural equations)                 ----                        
#------------------------------------------------------------------------------# 

full.1 <- full %>% 
  mutate(across(c(PercentAg, DaysSinceLastPrecipitation_5mm), 
                ~ as.numeric(scale(.))))

m1 <- glmer(EnvDetection ~ PercentAg +
            DaysSinceLastPrecipitation_5mm + (1 | Site), 
          family = "binomial", data = full.1)

#---

full.2 <- full %>% 
  mutate(across(c(seconds_since_midnight, PercentAg, Julian), 
                ~ as.numeric(scale(.))))

m2 <- glmer(PlasmaDetection ~ EnvDetection + seconds_since_midnight + 
                PercentAg + Julian + (1 | Site), family = "binomial", 
              data = full.2)

#---

full.3 <- full %>% 
  mutate(across(c(PercentLocalVeg_50m, PercentAg), 
                ~ as.numeric(scale(.))))

m3 <- lmer(Biomass ~ PercentAg + EnvDetection +
            PercentLocalVeg_50m + WaterQuality + (1 | Site),
            data = full.3)

#---

birds.4 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, Julian, SPEI), 
                ~ as.numeric(scale(.))))

m4 <- lmer(BodyCondition ~ seconds_since_midnight + Biomass + Julian + SPEI +
             PlasmaDetection + (1 | Site), 
           data = birds.4)

#---

invert.5 <- invert %>% 
  mutate(across(c(PercentAg, Julian, AnnualSnowfall_in), 
                ~ as.numeric(scale(.))))

m5 <- lm(PercentLocalVeg_50m ~ PercentAg + Julian + Permanence +
           AnnualSnowfall_in, 
         data = invert.5)

#---

birds.6 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, SPEI, Julian),
                ~ as.numeric(scale(.))))


m6 <- lmer(FatteningIndex ~ seconds_since_midnight + Biomass + PlasmaDetection +
             SPEI + Julian + (1 | Site), data = birds.6)

#---

wetland.7 <- wetland %>% 
  mutate(across(c(AnnualSnowfall_in, Julian),
                ~ as.numeric(scale(.))))

m7 <- lm(SPEI ~ AnnualSnowfall_in + Julian, data = wetland.7)


#---

wetland.8 <- wetland %>% 
  mutate(across(c(AnnualSnowfall_in, Julian, SPEI),
                ~ as.numeric(scale(.))))

m8 <- glm(Permanence ~ AnnualSnowfall_in + Julian + SPEI, data = wetland.8,
          family = "poisson")

#---

invert.9 <- invert %>% 
  mutate(across(c(PercentAg, PercentLocalVeg_50m),
                ~ as.numeric(scale(.))))

m9 <- lm(WaterQuality ~ PercentLocalVeg_50m + PercentAg, data = invert.9)

#---

model <- psem(m1, m2, m3, m4, m5, m6, m7, m8, m9)
summary(model)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# ...environmental detection model (ISSUES) ----
DHARM <- simulateResiduals(fittedModel = m1)
plot(DHARM)
plotResiduals(DHARM, full$DaysSinceLastPrecipitation_5mm)
testDispersion(DHARM) # no overdispersion
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation

# ...plasma detection model (GOOD) ----
DHARM <- simulateResiduals(fittedModel = m2)
plot(DHARM)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation


# ...biomass model (GOOD) ----
plot(residuals(m3), main = "Deviance Residuals", ylab = "Residuals", 
     xlab = "Index") # no clear pattern

#...body condition model ----
DHARM <- simulateResiduals(fittedModel = m4)
plot(DHARM) # combined adjusted quantile test significant
plotResiduals(DHARM, form = birds.4$Julian)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation






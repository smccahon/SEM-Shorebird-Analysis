#----------------------------------------#
#           SEM Model Building           #
#          Invertebrate Dataset          #
#  Created by Shelby McCahon on 10/31/25 #
#         Modified on 11/05/2025         #
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

invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))

options(tibble.print_max = Inf)
options(digits = 3)

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

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
      Permanence == "Temporary" ~ 1,
      Permanence == "Seasonal" ~ 2,
      Permanence == "Semipermanent" ~ 3,
      Permanence == "Permanent" ~ 4,
      TRUE ~ NA_real_))

#------------------------------------------------------------------------------#
#                         identify correlations                             ----                        
#------------------------------------------------------------------------------# 

cor_mat <- cor(invert[sapply(invert, is.numeric)], 
               use = "pairwise.complete.obs")

# Convert matrix to a data frame of pairs
cor_df <- as.data.frame(as.table(cor_mat))

# Rename columns for clarity
names(cor_df) <- c("Var1", "Var2", "Correlation")

# Convert factors to characters
cor_df$Var1 <- as.character(cor_df$Var1)
cor_df$Var2 <- as.character(cor_df$Var2)

# Keep only unique pairs (upper triangle)
cor_df <- cor_df[cor_df$Var1 < cor_df$Var2, ]

# Filter by threshold
high_corr <- subset(cor_df, abs(Correlation) > 0.6)

# Sort by correlation strength
high_corr <- high_corr[order(-abs(high_corr$Correlation)), ]

high_corr

#                      Var1                     Var2 Correlation
# 552    Conductivity_uS.cm             Salinity_ppt       0.999
# 287          Salinity_ppt             WaterQuality       0.998
# 288    Conductivity_uS.cm             WaterQuality       0.997
# 286              TDS_mg.L             WaterQuality       0.991
# 527          Salinity_ppt                 TDS_mg.L       0.981
# 528    Conductivity_uS.cm                 TDS_mg.L       0.979
# 9                  Julian                   Season      -0.970
# 333          EnvDetection InvertPesticideDetection       0.854
# 227              Buffered                PercentAg      -0.748
# 131              Buffered    NearestCropDistance_m       0.731
# 222 NearestCropDistance_m                PercentAg      -0.725
# 15             Permanence                   Season      -0.630
# 13              Diversity                   Season      -0.626
# 345                Julian               Permanence       0.615

#------------------------------------------------------------------------------#
#                 pick which water quality metric to use                    ----                        
#------------------------------------------------------------------------------# 

# invert <- invert %>% 
#   mutate(Biomass.m = Biomass + 0.001) # n = 77

# # which water quality metric should I use? I will have to log transform
# # response so water quality index (PCA) is not going to work
# m1 <- glm(Biomass.m ~ Conductivity_uS.cm, invert,
#           family = Gamma(link = "log"))
# m2 <- glm(Biomass.m ~ Salinity_ppt, invert,
#           family = Gamma(link = "log")) # the best, but all pretty much the same
# m3 <- glm(Biomass.m ~ TDS_mg.L, invert,
#           family = Gamma(link = "log"))
# 
# # what about pH? important variable
# m4 <- glm(Biomass.m ~ pH_probe, invert,
#           family = Gamma(link = "log"))
# 
# model_names <- paste0("m", 1:4)
# models <- mget(model_names)
# aictab(models, modnames = model_names)
# 
# m1 <- glm(Diversity ~ Conductivity_uS.cm, invert,
#           family = poisson())
# m2 <- glm(Diversity ~ Salinity_ppt, invert,
#           family = poisson())
# m3 <- glm(Diversity ~ TDS_mg.L, invert,
#           family = poisson()) # the best, but all within 2 delta AICc
# 
# # what about pH? important variable
# m4 <- glm(Diversity ~ pH_probe, invert,
#           family = poisson())
# 
# model_names <- paste0("m", 1:4)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


# biological relevance: choose conductivity (directly tied to ion concentration
# and osmotic stress on inverts; more commonly seen in literature)

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# add a small constant to allow gamma distribution use
# I chose the smallest constant I could without causing model to fail
# piecewiseSEM does not support tweedie distribution which is ideal dist.
# next best is gamma
invert <- invert %>% 
  mutate(Biomass.m = Biomass + 0.001) # n = 77

invert <- invert %>% 
  mutate(LogConductivity = log(Conductivity_uS.cm))

# need to simplify model to avoid model convergence issues
# removed vegetation cover and buffer presence 
m1 <- glm(Biomass.m ~ PercentAg + Season + EnvDetection + 
            LogConductivity + Permanence + pH_probe, 
          data = invert,
          family = Gamma(link = "log"))

# extract standardized coefficients manually from gamma distribution
# % surrounding cropland
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.239 is standardized estimate

# season
beta <- coef(m1)["Season"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std #-0.103  is standardized estimate

# environmental detection
beta <- coef(m1)["EnvDetection"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # 0.0396 is standardized estimate

# conductivity
beta <- coef(m1)["LogConductivity"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.0316 is standardized estimate

# pH
beta <- coef(m1)["pH_probe"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.0146 is standardized estimate

# permanence
beta <- coef(m1)["Permanence"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # 0.045 is standardized estimate



# ...diversity model ----
m2 <- glm(Diversity ~ PercentAg + Season + LogConductivity + 
            EnvDetection + Permanence + pH_probe, 
          data = invert,
          family = poisson())

# ...pesticide detection model ----
m3 <- glm(EnvDetection ~ PercentAg + Season + Buffered + Permanence,
          data = invert, family = binomial())

# should I add precipitation/snowfall?
# annual snowfall not significant, no spurious results (AIC = 595) -> 593 without
# days since last precipitation only has an effect on pH (no subsequent effects)
# precipitation amount results in lower pH and lower detections...?
# m3 <- glm(EnvDetection ~ PercentAg + Season + Buffered + Permanence +
#           DaysSinceLastPrecipitation_5mm,
#           data = invert, family = binomial())

# ...water quality index model ----
# must log transform to account for heteroscedasticity in response
m4 <- lm(LogConductivity ~ PercentAg + Buffered + Season + Permanence, 
         data = invert)

# view individual relationships
# # no clear relationship
# ggplot(invert, aes(x = PercentAg, y = WaterQuality)) +
#   geom_point() + geom_hline(yintercept = 0) + my_theme
# 
# # no clear relationship
# ggplot(invert, aes(x = as.factor(Buffered), y = WaterQuality)) +
#   geom_boxplot() + geom_hline(yintercept = 0) + my_theme
# 
# # water quality much lower in the spring (23 in spring, 54 in fall)
# ggplot(invert, aes(x = as.factor(Season), y = WaterQuality)) +
#   geom_boxplot() + geom_hline(yintercept = 0) + my_theme
# 
# # higher in semipermanent and permanent wetlands
# ggplot(invert, aes(x = as.factor(Permanence), y = WaterQuality)) +
#   geom_point() + geom_hline(yintercept = 0) + my_theme


# ...pH model ----
m5 <- lm(pH_probe ~ PercentAg + Buffered + Season + Permanence +
           PrecipitationAmount_7days, 
         data = invert)

# good
# plot(m5)


# run piecewiseSEM
model <- psem(m1,m2,m3,m4,m5)
summary(model, conserve = TRUE)

model2 <- update(model, 
                 Permanence %~~% Season,
                 PercentAg %~~% Buffered)
summary(model2, conserve = TRUE)

ggplot(invert, aes(x = as.factor(EnvDetection), y = PrecipitationAmount_7days)) +
  geom_boxplot() + my_theme


#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 


# m1 ---  good, no severe violations
# significant pattern in quantiles
simulationOutput <- simulateResiduals(fittedModel = m1) # quantile test sign.
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m1)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m1)$Season) # significant (maybe due to correlation)
plotResiduals(simulationOutput, form = model.frame(m1)$Permanence) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Conductivity_uS.cm) # significant
plotResiduals(simulationOutput, form = model.frame(m1)$pH_probe) # good

# m2 --- all good, no violations at all
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m2)$WaterNeonic.m) 
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection) 
plotResiduals(simulationOutput, form = model.frame(m2)$Season) 
plotResiduals(simulationOutput, form = model.frame(m2)$Permanence) 
plotResiduals(simulationOutput, form = model.frame(m2)$WaterQuality) 

# m3 --- all good, no violations
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m3)$Season) 
plotResiduals(simulationOutput, form = model.frame(m3)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m3)$Permanence) # slight issue

# m4 ---  good, no severe violations with log transformation of response
# patterns, but n.s.
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # pattern but n.s.

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m4)$Season) 
plotResiduals(simulationOutput, form = model.frame(m4)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m4)$Permanence) # pattern but n.s.

# m5 --- good, no severe violations
# patterns in residuals but not significant except for permanence
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m5)$Season) 
plotResiduals(simulationOutput, form = model.frame(m5)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m5)$Permanence) # sign.


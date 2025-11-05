#----------------------------------------#
#           SEM Model Building           #
#      Spring Invertebrate Dataset       #
#  Created by Shelby McCahon on 11/04/25 #
#         Modified on 11/04/2025         #
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

# subset to spring (n = 35)
invert <- subset(invert, Season == "Spring")

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

#                         Var1                     Var2 Correlation
# 286                 TDS_mg.L             WaterQuality       1.000
# 288       Conductivity_uS.cm             WaterQuality       1.000
# 287             Salinity_ppt             WaterQuality       1.000
# 528       Conductivity_uS.cm                 TDS_mg.L       1.000
# 527             Salinity_ppt                 TDS_mg.L       0.999
# 552       Conductivity_uS.cm             Salinity_ppt       0.999
# 131                 Buffered    NearestCropDistance_m       0.851
# 117             EnvDetection     WaterNeonicDetection       0.793
# 96        Conductivity_uS.cm     PesticideInvert_ng.g       0.779
# 268     PesticideInvert_ng.g             WaterQuality       0.772
# 532     PesticideInvert_ng.g             Salinity_ppt       0.772
# 508     PesticideInvert_ng.g                 TDS_mg.L       0.766
# 325                Diversity InvertPesticideDetection       0.758
# 222    NearestCropDistance_m                PercentAg      -0.739
# 227                 Buffered                PercentAg      -0.733
# 291                  Biomass                Diversity       0.717
# 350 InvertPesticideDetection               Permanence       0.683
# 110 InvertPesticideDetection     WaterNeonicDetection      -0.667

#------------------------------------------------------------------------------#
#                      which variables can go in model                      ----                        
#------------------------------------------------------------------------------# 

table(invert$Buffered) # 27 not buffered, 8 buffered -> just use % ag
table(invert$WaterNeonicDetection) # 14 with detections, 21 without detections
table(invert$InvertPesticideDetection) # 4 with detections, 6 without

table(invert$EnvDetection) # 18 with detections 
# this includes the 4 wetlands with invert pesticide detections 
# (that didn't have neonics in the wetland)

table(invert$Permanence) # very uneven, not worth including

# variables of interest and for simplicity:
# % cropland, biomass, diversity, water neonic conc

# SHOULD I INCLUDE DETECTION? IF I DO, IT's  SIGN. POSITIVE EFFECT WHICH
# MAKES NO SENSE (btoh for envdetection and waterneonicdetection)

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# add a small constant to allow gamma distribution use
# I chose the smallest constant I could without causing model to fail
# piecewiseSEM does not support tweedie distribution which is ideal dist.
# next best is gamma
invert <- invert %>% 
  mutate(Biomass.m = Biomass + 0.0001) # n = 77

# add a small constant so we can log transform
invert <- invert %>% 
  mutate(LogWaterNeonic.m = log(NeonicWater_ng.L + 0.000001))

invert <- invert %>% 
  mutate(WaterNeonic.m = NeonicWater_ng.L + 0.000001)


# need to simplify model to avoid model convergence issues
# removed vegetation cover, buffer presence, and water quality variables
m1 <- glm(Biomass.m ~ PercentAg + 
            LogWaterNeonic.m,
          data = invert,
          family = Gamma(link = "log"))

# view individual relationships
# ggplot(invert, aes(x = PercentAg, y = Biomass.m)) +
#   geom_point() + my_theme
# 
# ggplot(invert, aes(x = PercentAg, y = log(Biomass.m))) +
#   geom_point() + my_theme
# 
# ggplot(invert, aes(x = as.factor(WaterNeonicDetection), 
#                    y = Biomass.m)) +
#   geom_boxplot() + my_theme

# heavily right skewed still
# hist(invert$Biomass.m)

# extract standardized coefficients manually from gamma distribution
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.277 is standardized estimate

# ...diversity model ----
m2 <- glm(Diversity ~ PercentAg + 
            LogWaterNeonic.m,
          data = invert,
          family = poisson())

# ...water neonic conc. model ----
m3 <- lm(LogWaterNeonic.m ~ PercentAg,
         data = invert)

# run piecewiseSEM
model <- psem(m1,m2,m3)
summary(model, conserve = TRUE)

model2 <- update(model, 
                 Diversity %~~% Biomass.m)
summary(model2, conserve = TRUE)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 


# m1 ---  good, no violations
simulationOutput <- simulateResiduals(fittedModel = m1)
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m1)$LogWaterNeonic.m) 

# m2 --- overdispersed!!
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

# m5 --- ISSUES; issues with all variables
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput) # ISSUE
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m5)$Season) 
plotResiduals(simulationOutput, form = model.frame(m5)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m5)$Permanence) 


# m6 ---  good, no severe violations
# pattern but not significant
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg) # issue but n.s.
plotResiduals(simulationOutput, form = model.frame(m6)$Season) 
plotResiduals(simulationOutput, form = model.frame(m6)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m6)$Permanence) # issue, sign likely due to correlation

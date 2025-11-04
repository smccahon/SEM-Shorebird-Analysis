#----------------------------------------#
#           SEM Model Building           #
#          Invertebrate Dataset          #
#  Created by Shelby McCahon on 10/31/25 #
#         Modified on 11/03/2025         #
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

options(digits = 3)
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
# 9                  Julian                   Season      -0.970
# 294          EnvDetection InvertPesticideDetection       0.854
# 200              Buffered                PercentAg      -0.748
# 116              Buffered    NearestCropDistance_m       0.731
# 195 NearestCropDistance_m                PercentAg      -0.725
# 15             Permanence                   Season      -0.630
# 13              Diversity                   Season      -0.626
# 303                Julian               Permanence       0.615

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# lots of convergence issues...

# ...biomass model ----

# add a small constant to allow gamma distribution use
invert <- invert %>% 
  mutate(Biomass.m = Biomass + 0.001) # n = 77

invert <- invert %>% 
  mutate(WaterNeonic.m = NeonicWater_ng.L + 0.000000001)

# need to simplify model to avoid model convergence issues
# maybe remove vegetation and buffer presence
m1 <- glm(Biomass.m ~ PercentAg + Season + EnvDetection + Permanence + 
            log(WaterNeonic.m) + WaterQuality, 
          data = invert,
          family = Gamma(link = "log"))

# m1 <- glm(Biomass.m ~ log(WaterNeonic.m), 
#           data = invert,
#           family = Gamma(link = "log"))

summary(m1)

# issue with this variable
# ggplot(invert, aes(x = NeonicWater_ng.L, y = Biomass.m)) +
#   geom_point()

# extract standardized coefficients manually from gamma distribution
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.243 is standardized estimate

# ...diversity model ----
m2 <- glm(Diversity ~ PercentAg + Season + WaterQuality + 
            EnvDetection + Permanence + log(WaterNeonic.m), 
          data = invert,
          family = poisson())

summary(m2)

# ...pesticide detection model ----
m3 <- glm(EnvDetection ~ PercentAg + Season + Buffered + Permanence,
          data = invert, family = binomial())

summary(m3)


# ...water quality index model ----
m4 <- lm(WaterQuality ~ PercentAg + Buffered + Season +
           Permanence, data = invert)


# ...water neonic conc. model ----
m5 <- lm(log(WaterNeonic.m) ~ PercentAg + Season + Buffered + Permanence,
          data = invert)

# run piecewiseSEM
model <- psem(m1,m2,m3,m4,m5)
summary(model, conserve = TRUE)

model2 <- update(model, 
                 Permanence %~~% Season,
                 PercentAg %~~% Buffered,
                 log(WaterNeonic.m) %~~% EnvDetection)
summary(model2, conserve = TRUE)

#------------------------------------------------------------------------------#
#                             model diagnostics                             ----                        
#------------------------------------------------------------------------------# 


# m1 --- ISSUES
simulationOutput <- simulateResiduals(fittedModel = m1) # hmm some issues
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # significant 

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m1)$NeonicWater_ng.L)  # ISSUE
plotResiduals(simulationOutput, form = model.frame(m1)$EnvDetection)  # good
plotResiduals(simulationOutput, form = model.frame(m1)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m1)$Permanence) # good
plotResiduals(simulationOutput, form = model.frame(m1)$WaterQuality) # ISSUE

# m2 --- all good, no violations
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

# m4 --- ISSUES; issues with all variables
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput) # ISSUE
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) 
plotResiduals(simulationOutput, form = model.frame(m4)$Season) 
plotResiduals(simulationOutput, form = model.frame(m4)$Buffered) 
plotResiduals(simulationOutput, form = model.frame(m4)$Permanence) 

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
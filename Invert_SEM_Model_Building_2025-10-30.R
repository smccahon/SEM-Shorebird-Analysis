#----------------------------------------#
#           SEM Model Building           #
#          Invertebrate Dataset          #
#  Created by Shelby McCahon on 10/31/25 #
#         Modified on 10/31/2025         #
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
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# lots of convergence issues...

# ...biomass model ----

# model wetlands with biomass > 0 only, given tweedie is not supported (n = 66)
invert.pos <- subset(invert, Biomass > 0)

m1 <- glm(Biomass ~ PercentAg + Season + EnvDetection, data = invert.pos,
          family = Gamma(link = "log"))

summary(m1)

# extract standardized coefficients manually from gamma distribution
beta <- coef(m1)["PercentAg"]
sd_y <- sqrt(var(predict(m1, type = "link")) + # variance (y)
               trigamma(1 / summary(m1)$dispersion)) # observation-level variance
sd_x <- sd(invert.pos$PercentAg)
beta_std <- beta * (sd_x / sd_y)
beta_std # -0.374 is standardized estimate

# ...diversity model ----
m2 <- glm(Diversity ~ PercentAg + Season + WaterQuality +
            PercentLocalVeg_50m + EnvDetection, data = invert.pos,
          family = poisson())

summary(m2)

# ...pesticide detection model ----
m3 <- glm(EnvDetection ~ PercentAg + Season + Buffered + Permanence,
          data = invert.pos, family = binomial())

summary(m3)


# ...water quality index model ----
m4 <- lm(WaterQuality ~ PercentAg + PercentLocalVeg_50m + Buffered + Season +
           Permanence, data = invert.pos)

model <- psem(m1,m2,m3,m4)
summary(model, conserve = TRUE)






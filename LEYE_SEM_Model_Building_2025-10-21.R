#-----------------------------------------#
#      SEM Model Building for LEYE        #
# Created by Shelby McCahon on 10/21/2025 #
#         Modified on 10/21/2025          #
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

leye <- read.csv("cleaned_data/leye_fulldata_cleaned_2025-10-21.csv")
leye.2023 <- read.csv("cleaned_data/leye_2023data_cleaned_2025-10-21.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("cleaned_data/wetland_data_cleaned_2025-09-30.csv")

# filter to only include 2023 data
wetland.2023 <- wetland %>% 
  filter(Year == "2023") # 79 wetland surveys surveys

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

leye.2023 <- leye.2023 %>% 
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
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

leye <- leye %>% 
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
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

wetland.2023 <- wetland.2023 %>% 
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
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...biomass model ----

# extract one row per site to avoid pseudoreplication in analysis (n = 6 wetlands)
site_data <- leye.2023 %>%
  distinct(Site, Biomass, PercentAg)

#### ...final biomass model----
m1 <- lm(Biomass ~ PercentAg, data = site_data)

# # gamma model did not converge...
# m1 <- lm(Biomass ~ PercentAg, data = site_data)
# m2 <- lm(log(Biomass) ~ PercentAg, data = site_data)
# 
# summary(m1)
# 
# # explore individual relationships
# ggplot(site_data, aes(x = PercentAg, y = log(Biomass))) + geom_point() +
#   my_theme
# 
# ggplot(site_data, aes(x = PercentAg, y = Biomass)) + geom_point() +
#   my_theme
# 
# # examine residuals -- not much difference, choose simpler raw model
# par(mfrow=c(2,2))
# plot(m1, main = "Raw Model")
# plot(m2, main = "Log Model")

#---

# ...plasma detection model ----

# model with complete dataset -- SPEI removed because only 1 year of data
#### final plasma detection model ----
m2 <- glm(PlasmaDetection ~ PercentAg + EnvDetection + time_hours,
          data = leye.2023,
          na.action = na.omit,
          family = binomial(link = "logit"))

# summary(m2)

#---

# ...body condition model ----
#### ...final body condition index model ----
m3 <- lm(BCI ~ Biomass + PercentAg + PlasmaDetection + Julian,
         data = leye.2023,
         na.action = na.omit)

# summary(m3)

# does sex and age matter? have to test separately due to low sample size
# answer is no
# m3 <- lm(BCI ~ Sex + Age,
#          data = leye.2023,
#          na.action = na.omit)


#---

# ...environmental detection model ----
# extract one row per site to avoid pseudoreplication in analysis (n = 6 wetlands)
site_data_wetland <- leye.2023 %>%
  distinct(PercentAg, AnnualSnowfall_in,
           DaysSinceLastPrecipitation_5mm, EnvDetection)

#  model with complete dataset with main variable of interest (PercentAg)
#### ...final environmental detection model ----
m4 <- glm(EnvDetection ~ PercentAg,
          family = binomial(link = "logit"),
          data = site_data_wetland,
          na.action = na.omit)
# 
#---

# ...fattening index (FI) model ----

#### ...final fattening index model ----
m5 <- lm(FatteningIndex ~ time_hours + log(Biomass) + PercentAg + BCI,
        data = leye.2023)

# cyclical transformation? sin and cos not significant and psem doesn't support
# m <- lm(FatteningIndex ~ time_hours +
#            sin(2 * pi * time_hours / 24) +
#            cos(2 * pi * time_hours / 24),
#          data = leye.2023)
# 
# d <- expand.grid(
#   PercentAg = mean(leye.2023$PercentAg),
#   time_hours = seq(0, 24, length.out = 1000)
# )
# 
# predictions <- predict(m, newdata = d)
# 
# d$yhat <- predictions
# 
# d$yhat <- predict(m, newdata = d)
# 
# d$se <- predict(m, newdata = d, se.fit = TRUE)$se.fit
# d$lower <- d$yhat - 2*d$se
# d$upper <- d$yhat + 2*d$se
# 
# ggplot(d, aes(x = time_hours, y = yhat)) +
#   geom_line(size = 1) +  
#   theme_classic() +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#   labs(x = "Capture Time (24:00)", 
#        y = "Lesser Yellowlegs Fattening Index") +
#   theme(axis.title.x = element_text(size = 21,
#                                     margin = margin(t = 12)),
#         axis.title.y = element_text(size = 21,
#                                     margin = margin(r = 12)),
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size = 18),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 15),
#         legend.position = "none") +
#   geom_hline(yintercept = 0, linetype = "twodash", color = "red",
#              size = 1) +
#   scale_color_viridis_d(alpha = 1, begin = 1, end = 0) +
#   geom_point(data = leye.2023, aes(x = time_hours, y = FatteningIndex))



# biomass outlier was VERY influential...log transformed to satisfy assumptions
# all had similar AIC results but log was the first
# m1 <- lm(FatteningIndex ~ log(Biomass),
#           data = leye.2023,
#           na.action = na.omit)
# 
# m2 <- lm(FatteningIndex ~ Biomass,
#          data = leye.2023,
#          na.action = na.omit)
# 
# m3 <- lm(FatteningIndex ~ sqrt(Biomass),
#          data = leye.2023,
#          na.action = na.omit)
# 
# model_names <- paste0("m", 1:3)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


# # examine relationships
# # negative effect of % cropland on LEYE fattening index
# ggplot(leye.2023, aes(x = PercentAg, y = FatteningIndex)) + geom_point() +
# my_theme
# 
# # lower fattening index in wetlands with higher invert biomass
# ggplot(leye.2023, aes(x = Biomass, y = FatteningIndex)) + geom_point() +
#    my_theme + geom_hline(yintercept = 0, color = "red")

# 
# ggplot(leye.2023, aes(x = as.factor(PlasmaDetection), y = FatteningIndex)) + 
#   geom_boxplot() +
#   my_theme

# ggplot(leye.2023, aes(x = time_hours, y = FatteningIndex)) + geom_point()

# what if i removed outlier? effect is still trending negative but no longer sign.
# leye.2023.or <- subset(leye.2023, Biomass < 3)
# m <- lm(FatteningIndex ~ Biomass + PercentAg + time_hours + 
#            PlasmaDetection,
#          data = leye.2023.or,
#          na.action = na.omit)
# 
# ggplot(leye.2023.or, aes(x = Biomass, y = FatteningIndex)) + geom_point()
# 
# summary(m)
# confint(m)

#------------------------------------------------------------------------------#
#                            run piecewise SEMs                             ----                        
#------------------------------------------------------------------------------# 

model <- psem(m1,m2,m3,m4,m5)
summary(model, conserve = TRUE)

#------------------------------------------------------------------------------#
#                            model diagnostics                              ----                        
#------------------------------------------------------------------------------# 

# m5 ISSUES ----
simulationOutput <- simulateResiduals(fittedModel = m) 
plot(simulationOutput)
testDispersion(m) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m)$PercentAg) # pattern
plotResiduals(simulationOutput, form = model.frame(m)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m5)$time_hours) # pattern
plotResiduals(simulationOutput, form = model.frame(m5)$Biomass) # pattern
plotResiduals(simulationOutput, form = model.frame(m5)$BCI) # pattern

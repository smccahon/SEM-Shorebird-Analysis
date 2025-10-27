#-----------------------------------------#
#           SEM Model Building            #
#   Drivers of Shorebird Body Condition   #
# Created by Shelby McCahon on 10/27/2025 #
#         Modified on 10/27/2025          #
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
library(statmod)
library(MuMIn)

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
  distinct(Site, SPEI, PercentAg, Season, EnvDetection)

# no missing data in body condition:
birds.c <- birds %>%
  filter(complete.cases(BCI, Fat, PecSizeBest, PlasmaDetection))


#------------------------------------------------------------------------------#
#         fit individual models to full dataset (structural equations)      ----                        
#------------------------------------------------------------------------------# 

# ...wetland pesticide detection model ----

m1 <- glm(EnvDetection ~ PercentAg + SPEI + Season, data = site_data, 
          family = "binomial")

# view relationships
# some relationship apparent
# ggplot(site_data, aes(y = PercentAg, x = as.factor(EnvDetection))) +
#   geom_boxplot() + my_theme + geom_point()
# 
# # no relationship apparent
# ggplot(site_data, aes(y = SPEI, x = as.factor(EnvDetection))) +
#   geom_boxplot() + my_theme + geom_point()


# ...plasma detection model ----
# julian is a better fit than season (AICc wt = 0.81), but I used season
# for consistency with wetland pesticide detection model. Julian has high
# correlation with Season anyways (r = 0.95)
m2 <- glm(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
            MigStatus,
          data = birds.c, family = "binomial")

# ...body condition index model ----
# site as a random effect caused model to fail to converge
# model explains almost nothing so you get an error (but model is still valid)
m3 <- lm(BCI ~ time_hours + PlasmaDetection + SPEI + PercentAg,
         data = birds.c)


# ...fattening index model ----
# species performs better than migratory status alone (AICc wt = 0.96)
# species also performs better than mig. status + (1|Species) (wt = 0.80)
# issue is, I can't include species as a factor in SEM
# best I can do is species as a random effect and migratory status as fixed
m4 <- glmmTMB(FatteningIndex ~ BCI + MigStatus + Season + SPEI + PercentAg +
           time_hours + PlasmaDetection + (1|Species), data = birds.c)


# view Fattening Index ~ Plasma Detection
# birds %>%
#   filter(!is.na(PlasmaDetection) & !is.na(FatteningIndex)) %>%
#   ggplot(aes(x = as.factor(PlasmaDetection), y = FatteningIndex)) +
#   geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed",
#                               color = "red") +
#   labs(x = "Plasma Neonic Detection",
#        y = "Fattening Index") + 
#   my_theme

# view Fattening Index ~ Species
# there is variation within and across species
# birds %>%
#   ggplot(aes(x = Species, y = FatteningIndex)) +
#   geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed",
#                               color = "red") +
#   labs(x = "Species",
#        y = "Fattening Index") + 
#   my_theme


# interestingly, not significant by itself
# m4 <- lm(FatteningIndex ~ PlasmaDetection, data = birds)
# summary(m4)

# ...fat model ----
# tried species as a random effect, 
m5 <- glmmTMB(Fat ~ time_hours + MigStatus + PercentAg + Season + SPEI +
                PlasmaDetection + (1|Species), data = birds.c)


m6 <- glmmTMB(PecSizeBest ~ time_hours + MigStatus + PercentAg + Season + SPEI +
                PlasmaDetection + (1|Species), data = birds.c)


#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 


# run piecewiseSEM (fattening index and BCI only)
model <- psem(m1,m2,m3,m4)

model <- psem(m1,m2,m3,m4)
summary(model, conserve = TRUE)

# run piecewiseSEM (add in fat)
model <- psem(m1,m2,m3,m4,m5,
              Fat %~~% BCI)
summary(model, conserve = TRUE) # different number of rows...


# run piecewiseSEM (add in pectoral muscle size)
model <- psem(m1,m2,m3,m4,m6,
              PecSizeBest %~~% BCI,
              BCI %~~% MigStatus,
              PecSizeBest %~~% FatteningIndex,
              EnvDetection %~~% time_hours)
summary(model, conserve = TRUE) # different number of rows...



# ...model comparison
model_names <- paste0("m", 1:7)
models <- mget(model_names)
aictab(models, modnames = model_names)



# Generate example data
dat <- data.frame(x1 = runif(50),
                  x2 = runif(50), y1 = runif(50),
                  y2 = runif(50))
# Create list of structural equations
sem <- psem(
  lm(y1 ~ x1 + x2, dat),
  lm(y2 ~ y1 + x1, dat)
)
summary(sem)
# Look at correlated error between x1 and x2
# (exogenous)
cerror(x1 %~~% x2, sem, dat)
# Same as cor.test
with(dat, cor.test(x1, x2))
# Look at correlated error between x1 and y1
# (endogenous)
cerror(y1 %~~% x1, sem, dat)
# Not the same as cor.test
# (accounts for influence of x1 and x2 on y1)
with(dat, cor.test(y1, x1))
# Specify in psem
sem <- update(sem, x1 %~~% y1)
coefs(sem)



#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# I think issues here are related to small sample size...
# wetland pesticide detection model
simulationOutput <- simulateResiduals(fittedModel = m1) # overall good
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # pattern

plotResiduals(simulationOutput, form = model.frame(m1)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m1)$Season) # good
plotResiduals(simulationOutput, form = model.frame(m1)$SPEI) # some pattern

# plasma detection model --- all great, no patterns
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
plotResiduals(simulationOutput, form = model.frame(m2)$SPEI) # good

# body condition index model --- all good, no patterns
simulationOutput <- simulateResiduals(fittedModel = m3) 
plot(simulationOutput)
testDispersion(m3) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m3)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m3)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m3)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m3)$SPEI) # good


# fattening index model --- overall good, some pattern
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$BCI) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # some pattern



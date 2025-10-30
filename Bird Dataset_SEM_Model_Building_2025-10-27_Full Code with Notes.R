#-----------------------------------------#
#           SEM Model Building            #
#   Drivers of Shorebird Body Condition   #
# Created by Shelby McCahon on 10/27/2025 #
#         Modified on 10/30/2025          #
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

# Only include species with at least three individuals
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# add a 1 to each value of Fat to avoid zeros (gamma can now be used)
# birds <- birds %>% 
#   mutate(Fat = Fat + 0.001)

# extract one row per site to avoid pseudoreplication in analysis (n = 24 wetlands)
site_data <- birds %>%
  distinct(Site, SPEI, PercentAg, Season, EnvDetection)

birds <- birds %>%
   filter(complete.cases(BCI.NoEvent, Standardized.Pec.NoEvent, 
                         PlasmaDetection))


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

# removed season to avoid overfitting the data
# model without season is also preferred (AICcwt = 0.84), even though
# that claim is significant in dSEP test
m1 <- glm(EnvDetection ~ PercentAg + SPEI, data = site_data, 
          family = "binomial")

# m2 <- glm(EnvDetection ~ PercentAg + SPEI + Season, data = site_data, 
#            family = "binomial")
#  
#  model_names <- paste0("m", 1:2)
#  models <- mget(model_names)
#  aictab(models, modnames = model_names)


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
# is random effect needed? no (AICc wt in model without random effect = 0.59)
# but there is enough variance where it could be necessary

# WITHOUT RANDOM EFFECT
# m2 <- glm(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
#                 MigStatus,
#               data = birds, family = "binomial")

# 
# # WITH RANDOM EFFECT
m2 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
             MigStatus + (1|Species),
           data = birds, family = "binomial")
# 
# # WITHOUT RANDOM EFFECT
# m2 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
#                  MigStatus,
#                data = birds, family = "binomial")
 
 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

# random effect of species needed? no, but results differ...
# m1 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
#                 MigStatus + (1|Species),
#               data = birds, family = "binomial")
# 
# m2 <- glmmTMB(PlasmaDetection ~ PercentAg + SPEI + EnvDetection + Season +
#                 MigStatus,
#               data = birds, family = "binomial")
# 
# model_names <- paste0("m", 1:2)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

# ...body condition index model ----
# site as a random effect caused model to fail to converge
# model explains almost nothing so you get an error (but model is still valid)
# no random effect of species because species is already considered in 
# residual analysis (log(Mass)~log(Wing) + Species)

# adding in pec index gave nonsensical results...decided to omit here
# d sep did not flag it as necessary to include anyways
# adding in fattening index also made model significantly worse...
m3 <- lm(BCI.NoEvent ~ time_hours + PlasmaDetection + SPEI + PercentAg + 
           Season,
         data = birds)


# ...fattening index model ----
# species performs better than migratory status alone (AICc wt = 0.96)
# species also performs better than mig. status + (1|Species) (wt = 0.80)
# issue is, I can't include species as a factor in SEM
# best I can do is species as a random effect and migratory status as fixed
# effect of plasma detection on fattening index is still significant with
# species as a random or fixed effect!
# including body condition in fattening index model significantly improved fit

# WITH RANDOM EFFECT
m4 <- glmmTMB(FatteningIndex ~ MigStatus + Season + SPEI + 
                PercentAg + time_hours + PlasmaDetection + EnvDetection + 
                # BCI.NoEvent +
                (1|Species), data = birds)

# m4 <- lm(FatteningIndex ~ BCI.NoEvent + Season + SPEI + 
#                 PercentAg + time_hours + PlasmaDetection + EnvDetection + 
#                 Species, data = birds)

# WITHOUT RANDOM EFFECT
# m4 <- lm(FatteningIndex ~ BCI.NoEvent + MigStatus + Season + SPEI + PercentAg +
#                  time_hours + PlasmaDetection + EnvDetection, data = birds)


# view Fattening Index ~ Plasma Detection
# birds %>%
#   filter(!is.na(PlasmaDetection) & !is.na(FatteningIndex)) %>%
#   ggplot(aes(x = as.factor(PlasmaDetection), y = FatteningIndex,
#              color = as.factor(MigStatus))) +
#    geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed",
#                                color = "red") +
#    labs(x = "Plasma Neonic Detection",
#         y = "Fattening Index") + 
#   geom_point() +
#    my_theme
# 
# # view Fattening Index ~ PercentAg
# birds %>%
#   filter(!is.na(PlasmaDetection) & !is.na(PercentAg)) %>%
#   ggplot(aes(x = PercentAg, y = FatteningIndex)) +
#   geom_hline(yintercept = 0, linetype = "dashed",
#                               color = "red") +
#   labs(x = "% Surrounding Cropland",
#        y = "Fattening Index") + 
#   geom_point() +
#   my_theme
# 
# 
# birds %>% select(FatteningIndex,Species)
# 
# cor(birds$Mass, birds$FatteningIndex, use = "complete") # low-moderate cor.

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
# can't get a good model fit...what to do with 0s?
# results change quite a bit with the distribution
# heavily right skewed
# m5 <- glm(Fat ~ PercentAg + Season + SPEI +
#                 PlasmaDetection + MigStatus, data = birds,
#               family = Gamma(link = "log"))

# # Deviance residuals vs fitted values
# plot(m5$fitted.values, residuals(m5, type = "deviance"),
#      xlab = "Fitted values", ylab = "Deviance residuals")
# abline(h = 0, col = "red")
# 
# 

# # view relationships
# ggplot(birds, aes(x = SPEI, y = Fat, color = MigStatus)) + 
#   geom_point() + my_theme
# 
# table(birds$Fat)
# table(birds$MigStatus)
# table(birds$Fat, birds$MigStatus) # a lot more migrants (more than double)
# table(birds$Fat, birds$Season) # very low fat in spring...
# table(birds$Fat, birds$PlasmaDetection)
# 
# # plot model
# d <- expand.grid(PercentAg = mean(birds$PercentAg),
#                  SPEI = seq(min(birds$SPEI),
#                                  max(birds$SPEI),
#                                  length = 1000),
#                  MigStatus = unique(birds$MigStatus)) 
# 
# pred <- predict(m5, newdata = d, type = "link", se.fit = TRUE)
# 
# d$fit <- exp(pred$fit)
# d$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
# d$upr <- exp(pred$fit + 1.96 * pred$se.fit)
# 
# ggplot(d, aes(x = SPEI, y = fit)) +
#   geom_line(linewidth = 0.8) +
#   geom_ribbon(aes(ymin = lwr, ymax = upr),
#               alpha = 0.25, color = NA, show.legend = FALSE) +
#   theme_classic() +
#   labs(y = "Predicted Fat", x = "SPEI")

# ...pectoral muscle size model ----
# corrected for size and species -- random effect and mig. status not needed

# fattening index not significant here
m6 <- lm(Standardized.Pec.NoEvent ~ time_hours + PercentAg + SPEI + 
                Season + PlasmaDetection, 
              data = birds)

# model_names <- paste0("m", 6:7)
# models <- mget(model_names)
# aictab(models, modnames = model_names)


# ...uric acid level model ----
# not enough variation for random effect, species was removed
# migratory status is an important predictor --> keep
# AICc favors model with migratory status and NO random effect of species
# results don't differ much with or without random effect of species
# keep random effect in for consistency with model (species is now always
# accounted for)

# WITH RANDOM EFFECT OF SPECIES
m7 <- glmmTMB(Uric ~ time_hours + PercentAg + Season + SPEI +
                PlasmaDetection + MigStatus + (1|Species), data = birds)
# 
# m7 <- glmmTMB(Uric ~ time_hours + PercentAg + Season + SPEI +
#              PlasmaDetection + MigStatus, data = birds)
# #  
#   m8 <- glmmTMB(Uric ~ time_hours + MigStatus + PercentAg + Season + SPEI +
#              PlasmaDetection + (1|Species), data = birds)
# #  
#   m9 <- glmmTMB(Uric ~ time_hours + PercentAg + Season + SPEI +
#                   PlasmaDetection + (1|Species), data = birds)
# #  
#   model_names <- paste0("m", 7:9)
#   models <- mget(model_names)
#   aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#                           run piecewise SEM                               ----                        
#------------------------------------------------------------------------------# 

# run piecewiseSEM
# not including fat due to model diagnostics issues
model <- psem(m1,m2,m3,m4,m6,m7)
summary(model, conserve = TRUE)
# sem <- update(model, EnvDetection %~~% Uric)

# summary(sem)

# ...model comparison
# model_names <- paste0("m", 6:7)
# models <- mget(model_names)
# aictab(models, modnames = model_names)

#------------------------------------------------------------------------------#
#                           model diagnostics                               ----                        
#------------------------------------------------------------------------------# 

# wetland pesticide detection model 
# (DHARMa diagnostics unreliable with effect size of 11)
plot(residuals(m1, type = "deviance") ~ fitted(m1))
abline(h = 0, col = "red") # reasonable fit

# plasma detection model --- all great, no patterns
simulationOutput <- simulateResiduals(fittedModel = m2) 
plot(simulationOutput)
testDispersion(m2) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m2)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m2)$MigStatus)  # good
plotResiduals(simulationOutput, form = model.frame(m2)$EnvDetection)  # some pattern
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
plotResiduals(simulationOutput, form = model.frame(m3)$FatteningIndex) # good

# fattening index model --- overall good, some pattern
simulationOutput <- simulateResiduals(fittedModel = m4) 
plot(simulationOutput)
testDispersion(m4) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) # some pattern

plotResiduals(simulationOutput, form = model.frame(m4)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m4)$MigStatus)  # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$time_hours)  # cyclical but okay
plotResiduals(simulationOutput, form = model.frame(m4)$Season) # some pattern
plotResiduals(simulationOutput, form = model.frame(m4)$BCI.NoEvent) # good
plotResiduals(simulationOutput, form = model.frame(m4)$SPEI) # some pattern

plot(residuals(m4, type = "deviance") ~ fitted(m4))
abline(h = 0, col = "red") # residuals look just fine though

# fat model --- lots of issues: can't get a good model fit right now
simulationOutput <- simulateResiduals(fittedModel = m5) 
plot(simulationOutput)
testDispersion(m5) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m5)$PercentAg) # some pattern
plotResiduals(simulationOutput, form = model.frame(m5)$PlasmaDetection) # bad
plotResiduals(simulationOutput, form = model.frame(m5)$MigStatus)  # bad
plotResiduals(simulationOutput, form = model.frame(m5)$time_hours)  # good
plotResiduals(simulationOutput, form = model.frame(m5)$Season) # very bad
plotResiduals(simulationOutput, form = model.frame(m5)$BCI) # good
plotResiduals(simulationOutput, form = model.frame(m5)$SPEI) # some pattern

# pectoral muscle model --- all good, no patterns
simulationOutput <- simulateResiduals(fittedModel = m6) 
plot(simulationOutput)
testDispersion(m6) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m6)$PercentAg)
plotResiduals(simulationOutput, form = model.frame(m6)$PlasmaDetection)
plotResiduals(simulationOutput, form = model.frame(m6)$time_hours)  
plotResiduals(simulationOutput, form = model.frame(m6)$SPEI)

# uric acid model --- all good, no patterns
simulationOutput <- simulateResiduals(fittedModel = m7) 
plot(simulationOutput)
testDispersion(m7) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = model.frame(m7)$PercentAg) # good
plotResiduals(simulationOutput, form = model.frame(m7)$PlasmaDetection) # good
plotResiduals(simulationOutput, form = model.frame(m7)$time_hours) # good
plotResiduals(simulationOutput, form = model.frame(m7)$SPEI) # good
plotResiduals(simulationOutput, form = model.frame(m7)$MigStatus) # good
plotResiduals(simulationOutput, form = model.frame(m7)$Season) # good


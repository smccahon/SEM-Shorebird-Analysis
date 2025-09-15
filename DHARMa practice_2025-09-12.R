#----------------------------------------#
#             DHARMa Practice            #
# Created by Shelby McCahon on 9/12/2025 #
#         Modified on 9/15/2025          #
#----------------------------------------#

# load packages ----

library(tidyverse)
library(lme4)
library(DHARMa)
library(car)
library(glmmTMB)

# load data ----

birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")

# convert factors to numeric ----

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


# build models ----
birds <- na.omit(birds[, c("PlasmaDetection", "EnvDetection", 
                                 "PercentAg")])

fittedModel <- glm(PlasmaDetection ~ PercentAg + EnvDetection, 
                   family = "binomial", data = birds)

#---

# model diagnostics (DHARMa) ----

# test for dispersion
testDispersion(fittedModel)

# calculate scaled residuals
simulationOutput <- simulateResiduals(fittedModel = fittedModel) 

# access scaled residuals
residuals(simulationOutput) 

# plot scaled residuals
plot(simulationOutput)

# left plot only -- detects overall deviations from the expected distribution
# KS test = test for correct distribution
plotQQunif(simulationOutput)

# right plot only -- residuals against the predicted value/variable
# red stars: simulation outliers (points outside of range of simulated values)
plotResiduals(simulationOutput)

# plot residuals against other predictors (highly recommended to do)
# need to remove all NAs prior to model fitting (same results though)
# ideal to see flat black line (red is violation)
plotResiduals(simulationOutput, form = birds$PercentAg)

#---

# goodness of fit tests on the scaled residuals ----

# tests if the overall distribution conforms to expectations
testUniformity(simulationOutput) 

# tests if there are more simulation outliers than expected
testOutliers(simulationOutput)

# tests if the simulated dispersion is equal to the observed dispersion
testDispersion(simulationOutput)

# fits residuals against a predictor
testQuantiles(simulationOutput)

# tests residuals against a categorical predictor
testCategorical(simulationOutput, catPred = birds$EnvDetection)

# tests if there are more zeros in the data than expected from simulations
testZeroInflation(simulationOutput)

#---

# zero-inflation (glmmTMB package) ----

fittedModelZERO <- glmmTMB(PlasmaDetection ~ PercentAg + EnvDetection + (1|Site), 
                   family = "binomial", data = birds)

summary(fittedModelZERO)

simulationOutputZERO <- simulateResiduals(fittedModel = fittedModelZERO)

plot(simulationOutputZERO)

#---

# heteroscedasticity ----

testQuantiles(simulationOutput) # exact p-values for the quantile lines

testQuantiles(simulationOutput, predictor = birds$PercentAg)


# test for categorical predictors
testCategorical(simulationOutput, catPred = birds$EnvDetection)

testCategorical(simulationOutput, catPred = birds$PlasmaDetection)


#---

# detecting missing predictors or wrong functional assumptions ----

# retrieve residuals
simulationOutput$scaledResiduals

plotResiduals(simulationOutput, birds$PercentAg)


# more help ----
?simulateResiduals





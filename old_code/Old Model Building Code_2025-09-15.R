#----------------------------------------#
#        Original Code Documentation     #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/17/2025          #
#----------------------------------------#

#------------------------------------------------------------------------------#
#                        code from August 2025                              ----                        
#------------------------------------------------------------------------------# 

#---

full.1 <- full %>% 
  mutate(across(c(PercentAg, DaysSinceLastPrecipitation_5mm), 
                ~ as.numeric(scale(.))))

# site as a random effect
m1 <- glmer(EnvDetection ~ PercentAg +
              DaysSinceLastPrecipitation_5mm + (1 | Site), 
            family = "binomial", data = full.1)

# site not accounted for
m1 <- glm(EnvDetection ~ PercentAg +
            DaysSinceLastPrecipitation_5mm, 
          family = "binomial", data = full.1)


#---

full.2 <- full %>% 
  mutate(across(c(seconds_since_midnight, PercentAg, Julian), 
                ~ as.numeric(scale(.))))

# site as a random effect
m2 <- glmer(PlasmaDetection ~ EnvDetection + seconds_since_midnight + 
              PercentAg + Julian + (1 | Site), family = "binomial", 
            data = full.2)

# site not accounted for
m2 <- glm(PlasmaDetection ~ EnvDetection + seconds_since_midnight + 
            PercentAg + Julian, family = "binomial", 
          data = full.2)

#---

full.3 <- full %>% 
  mutate(across(c(PercentLocalVeg_50m, PercentAg), 
                ~ as.numeric(scale(.))))

# site as a random effect
m3 <- lmer(Biomass ~ PercentAg + EnvDetection +
             PercentLocalVeg_50m + WaterQuality + (1 | Site),
           data = full.3)

# site not accounted for
m3 <- lm(Biomass ~ PercentAg + EnvDetection +
           PercentLocalVeg_50m + WaterQuality,
         data = full.3)

#---

birds.4 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, Julian, SPEI), 
                ~ as.numeric(scale(.))))

# site as a random effect
m4 <- lmer(BodyCondition ~ seconds_since_midnight + Biomass + Julian + SPEI +
             PlasmaDetection + (1 | Site), 
           data = birds.4)

# site not accounted for
m4 <- lm(BodyCondition ~ seconds_since_midnight + Biomass + Julian + SPEI +
           PlasmaDetection, 
         data = birds.4)


#---

invert.5 <- invert %>% 
  mutate(across(c(PercentAg, Julian, AnnualSnowfall_in), 
                ~ as.numeric(scale(.))))

# site as a random effect
m5 <- lm(PercentLocalVeg_50m ~ PercentAg + Julian + Permanence +
           AnnualSnowfall_in, 
         data = invert.5)

# site not accounted for
m5 <- lm(PercentLocalVeg_50m ~ PercentAg + Julian + Permanence +
           AnnualSnowfall_in, 
         data = invert.5)

#---

birds.6 <- birds %>% 
  mutate(across(c(seconds_since_midnight, Biomass, SPEI, Julian),
                ~ as.numeric(scale(.))))

# site as a random effect
m6 <- lmer(FatteningIndex ~ seconds_since_midnight + Biomass + PlasmaDetection +
             SPEI + Julian + (1 | Site), data = birds.6)

# site not accounted for
m6 <- lm(FatteningIndex ~ seconds_since_midnight + Biomass + PlasmaDetection +
           SPEI + Julian, data = birds.6)

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

model <- psem(m1, m2, m3, m4, m5, m6)
summary(model)

# model diagnostics

# ...environmental detection model (ISSUES)
DHARM <- simulateResiduals(fittedModel = m1)
plot(DHARM)
plotResiduals(DHARM, full$DaysSinceLastPrecipitation_5mm)
testDispersion(DHARM) # no overdispersion
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation

# ...plasma detection model (GOOD) 
DHARM <- simulateResiduals(fittedModel = m3)
plot(DHARM)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation


# ...biomass model (GOOD) 
plot(residuals(m3), main = "Deviance Residuals", ylab = "Residuals", 
     xlab = "Index") # no clear pattern

#...body condition model 
DHARM <- simulateResiduals(fittedModel = m4)
plot(DHARM) # combined adjusted quantile test significant
plotResiduals(DHARM, form = birds.4$Julian)
testDispersion(DHARM) # overdispersion: 1.00
testOutliers(DHARM) # 0 outliers
testUniformity(DHARM) # no evidence of non-uniformity of residuals
testZeroInflation(DHARM) # no evidence of zero-inflation







#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
#                 initial biomass models for sem (9/17/2025)                ----
#------------------------------------------------------------------------------#

# transforming 0's to small numbers to allow the use of a Gamma distribution
invert$Biomass_adj <- ifelse(invert$Biomass == 0, 0.0001, 
                              invert$Biomass)
 
birds$Biomass_adj <- ifelse(birds$Biomass == 0, 0.0001, 
                             birds$Biomass)
 
m1 <- glm(Biomass_adj ~ PercentAg,
               family = Gamma(link = "log"), data = invert)
 
m2 <- lm(BodyCondition ~ Biomass_adj, data = birds)

# psem expects a random effect in glmmTMB so I created a dummy random effect
invert$dummy <- 1

# tweedie model (psem does not support)
invert$log_PercentAg <- log(invert$PercentAg)
m1 <- glmmTMB(Biomass ~ log_PercentAg + Julian + (1 | dummy), 
               data = invert, 
               family = glmmTMB::tweedie(link = "log"))

# zero-inflated Gamma model (psem does not support)
m1 <- glmmTMB(Biomass ~ log_PercentAg + Julian + (1 | dummy),
                 ziformula = ~ log_PercentAg + Julian + (1 | dummy),
                 family = ziGamma(link = "log"),
                 data = invert)





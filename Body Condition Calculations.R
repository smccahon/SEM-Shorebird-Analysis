# load packages
library(tidyverse)


# load data
birds <- read.csv("original_data/shorebird_body_condition_data_2025-05-29.csv")

# subset by species
leye <- subset(birds, Species == "Lesser Yellowlegs") 
lesa <- subset(birds, Species == "Least Sandpiper") 
sesa <- subset(birds, Species == "Semipalmated Sandpiper") 
will <- subset(birds, Species == "Willet")
kill <- subset(birds, Species == "Killdeer")
pesa <- subset(birds, Species == "Pectoral Sandpiper")
amav <- subset(birds, Species == "American Avocet")
wiph <- subset(birds, Species == "Wilsons Phalarope")
lbdo <- subset(birds, Species == "Longbilled Dowitcher")
grye <- subset(birds, Species == "Greater Yellowlegs")

# regression
m <- lm(Mass ~ Fat, data = grye)
summary(m)

plot(lbdo$Mass, lbdo$Fat)


# correlation summary by species (mass ~ fat)
# LEYE: (+)
# LESA: (+)
# SESA: (+)
# WILL: (+)
# KILL: no correlation (only one bird with a fat score > 0)
# PESA: no correlation
# AMAV: no correlation (all had fat of 0)
# WIPH: no correlation (most had fat of 0)
# LBDO: no correlation (most had fat of 0)

# ---

# pectoral muscle and mass
m <- lm(Mass ~ PecSizeBest, data = grye)
summary(m)

plot(grye$Mass, grye$PecSizeBest) 

# correlation summary by species (mass ~ pectoral muscle size)
# LEYE: (+)
# LESA: (marginally positive)
# SESA: no correlation (but looks positive)
# WILL: no correlation (but looks positive)
# KILL: no correlation (but looks positive)
# PESA: (marginally positive)
# AMAV: no correlation
# WIPH: (+)
# LBDO: no correlation but looks positive



# PCA to figure out which structural trait to use (all about the same, so use wing?)
# R example
structural_traits <- na.omit(birds[, c("Wing", "DiagTarsus", "Culmen", "Head")])
pca_struct <- prcomp(structural_traits, scale. = TRUE)

summary(pca_struct)         # Variance explained
pca_struct$rotation         # Loadings of each trait on PC1



# Option 1. Regression of mass and wing chord

# 1. Fit model
birds.complete <- na.omit(birds[, c("Mass", "Species", "Wing", "Julian",
                                    "Event")])

birds.complete$Species <- as.factor(birds.complete$Species)

m1 <- lmer(Mass ~ as.factor(Species) + Wing + Julian + (1|Event), data = birds.complete,
           na.action = na.omit)

m1 <- lm(Mass ~ Wing, data = birds)

# 2. Create a new residual column with NAs
birds$residual <- NA

# 3. Get the rows that were used in the model
used_rows <- as.numeric(rownames(model.frame(m)))

# 4. Assign residuals only to those rows
birds$residual[used_rows] <- residuals(m)

simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testUniformity(simulationOutput)
testOutliers(simulationOutput) 
testQuantiles(simulationOutput) 

plotResiduals(simulationOutput, form = birds.complete$Wing)


# Option 2: Scaled Mass Index (lumping all birds together)

# 1. Calculate log-transformed variables
birds$logMass <- log(birds$Mass)
birds$logWing <- log(birds$Wing)

# 2. Fit regression to get scaling exponent 'b'
model <- lm(logMass ~ logWing, data = birds)
b <- coef(model)["logWing"]

# 3. Calculate mean reference size L0
L0 <- mean(birds$Wing, na.rm = TRUE)

# 4. Calculate Scaled Mass Index
birds$SMI <- birds$Mass * (L0 / birds$Wing)^b

plot(birds$Wing, birds$Mass, main = "Mass vs Wing Length")
points(birds$Wing, birds$SMI, col = "red", pch = 19)
legend("topleft", legend = c("Original Mass", "SMI"), col = c("black", "red"), pch = c(1,19))

ggplot(birds, aes(x = Species, y = SMI)) + geom_boxplot()

m1 <- lm(SMI ~ Species, data = birds)


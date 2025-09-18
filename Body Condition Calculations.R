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


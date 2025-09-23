#----------------------------------------#
#          SEM Data Preparation          #
# Created by Shelby McCahon on 8/01/2025 #
#         Modified on 9/18/2025          #
#----------------------------------------#

# load packages
library(tidyverse)
library(purrr)
library(lme4)
library(dplyr)


# load data
birds <- read.csv("original_data/shorebird_body_condition_data_2025-05-29.csv")
invert <- read.csv("original_data/macroinvertebrate_data_2025-08-22.csv")
full <- read.csv("original_data/full_model_dataset_2025-08-29.csv")

#------------------------------------------------------------------------------#
#                         Shorebird Data Preparation                        ----                        
#------------------------------------------------------------------------------#   

### ...standardize time to something more simple -------------------------------
birds$Time <- strptime(birds$Time, format = "%H:%M")
birds$Time <- as.POSIXct(birds$Time, tz = "America/Chicago")
attributes(birds$Time)$tzone

# Calculate seconds since midnight
birds$seconds_since_midnight <- as.numeric(difftime(birds$Time, 
                                                    floor_date(birds$Time, "day"), 
                                                    units = "secs"))

birds$sin <- sin(2 * pi * birds$seconds_since_midnight / (24 * 3600))
birds$cos <- cos(2 * pi * birds$seconds_since_midnight / (24 * 3600))

birds$formatted_time <- format(as.POSIXct(birds$seconds_since_midnight, 
                                          origin = "1970-01-01", tz = "UTC"), 
                               "%H:%M")

# fix numeric instability issue in SEM modeling
birds$time_hours <- birds$seconds_since_midnight / 3600 


### ...create body condition index (accounting for structural size ) -----------

# calculate species-specific SMI
birds <- birds %>%
  group_by(Species) %>%
  mutate(
    logMass = log(Mass),
    logWing = log(Wing)
  ) %>%
  group_modify(~ {
    mod <- lm(logMass ~ logWing, data = .x)
    b <- coef(mod)["logWing"]
    L0 <- mean(.x$Wing, na.rm = TRUE)
    .x %>%
      mutate(SMI = Mass * (L0 / Wing)^b)
  }) %>%
  ungroup()

# standardize SMI within each species
birds <- birds %>%
  group_by(Species) %>%
  mutate(
    SMI = (SMI - mean(SMI, na.rm = TRUE)) / sd(SMI, na.rm = TRUE)
  ) %>%
  ungroup()



# view results
excluded_species <- c("Marbled Godwit", "Shortbilled Dowitcher")

birds_subset <- birds %>%
  filter(!(Species %in% excluded_species))
  
ggplot(birds_subset, aes(x = Species, y = SMI)) + geom_boxplot() +
  theme_classic() + geom_hline(yintercept = 0, col = "red", size = 1,
                               linetype = "dashed") +
  labs(x = NULL, y = "Scaled Mass Index (SMI)") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x))


# standardize fat and pectoral score
birds <- birds %>%
  group_by(Species) %>%
  mutate(
    Fat_z = scale(Fat),
    Pec_z = scale(PecSizeBest)
  ) %>%
  ungroup()

condition_data <- birds %>%
  select(SMI, Fat_z, Pec_z) %>%
  drop_na()

pca <- prcomp(condition_data, scale. = TRUE)

# Add PC1 back to original birds dataset
birds$Condition_PC1 <- NA
birds$Condition_PC1[as.numeric(rownames(condition_data))] <- pca$x[, 1]


m1 <- lm(SMI ~ 1, data = birds)
m2 <- lm(Condition_PC1 ~ Biomass, data = birds)
summary(m2)

simulationOutput <- simulateResiduals(fittedModel = m1) 
plot(simulationOutput)
testDispersion(m1) 
testZeroInflation(m1)
testUniformity(simulationOutput) 
testOutliers(simulationOutput) 
testQuantiles(simulationOutput)

### ...create fattening index with PCA -----------------------------------------
# high tri and low beta for high fattening

# 1. Subset the dataset to only include relevant variables
birds_subset <- birds[, c("Beta", "Tri")]

# 2. Remove rows with NAs
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ]  # n = 89

# 3. Check correlation (optional)
correlation <- cor(birds_subset_clean)
print(correlation)  # -0.23

# 4. Run PCA on the cleaned data (centered and scaled)
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# 5. Check PCA summary
print(summary(pca_result))

# 6. Extract PCA scores (principal components)
pca_scores <- pca_result$x

# 7. Flip the sign of PC1 so higher scores indicate fattening
pca_scores[, 1] <- -pca_scores[, 1]

# 8. Initialize FatteningIndex column with NA in original dataset
birds$FatteningIndex <- NA

# 9. Add the flipped PC1 scores for rows without missing Beta or Tri
birds[complete.cases(birds_subset), "FatteningIndex"] <- pca_scores[, 1]

# 10. View PCA rotation/loadings
print(pca_result$rotation)

tri_cor <- cor(pca_scores[, 1], birds_subset_clean$Tri)
beta_cor <- cor(pca_scores[, 1], birds_subset_clean$Beta)

cat("Correlation with Tri:", round(tri_cor, 3), "\n")   # Should be positive
cat("Correlation with Beta:", round(beta_cor, 3), "\n") # Should be negative

m <- lm(FatteningIndex ~ Species, data = birds)


### ...trim down data and export file for analysis -----------------------------
birds_cleaned <- birds %>% 
  select(Individual, Site, Event, Season, Age, Fat, PecSizeBest, Beta, Permanence, 
         NearestCropDistance_m, Max_Flock_Size,
         PlasmaDetection, seconds_since_midnight, Species, PercentAg, SMI,
         Percent_Total_Veg, Julian, Mass, OverallNeonic,
         Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
         AnyDetection, MigStatus, AgCategory, SPEI, Sex, Tri, 
         Diversity, Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
         InvertPesticideDetection, EnvDetection, FatteningIndex, time_hours,
         AnnualSnowfall_in, PrecipitationAmount_7days, DaysSinceLastPrecipitation_5mm)

write.csv(birds_cleaned, "cleaned_data/shorebird_data_cleaned_2025-08-11.csv",
          row.names = FALSE)

#------------------------------------------------------------------------------#
#                         Invert Data Preparation                           ----                        
#------------------------------------------------------------------------------#   

### ...create water quality index with PCA -------------------------------------

# subset the dataset to only include relevant variables
invert_subset <- invert[, c("Conductivity_uS.cm", "Salinity_ppt", 
                               "TDS_mg.L")]

# remove rows with NAs
invert_subset_clean <- invert_subset[complete.cases(invert_subset), ] # n = 152
cor(invert_subset_clean) # ranges from 0.98-0.99

# Run PCA on the cleaned data
pca_result <- prcomp(invert_subset_clean, center = TRUE, scale. = TRUE)

# 99% explained by PC1
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
invert$WaterQuality <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
invert[complete.cases(invert_subset), "WaterQuality"] <- pca_scores[, 1]

### ...trim down data and export file for analysis -----------------------------
invert_cleaned <- invert %>% 
  select(Wetland, Season, pH_probe, Biomass, PesticideInvert_ng.g, 
         WaterNeonicDetection, DominantCrop, NearestCropDistance_m,
         Dist_Closest_Wetland_m, NearestCropType, 
         PercentLocalVeg_50m, Julian, PercentAg, Buffered,
         WaterQuality, Diversity, InvertPesticideDetection, Permanence,
         NearestCropDistance_m, AgCategory, ShorebirdsSeen,
         NeonicInvert_ng.g, NeonicWater_ng.L, DominantCrop, NearestCropType,
         PrecipitationAmount_7days, DaysSinceLastPrecipitation_5mm, 
         AnnualSnowfall_in, EnvDetection)

write.csv(invert_cleaned, "cleaned_data/invert_data_cleaned_2025-08-11.csv",
          row.names = FALSE)






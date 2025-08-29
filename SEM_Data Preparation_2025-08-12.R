#----------------------------------------#
#          SEM Data Preparation          #
# Created by Shelby McCahon on 8/01/2025 #
#         Modified on 8/29/2025          #
#----------------------------------------#

# load packages
library(tidyverse)

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

### ...create body condition index with PCA ------------------------------------

# subset the dataset to only include relevant variables
birds_subset <- birds[, c("Fat", "Mass", "PecSizeBest")]

# remove rows with NAs
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ] # n = 152
cor(birds_subset_clean) # ranges from -0.13 to 0.36

# Run PCA on the cleaned data
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# 52% explained by PC1, 29% explained by PC2, 19% explained by PC3
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds$BodyCondition <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds[complete.cases(birds_subset), "BodyCondition"] <- pca_scores[, 1]

### ...create fattening index with PCA -----------------------------------------
# subset the dataset to only include relevant variables
birds_subset <- birds[, c("Beta", "Tri")]

# remove rows with NAs
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ] # n = 89
cor(birds_subset_clean) # -0.23

# Run PCA on the cleaned data
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# 61% explained by PC1, 39% explained by PC2
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
birds$FatteningIndex <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
birds[complete.cases(birds_subset), "FatteningIndex"] <- pca_scores[, 1]

### ...trim down data and export file for analysis -----------------------------
birds_cleaned <- birds %>% 
  select(Individual, Site, Event, Season, Age, Fat, PecSizeBest, Beta, Permanence, 
         NearestCropDistance_m, Max_Flock_Size,
         PlasmaDetection, seconds_since_midnight, Species, PercentAg, 
         Percent_Total_Veg, Julian, Mass, OverallNeonic,
         Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
         AnyDetection, BodyCondition, MigStatus, AgCategory, SPEI, Sex, Tri, 
         Diversity, Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
         InvertPesticideDetection, EnvDetection, FatteningIndex,
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


#------------------------------------------------------------------------------#
#                         Full Dataset Preparation                          ----                        
#------------------------------------------------------------------------------#   

### ...standardize time to something more simple -------------------------------
full$Time <- strptime(full$Time, format = "%H:%M")
full$Time <- as.POSIXct(full$Time, tz = "America/Chicago")
attributes(full$Time)$tzone

# Calculate seconds since midnight
full$seconds_since_midnight <- as.numeric(difftime(full$Time, 
                                                    floor_date(full$Time, "day"), 
                                                    units = "secs"))

full$sin <- sin(2 * pi * full$seconds_since_midnight / (24 * 3600))
full$cos <- cos(2 * pi * full$seconds_since_midnight / (24 * 3600))

full$formatted_time <- format(as.POSIXct(full$seconds_since_midnight, 
                                          origin = "1970-01-01", tz = "UTC"), 
                               "%H:%M")

### ...create body condition index with PCA ------------------------------------

# subset the dataset to only include relevant variables
full_subset <- full[, c("Fat", "Mass", "PecSizeBest")]

# remove rows with NAs
full_subset_clean <- full_subset[complete.cases(full_subset), ] # n = 123
cor(full_subset_clean) # ranges from -0.45 to 0.38

# Run PCA on the cleaned data
pca_result <- prcomp(full_subset_clean, center = TRUE, scale. = TRUE)

# 56% explained by PC1, 27% explained by PC2, 17% explained by PC3
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
full$BodyCondition <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
full[complete.cases(full_subset), "BodyCondition"] <- pca_scores[, 1]

### ...create fattening index with PCA -----------------------------------------

# subset the dataset to only include relevant variables
full_subset <- full[, c("Beta", "Tri")]

# remove rows with NAs
full_subset_clean <- full_subset[complete.cases(full_subset), ] # n = 89
cor(full_subset_clean) # -0.25

# Run PCA on the cleaned data
pca_result <- prcomp(full_subset_clean, center = TRUE, scale. = TRUE)

# 63% explained by PC1, 37% explained by PC2
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
full$FatteningIndex <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
full[complete.cases(full_subset), "FatteningIndex"] <- pca_scores[, 1]

### ...create water quality index with PCA -------------------------------------

# subset the dataset to only include relevant variables
# salinity and pH did not contribute anything meaningful and should prob.
# be treated separately
full_subset <- full[, c("Conductivity_uS.cm", "TDS_mg.L")]

# remove rows with NAs
full_subset_clean <- full_subset[complete.cases(full_subset), ] # n = 126
cor(full_subset_clean) # 0.84

# Run PCA on the cleaned data
pca_result <- prcomp(full_subset_clean, center = TRUE, scale. = TRUE)

# 92% explained by PC1
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
full$WaterQuality <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
full[complete.cases(full_subset), "WaterQuality"] <- pca_scores[, 1]

### ...trim down data and export file for analysis -----------------------------
full_cleaned <- full %>% 
  select(Individual, Site, Percent_Total_Veg, 
         PercentLocalVeg_50m, Season, Time,
         Age, Fat, PecSizeBest, Beta, Permanence, 
         NearestCropDistance_m, Max_Flock_Size,
         PlasmaDetection, seconds_since_midnight, Species, PercentAg, 
         Julian, Mass, OverallNeonic, WaterQuality, pH_probe, Salinity_ppt,
         Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
         AnyDetection, BodyCondition, MigStatus, AgCategory, SPEI, Sex, Tri, 
         Diversity, Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
         InvertPesticideDetection, EnvDetection, FatteningIndex,
         AnnualSnowfall_in, PrecipitationAmount_7days, 
         DaysSinceLastPrecipitation_5mm)

write.csv(full_cleaned, "cleaned_data/full_data_cleaned_2025-08-29.csv",
          row.names = FALSE)






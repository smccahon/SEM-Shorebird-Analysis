#----------------------------------------#
#          SEM Data Preparation          #
# Created by Shelby McCahon on 8/01/2025 #
#         Modified on 8/26/2025          #
#----------------------------------------#

# load packages
library(tidyverse)

# load data
birds <- read.csv("original_data/shorebird_body_condition_data_2025-05-29.csv")
invert <- read.csv("original_data/macroinvertebrate_data_2025-08-22.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")

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



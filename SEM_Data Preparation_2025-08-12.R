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

### ...create body condition index (accounting for structural size ) -----------

# 1. Define variable groups
size_vars <- c("Wing", "DiagTarsus", "Culmen")
condition_vars <- c("Mass", "Fat", "PecSizeBest")
all_vars <- c(size_vars, condition_vars)

# 2. Subset to complete cases only
birds_cc <- birds[complete.cases(birds[, all_vars]), ]

# 3. Create size index from PCA on structural size traits
size_pca <- prcomp(birds_cc[, size_vars], scale. = TRUE)
birds_cc$SizeIndex <- size_pca$x[, 1]  # PC1 = size axis

# 4. Regress each condition-related variable on size, extract residuals
birds_cc$resid_mass <- residuals(lm(Mass ~ SizeIndex + Species, data = birds_cc))
birds_cc$resid_fat <- residuals(lm(Fat ~ SizeIndex + Species, data = birds_cc))
birds_cc$resid_pec <- residuals(lm(PecSizeBest ~ SizeIndex + Species, data = birds_cc))

# 5. PCA on the residuals — gives size-corrected condition index
condition_pca <- prcomp(birds_cc[, c("resid_mass", "resid_fat", "resid_pec")], scale. = TRUE)

# 6. Extract PC1 as BodyCondition (higher = better condition)
birds_cc$BodyCondition <- condition_pca$x[, 1]

# 7. Add BodyCondition back to full birds dataset (set NA where missing)
birds$BodyCondition <- NA
birds[rownames(birds_cc), "BodyCondition"] <- birds_cc$BodyCondition

print(condition_pca$rotation)
summary(condition_pca)
plot(condition_pca)


# ### ...create body condition index (NOT accounting for structural size)
# # subset the dataset to only include relevant variables
# birds_subset <- birds[, c("Fat", "Mass", "PecSizeBest")]
# 
# # remove rows with NAs
# birds_subset_clean <- birds_subset[complete.cases(birds_subset), ] # n = 152
# cor(birds_subset_clean) # ranges from -0.13 to 0.36
# 
# # Run PCA on the cleaned data
# pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)
# 
# # 52% explained by PC1, 29% explained by PC2, 19% explained by PC3
# summary(pca_result)
# 
# # View the PCA scores (principal components)
# pca_scores <- pca_result$x
# 
# # Merge PCA scores back to the original dataset
# # Create a data frame to store the PCA scores for the rows with no missing data
# birds$BodyCondition <- NA  # Initialize with NA values
# 
# # Add the principal component scores for the rows without NA values
# birds[complete.cases(birds_subset), "BodyCondition"] <- pca_scores[, 1]

### ...create fattening index with PCA -----------------------------------------

# 1. Subset the dataset to only include relevant variables
birds_subset <- birds[, c("Beta", "Tri")]

# 2. Remove rows with NAs
birds_subset_clean <- birds_subset[complete.cases(birds_subset), ]  # n = 89

# 3. Check correlation (optional)
correlation <- cor(birds_subset_clean)
print(correlation)  # -0.23

# 4. Run PCA on the cleaned data (centered and scaled)
pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)

# 5. Check PCA summary (optional)
print(summary(pca_result))

# 6. Extract PCA scores (principal components)
pca_scores <- pca_result$x

# 7. Flip the sign of PC1 so higher scores indicate fattening
pca_scores[, 1] <- -pca_scores[, 1]

# 8. Initialize FatteningIndex column with NA in original dataset
birds$FatteningIndex <- NA

# 9. Add the flipped PC1 scores for rows without missing Beta or Tri
birds[complete.cases(birds_subset), "FatteningIndex"] <- pca_scores[, 1]

# 10. View PCA rotation/loadings (optional)
print(pca_result$rotation)


ggplot(birds_subset_clean, aes(x = Beta, y = Tri)) +
  geom_point(aes(color = pca_scores[, 1])) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  labs(color = "FatteningIndex", title = "PCA-based Fattening Index") +
  theme_minimal()


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






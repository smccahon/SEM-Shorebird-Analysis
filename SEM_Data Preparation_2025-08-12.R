#----------------------------------------#
#          SEM Data Preparation          #
# Created by Shelby McCahon on 8/01/2025 #
#         Modified on 10/21/2025          #
#----------------------------------------#

# load packages
library(tidyverse)
library(purrr)
library(lme4)
library(dplyr)
library(AICcmodavg)
library(lubridate)


# load data
birds <- read.csv("original_data/shorebird_body_condition_data_2025-05-29.csv")
invert <- read.csv("original_data/macroinvertebrate_data_2025-08-22.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")
full <- read.csv("original_data/full_model_dataset_2025-08-29.csv")

# subset to lesser yellowlegs
leye <- subset(birds, Species == "Lesser Yellowlegs") # n = 54

# subset to birds in 2023
leye.2023 <- leye %>%
  filter(Event %in% c("Fall 2023")) # all birds are in the fall 

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

# filter complete cases
birds <- birds %>%
  filter(complete.cases(Mass, Wing))

# fit the log-log linear regression model and compare model fit
# species-corrected
m1 <- lm(log(Mass) ~ log(Wing) + Species, data = birds)

# time-corrected
# m2 <- lm(log(Mass) ~ log(Wing) + Species + time_hours, data = birds)
# m3 <- lm(log(Mass) ~ log(Wing) + Species + Julian, data = birds)
# m4 <- lm(log(Mass) ~ log(Wing) + Species + Event, data = birds)
# m5 <- lm(log(Mass) ~ log(Wing) + Species + time_hours + Julian, data = birds)
# m6 <- lm(log(Mass) ~ log(Wing) + Species + time_hours + Event, data = birds)
# 
# # model comparison
# model_names <- paste0("m", 1:6)
# 
# models <- mget(model_names)
# 
# aictab(models, modnames = model_names) # m4 is best
# 
# # assess model fit
# plot(m4) # good

# extract residuals
# birds <- birds %>%
#   mutate(BCI = residuals(m4))

birds <- birds %>%
  mutate(BCI.NoEvent = residuals(m1))

# plot data
ggplot(birds, aes(x = Species, y = BCI)) + geom_boxplot() +
  theme_classic() +
  labs(x = NULL, y = "Relative Body Condition Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)



### ...create pectoral muscle index (accounting for structural size ) ----------

# does this make any sense?
# filter complete cases
birds.pec <- birds %>%
  filter(complete.cases(PecSizeBest, Wing))

# invert so we can log
birds.pec$PecSizeBest <- birds.pec$PecSizeBest * -1

# fit the log-log linear regression model and compare model fit
# species-corrected
# m1 <- lm(log(PecSizeBest) ~ log(Wing) + Species + Event, data = birds.pec)
# m2 <- lm(log(PecSizeBest) ~ log(Wing) * Species + Event, data = birds.pec)
m3 <- lm(log(PecSizeBest) ~ log(Wing) + Species, data = birds.pec)

# model comparison: is interaction with species needed? absolutely not
# model_names <- paste0("m", 1:3)
# 
# models <- mget(model_names)
# 
# aictab(models, modnames = model_names) 
# 
# plot(m1) # good fit

# extract residuals
# birds.pec <- birds.pec %>%
#   mutate(Standardized.Pec = residuals(m1))

birds.pec <- birds.pec %>%
  mutate(Standardized.Pec.NoEvent = residuals(m3))

# flip sign back
birds.pec$Standardized.Pec.NoEvent <- birds.pec$Standardized.Pec.NoEvent * -1

# plot data
ggplot(birds.pec, aes(x = Species, y = Standardized.Pec.NoEvent)) + geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  labs(x = NULL, y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

ggplot(birds.pec, aes(x = BCI.NoEvent, y = Standardized.Pec.NoEvent)) + 
  geom_point() +
  theme_classic() +
  labs(x = "Mass Index", y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

# invert back 
birds.pec$PecSizeBest <- birds.pec$PecSizeBest * -1

ggplot(birds.pec, aes(x = PecSizeBest, y = Standardized.Pec.NoEvent)) + 
  geom_point() +
  theme_classic() +
  labs(x = "Pectoral Muscle Sizes", y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)


# put back in final dataset
birds <- birds %>%
  left_join(
    birds.pec %>%
      select(Individual, Standardized.Pec.NoEvent),
    by = "Individual"
  )

birds %>% select(Standardized.Pec.NoEvent, PecSizeBest)

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

### ...trim down data and export file for analysis -----------------------------
birds_cleaned <- birds %>% 
  dplyr::select(Individual, Site, Event, Season, Age, Fat, PecSizeBest, Beta, Permanence, 
         NearestCropDistance_m, Max_Flock_Size,
         PlasmaDetection, seconds_since_midnight, Species, PercentAg,,
         Percent_Total_Veg, Julian, Mass, OverallNeonic,
         Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
         AnyDetection, MigStatus, AgCategory, SPEI, Sex, Tri, 
         BCI.NoEvent, Standardized.Pec.NoEvent,
         Diversity, Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
         InvertPesticideDetection, EnvDetection, FatteningIndex, time_hours,
         AnnualSnowfall_in, PrecipitationAmount_7days, 
         DaysSinceLastPrecipitation_5mm)

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
#                         Wetland Data Preparation                          ----                        
#------------------------------------------------------------------------------#   

### ...create water quality index with PCA -------------------------------------

# subset the dataset to only include relevant variables
wetland_subset <- wetland[, c("Conductivity_uS.cm", "Salinity_ppt", 
                            "TDS_mg.L")]

# remove rows with NAs
wetland_subset_clean <- wetland_subset[complete.cases(wetland_subset), ] # n = 100
cor(wetland_subset_clean) # ranges from 0.97-0.99

# Run PCA on the cleaned data
pca_result <- prcomp(wetland_subset_clean, center = TRUE, scale. = TRUE)

# 99% explained by PC1
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
wetland$WaterQuality <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
wetland[complete.cases(wetland_subset), "WaterQuality"] <- pca_scores[, 1]

### ...trim down data and export file for analysis -----------------------------
wetland_cleaned <- wetland %>% 
  select(Wetland, Julian, Season, Year, Event, PercentAg, AgCategory,
         DominantCrop, NearestCropDistance_m, NearestCropType, SPEI,
         Drought.Classification, Buffered, WaterNeonicDetection, 
         NeonicWater_ng.L, NeonicInvert_ng.g, InvertPesticideDetection,
         PesticideInvert_ng.g, EnvDetection, AnnualSnowfall_in, 
         PrecipitationAmount_7days, DaysSinceLastPrecipitation_5mm, Biomass,
         Diversity, pH_probe, WaterTemp, WaterQuality, Permanence,
         Dist_Closest_Wetland_m, PercentLocalVeg_50m, ShorebirdsPresent,
         Wetland.Survey, Invert.Survey, Shorebird.Capture, PlasmaNeonic,
         InvertPesticide, WaterNeonic)

write.csv(wetland_cleaned, "cleaned_data/wetland_data_cleaned_2025-09-30.csv",
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

# fix numeric instability issue in SEM modeling
full$time_hours <- full$seconds_since_midnight / 3600 

### ...create pectoral muscle index (accounting for structural size ) ----------

# does this make any sense?
# filter complete cases
full.pec <- full %>%
  filter(complete.cases(PecSizeBest, Wing))

# invert so we can log
full.pec$PecSizeBest <- full.pec$PecSizeBest * -1

# fit the log-log linear regression model and compare model fit
# species-corrected
m3 <- lm(log(PecSizeBest) ~ log(Wing) + Species, data = full.pec)

# plot(m3) # good fit

# extract residuals
full.pec <- full.pec %>%
  mutate(Standardized.Pec.NoEvent = residuals(m3))

# flip sign back
full.pec$Standardized.Pec.NoEvent <- full.pec$Standardized.Pec.NoEvent * -1

# plot data
ggplot(full.pec, aes(x = Species, y = Standardized.Pec.NoEvent)) + geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  labs(x = NULL, y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

# invert back 
birds.pec$PecSizeBest <- birds.pec$PecSizeBest * -1

# put back in final dataset
full <- full %>%
  left_join(
    full.pec %>%
      select(Individual, Standardized.Pec.NoEvent),
    by = "Individual"
  )

full %>% select(Standardized.Pec.NoEvent, PecSizeBest, Wing)

# plot data
ggplot(full, aes(x = Species, y = Standardized.Pec.NoEvent)) + geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  labs(x = NULL, y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

### ...create body condition index (accounting for structural size ) -----------

# filter complete cases
full <- full %>%
  filter(complete.cases(Mass, Wing))

# fit the log-log linear regression model and compare model fit
# species-corrected
m1 <- lm(log(Mass) ~ log(Wing) + Species, data = full)

# model comparison: is interaction with species needed? absolutely not

# time-corrected
# m2 <- lm(log(Mass) ~ log(Wing) + Species + time_hours, data = full)
# m3 <- lm(log(Mass) ~ log(Wing) + Species + Julian, data = full)
# m4 <- lm(log(Mass) ~ log(Wing) + Species + Event, data = full)
# m5 <- lm(log(Mass) ~ log(Wing) + Species + time_hours + Julian, data = full)
# m6 <- lm(log(Mass) ~ log(Wing) + Species + time_hours + Event, data = full)
# 
# # model comparison
# model_names <- paste0("m", 1:6)
# 
# models <- mget(model_names)
# 
# aictab(models, modnames = model_names) # m1 is best, but use m4 for consistency
# 
# # assess model fit
# plot(m4) # good
# 
# # extract residuals
# full <- full %>%
#   mutate(BCI = residuals(m4))

# extract residuals
 full <- full %>%
   mutate(BCI.NoEvent = residuals(m1))

# plot data
ggplot(full, aes(x = Species, y = BCI.NoEvent)) + geom_boxplot() +
  theme_classic() +
  labs(x = NULL, y = "Relative Body Condition Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

### ...create fattening index with PCA -----------------------------------------
# high tri and low beta for high fattening

# 1. Subset the dataset to only include relevant variables
full_subset <- full[, c("Beta", "Tri")]

# 2. Remove rows with NAs
full_subset_clean <- full_subset[complete.cases(full_subset), ]  # n = 89

# 3. Check correlation (optional)
correlation <- cor(full_subset_clean)
print(correlation)  # -0.23

# 4. Run PCA on the cleaned data (centered and scaled)
pca_result <- prcomp(full_subset_clean, center = TRUE, scale. = TRUE)

# 5. Check PCA summary
print(summary(pca_result))

# 6. Extract PCA scores (principal components)
pca_scores <- pca_result$x

# 7. Flip the sign of PC1 so higher scores indicate fattening
pca_scores[, 1] <- -pca_scores[, 1]

# 8. Initialize FatteningIndex column with NA in original dataset
full$FatteningIndex <- NA

# 9. Add the flipped PC1 scores for rows without missing Beta or Tri
full[complete.cases(full_subset), "FatteningIndex"] <- pca_scores[, 1]

# 10. View PCA rotation/loadings
print(pca_result$rotation)

tri_cor <- cor(pca_scores[, 1], full_subset_clean$Tri)
beta_cor <- cor(pca_scores[, 1], full_subset_clean$Beta)

cat("Correlation with Tri:", round(tri_cor, 3), "\n")   # Should be positive
cat("Correlation with Beta:", round(beta_cor, 3), "\n") # Should be negative

### ...create water quality index with PCA -------------------------------------

# subset the dataset to only include relevant variables
full_subset <- full[, c("Conductivity_uS.cm", "Salinity_ppt", 
                              "TDS_mg.L")]

# remove rows with NAs
full_subset_clean <- full_subset[complete.cases(full_subset), ] # n = 100
cor(full_subset_clean) # ranges from 0.97-0.99

# Run PCA on the cleaned data
pca_result <- prcomp(full_subset_clean, center = TRUE, scale. = TRUE)

# 99% explained by PC1
summary(pca_result)

# View the PCA scores (principal components)
pca_scores <- pca_result$x

# look at loadings
pca_result$rotation

# Merge PCA scores back to the original dataset
# Create a data frame to store the PCA scores for the rows with no missing data
full$WaterQuality <- NA  # Initialize with NA values

# Add the principal component scores for the rows without NA values
full[complete.cases(full_subset), "WaterQuality"] <- pca_scores[, 1]

### ...trim down data and export file for analysis -----------------------------
full_cleaned <- full %>% 
  dplyr::select(Individual, Site, Event, Season, Age, Fat, PecSizeBest, Beta, 
                Permanence, PercentLocalVeg_50m,
                NearestCropDistance_m, Max_Flock_Size,
                PlasmaDetection, seconds_since_midnight, Species, PercentAg, 
                BCI.NoEvent,
                Percent_Total_Veg, Julian, Mass, OverallNeonic,
                Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
                AnyDetection, MigStatus, AgCategory, SPEI, Sex, Tri, 
                PesticideInvert_ng.g, WaterNeonicConc,
                Diversity, Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
                InvertPesticideDetection, EnvDetection, FatteningIndex, time_hours,
                AnnualSnowfall_in, PrecipitationAmount_7days, 
                DaysSinceLastPrecipitation_5mm, Standardized.Pec.NoEvent)

write.csv(full_cleaned, "cleaned_data/full_data_cleaned_2025-10-14.csv",
          row.names = FALSE)


#------------------------------------------------------------------------------#
#                  Lesser Yellowlegs Dataset Preparation (Full)             ----                        
#------------------------------------------------------------------------------#  


### ...standardize time to something more simple -------------------------------
leye$Time <- strptime(leye$Time, format = "%H:%M")
leye$Time <- as.POSIXct(leye$Time, tz = "America/Chicago")
attributes(leye$Time)$tzone

# Calculate seconds since midnight
leye$seconds_since_midnight <- as.numeric(difftime(leye$Time, 
                                                   floor_date(leye$Time, "day"), 
                                                   units = "secs"))

leye$sin <- sin(2 * pi * leye$seconds_since_midnight / (24 * 3600))
leye$cos <- cos(2 * pi * leye$seconds_since_midnight / (24 * 3600))

leye$formatted_time <- format(as.POSIXct(leye$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")

# fix numeric instability issue in SEM modeling
leye$time_hours <- leye$seconds_since_midnight / 3600 

### ...create pectoral muscle index (accounting for structural size ) ----------

# filter complete cases
leye <- leye %>%
  filter(complete.cases(PecSizeBest, Wing))

# invert so we can log
leye$PecSizeBest <- leye$PecSizeBest * -1

# fit the log-log linear regression model
m1 <- lm(log(PecSizeBest) ~ log(Wing), data = leye)

plot(m1) # good fit

# extract residuals
leye <- leye %>%
  mutate(Standardized.Pec = residuals(m1))

# plot data
ggplot(leye, aes(x = Individual, y = Standardized.Pec)) + geom_point() +
  theme_classic() +
  labs(x = NULL, y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

### ...create body condition index (accounting for structural size ) -----------

# filter complete cases
leye <- leye %>%
  filter(complete.cases(Mass, Wing))

# fit the log-log linear regression model and compare model fit
m1 <- lm(log(Mass) ~ log(Wing), data = leye)

# assess model fit
plot(m1) # good

# extract residuals
leye <- leye %>%
  mutate(BCI = residuals(m1))

# plot data
ggplot(leye, aes(x = Individual, y = BCI)) + geom_point() +
  theme_classic() +
  labs(x = NULL, y = "Relative Body Condition Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

# compare BCI and standardized pectoral muscle (highly correlated)
cor(leye$BCI, leye$Standardized.Pec)
plot(leye$BCI, leye$Standardized.Pec)

m <- lm(BCI ~ Standardized.Pec, data = leye)
summary(m)

# compare BCI and mass
cor(leye$BCI, leye$Mass) # very highly correlated

# compare standardized pectoral muscle and mass
cor(leye$Mass, leye$Standardized.Pec) # highly correlated


### ...create fattening index with PCA -----------------------------------------
# high tri and low beta for high fattening

# 1. Subset the dataset to only include relevant variables
leye_subset <- leye[, c("Beta", "Tri")]

# 2. Remove rows with NAs
leye_subset_clean <- leye_subset[complete.cases(leye_subset), ]  # n = 27

# 3. Check correlation (optional)
correlation <- cor(leye_subset_clean)
print(correlation)  # -0.51

# 4. Run PCA on the cleaned data (centered and scaled)
pca_result <- prcomp(leye_subset_clean, center = TRUE, scale. = TRUE)

# 5. Check PCA summary
print(summary(pca_result))

# 6. Extract PCA scores (principal components)
pca_scores <- pca_result$x

# 8. Initialize FatteningIndex column with NA in original dataset
leye$FatteningIndex <- NA

# 9. Add the flipped PC1 scores for rows without missing Beta or Tri
leye[complete.cases(leye_subset), "FatteningIndex"] <- pca_scores[, 1]

# 10. View PCA rotation/loadings
print(pca_result$rotation)

tri_cor <- cor(pca_scores[, 1], leye_subset_clean$Tri)
beta_cor <- cor(pca_scores[, 1], leye_subset_clean$Beta)

cat("Correlation with Tri:", round(tri_cor, 3), "\n")   # Should be positive
cat("Correlation with Beta:", round(beta_cor, 3), "\n") # Should be negative

### ...trim down data and export file for analysis -----------------------------
leye_cleaned <- leye %>% 
  dplyr::select(Individual, Site, Age, Fat, PecSizeBest, Beta, 
                Permanence, 
                NearestCropDistance_m, Max_Flock_Size,
                PlasmaDetection, PercentAg, BCI,
                Percent_Total_Veg, Julian, Mass, OverallNeonic,
                Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
                AnyDetection, AgCategory, SPEI, Sex, Tri,Diversity, 
                Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
                InvertPesticideDetection, EnvDetection, FatteningIndex, time_hours,
                AnnualSnowfall_in, PrecipitationAmount_7days, 
                DaysSinceLastPrecipitation_5mm, Standardized.Pec)

write.csv(leye_cleaned, "cleaned_data/leye_fulldata_cleaned_2025-10-21.csv",
          row.names = FALSE)

#------------------------------------------------------------------------------#
#                  Lesser Yellowlegs Dataset Preparation (2023)             ----                        
#------------------------------------------------------------------------------#  

### ...standardize time to something more simple -------------------------------
leye.2023$Time <- strptime(leye.2023$Time, format = "%H:%M")
leye.2023$Time <- as.POSIXct(leye.2023$Time, tz = "America/Chicago")
attributes(leye.2023$Time)$tzone

# Calculate seconds since midnight
leye.2023$seconds_since_midnight <- as.numeric(difftime(leye.2023$Time, 
                                                   floor_date(leye.2023$Time, "day"), 
                                                   units = "secs"))

leye.2023$sin <- sin(2 * pi * leye.2023$seconds_since_midnight / (24 * 3600))
leye.2023$cos <- cos(2 * pi * leye.2023$seconds_since_midnight / (24 * 3600))

leye.2023$formatted_time <- format(as.POSIXct(leye.2023$seconds_since_midnight, 
                                         origin = "1970-01-01", tz = "UTC"), 
                              "%H:%M")

# fix numeric instability issue in SEM modeling
leye.2023$time_hours <- leye.2023$seconds_since_midnight / 3600 

### ...create pectoral muscle index (accounting for structural size ) ----------

# filter complete cases
leye.2023 <- leye.2023 %>%
  filter(complete.cases(PecSizeBest, Wing))

# invert so we can log
leye.2023$PecSizeBest <- leye.2023$PecSizeBest * -1

# fit the log-log linear regression model and compare model fit
# species-corrected
m1 <- lm(log(PecSizeBest) ~ log(Wing), data = leye.2023)

plot(m1) # okay fit

# extract residuals
leye.2023 <- leye.2023 %>%
  mutate(Standardized.Pec = residuals(m1))

# plot data
ggplot(leye.2023, aes(x = Individual, y = Standardized.Pec)) +
  geom_point() +
  theme_classic() +
  labs(x = NULL, y = "Relative Pectoral Muscle Size Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

### ...create body condition index (accounting for structural size ) -----------

# filter complete cases
leye.2023 <- leye.2023 %>%
  filter(complete.cases(Mass, Wing))

# fit the log-log linear regression model and compare model fit
m1 <- lm(log(Mass) ~ log(Wing), data = leye.2023)

# assess model fit
plot(m1) # good

# extract residuals
leye.2023 <- leye.2023 %>%
  mutate(BCI = residuals(m1))

# plot data
ggplot(leye.2023, aes(x = Individual, y = BCI)) + geom_point() +
  theme_classic() +
  labs(x = NULL, y = "Relative Body Condition Index") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  geom_hline(linetype = "dashed", color = "red", yintercept = 0,
             size = 1)

# compare BCI and standardized pectoral muscle
cor(leye.2023$BCI, leye.2023$Standardized.Pec) # highly correlated
plot(leye.2023$BCI, leye.2023$Standardized.Pec)

m <- lm(BCI ~ Standardized.Pec, data = leye.2023)
summary(m)

# compare BCI and mass
cor(leye.2023$BCI, leye.2023$Mass) # extremely correlated
cor(leye.2023$Standardized.Pec, leye.2023$Mass) # highly correlated


### ...create fattening index with PCA -----------------------------------------
# high tri and low beta for high fattening

# 1. Subset the dataset to only include relevant variables
leye.2023_subset <- leye.2023[, c("Beta", "Tri")]

# 2. Remove rows with NAs
leye.2023_subset_clean <- leye.2023_subset[complete.cases(leye.2023_subset), ]  # n = 89

# 3. Check correlation (optional)
correlation <- cor(leye.2023_subset_clean)
print(correlation)  # -0.73

# 4. Run PCA on the cleaned data (centered and scaled)
pca_result <- prcomp(leye.2023_subset_clean, center = TRUE, scale. = TRUE)

# 5. Check PCA summary
print(summary(pca_result))

# 6. Extract PCA scores (principal components)
pca_scores <- pca_result$x

# 7. Flip the sign of PC1 so higher scores indicate fattening
pca_scores[, 1] <- -pca_scores[, 1]

# 8. Initialize FatteningIndex column with NA in original dataset
leye.2023$FatteningIndex <- NA

# 9. Add the flipped PC1 scores for rows without missing Beta or Tri
leye.2023[complete.cases(leye.2023_subset), "FatteningIndex"] <- pca_scores[, 1]

# 10. View PCA rotation/loadings
print(pca_result$rotation)

tri_cor <- cor(pca_scores[, 1], leye.2023_subset_clean$Tri)
beta_cor <- cor(pca_scores[, 1], leye.2023_subset_clean$Beta)

cat("Correlation with Tri:", round(tri_cor, 3), "\n")   # Should be positive
cat("Correlation with Beta:", round(beta_cor, 3), "\n") # Should be negative

### ...trim down data and export file for analysis -----------------------------
leye.2023_cleaned <- leye.2023 %>% 
  dplyr::select(Individual, Site, Age, Fat, PecSizeBest, Beta, 
                Permanence, 
                NearestCropDistance_m, Max_Flock_Size,
                PlasmaDetection, PercentAg, BCI,
                Percent_Total_Veg, Julian, Mass, OverallNeonic,
                Uric, Biomass, DominantCrop, NearestCropType, WaterNeonicDetection, 
                AnyDetection, AgCategory, SPEI, Sex, Tri,Diversity, 
                Dist_Closest_Wetland_m, Percent_Exposed_Shoreline,
                InvertPesticideDetection, EnvDetection, FatteningIndex, time_hours,
                AnnualSnowfall_in, PrecipitationAmount_7days, 
                DaysSinceLastPrecipitation_5mm, Standardized.Pec)

write.csv(leye.2023_cleaned, "cleaned_data/leye_2023data_cleaned_2025-10-21.csv",
          row.names = FALSE)












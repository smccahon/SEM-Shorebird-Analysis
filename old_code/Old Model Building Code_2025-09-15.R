#----------------------------------------#
#        Original Code Documentation     #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/17/2025          #
#----------------------------------------#


#------------------------------------------------------------------------------#
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

#------------------------------------------------------------------------------#
#               initial body condition index calculation (9/15/2025)        ----
#------------------------------------------------------------------------------#

# ### ...create body condition index 
# ### (NOT accounting for structural size or Species!!)

# # subset the dataset to only include relevant variables
 birds_subset <- birds[, c("Fat", "Mass", "PecSizeBest")]
# 
# # remove rows with NAs
 birds_subset_clean <- birds_subset[complete.cases(birds_subset), ] # n = 152
 cor(birds_subset_clean) # ranges from -0.13 to 0.36
# 
# # Run PCA on the cleaned data
 pca_result <- prcomp(birds_subset_clean, center = TRUE, scale. = TRUE)
# 
# # 52% explained by PC1, 29% explained by PC2, 19% explained by PC3
 summary(pca_result)
# 
# # View the PCA scores (principal components)
 pca_scores <- pca_result$x
# 
# # Merge PCA scores back to the original dataset
# # Create a data frame to store the PCA scores for the rows with no missing data
 birds$BodyCondition <- NA  # Initialize with NA values
# 
# # Add the principal component scores for the rows without NA values
 birds[complete.cases(birds_subset), "BodyCondition"] <- pca_scores[, 1]
 
 #------------------------------------------------------------------------------#
 
 ### ...accounting for structural size but I don't think this is right
 
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
 
 #------------------------------------------------------------------------------#
 #                    Full Dataset Preparation (9/11/2025)                  ----                        
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
 ### NOT ACCOUNTING FOR SPECIES OR SIZE
 
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
 ### NOT ACCOUNT FOR SPECIES 
 
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
 
 ### ...extract standardized coefficients manually (long code)-----------------------------
 
 m2 <- glm(Biomass ~ PercentAg, data = invert.pos, family = Gamma(link = "log"))
 beta <- coef(m2)["PercentAg"] # get raw unstandardized coefficient
 phi <- summary(m2)$dispersion # extract dispersion parameter
 v <- 1/phi # calculate shape parameter
 sigma2_E <- trigamma(v) # distribution-specific error variance
 sd_x <- sd(invert.pos$PercentAg) # calculate SD of predictor
 sd_y <- sqrt(sigma2_E + var(predict(m2, type = "link"))) # calculate SD of response on link scale
 beta_std <- beta * (sd_x / sd_y) # compute standardized coefficient
 beta_std
 
 
 
 

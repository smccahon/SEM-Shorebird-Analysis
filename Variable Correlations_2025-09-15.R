#----------------------------------------#
#          Variable Correlations         #
# Created by Shelby McCahon on 9/15/2025 #
#         Modified on 9/15/2025          #
#----------------------------------------#

# load data
birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
invert <- read.csv("cleaned_data/invert_data_cleaned_2025-08-11.csv")
wetland <- read.csv("original_data/neonic_wetland_survey_data_2025-08-12.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-08-29.csv")

# filter to only include 2023 data
birds <- birds %>% 
  filter(Event %in% c("Fall 2023", "Spring 2023")) # 126 birds

wetland <- wetland %>% 
  filter(Year == "2023") # 79 wetland surveys surveys

#------------------------------------------------------------------------------#
#                        convert factors to numeric                         ----                        
#------------------------------------------------------------------------------# 

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


invert <- invert %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0),
    Buffered = case_when(
      Buffered == "Y" ~ 1,
      Buffered == "N" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))


wetland <- wetland %>% 
  mutate(WaterNeonicDetection = case_when(
    WaterNeonicDetection == "Y" ~ 1,
    WaterNeonicDetection == "N" ~ 0,
    TRUE ~ NA_real_),
    InvertPesticideDetection = case_when(
      InvertPesticideDetection == "Y" ~ 1,
      InvertPesticideDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    EnvDetection = case_when(
      EnvDetection == "Y" ~ 1,
      EnvDetection == "N" ~ 0,
      TRUE ~ NA_real_),
    Season = case_when(
      Season == "Spring" ~ 1,
      Season == "Fall" ~ 0),
    Buffered = case_when(
      Buffered == "Y" ~ 1,
      Buffered == "N" ~ 0,
      TRUE ~ NA_real_),
    Permanence = case_when(
      Permanence %in% c("Temporary", "Seasonal") ~ 1,
      Permanence == "Semipermanent" ~ 2,
      Permanence == "Permanent" ~ 3,
      TRUE ~ NA_real_))

full <- full %>% 
  mutate(PlasmaDetection = case_when(
    PlasmaDetection == "Y" ~ 1,
    PlasmaDetection == "N" ~ 0,
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

#------------------------------------------------------------------------------#
#                    continuous variable correlations                       ----                        
#------------------------------------------------------------------------------# 

# Select only numeric columns
numeric_vars <- invert %>% select(where(is.numeric))

# Calculate correlation matrix (Pearson by default)
cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")

# Set a correlation threshold
threshold <- 0.60

# Get upper triangle of correlation matrix without diagonal
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA

# Find pairs with absolute correlation >= threshold
strong_corr <- which(abs(cor_matrix) >= threshold, arr.ind = TRUE)

# Format output nicely
result <- data.frame(
  Var1 = rownames(cor_matrix)[strong_corr[, 1]],
  Var2 = colnames(cor_matrix)[strong_corr[, 2]],
  Correlation = cor_matrix[strong_corr]
)

print(result)

#------------------------------------------------------------------------------#
#              categorical variable correlations (INVERTS)                  ----                        
#------------------------------------------------------------------------------# 


# dominant crop type and % cropland cover
summary(aov(data = invert, PercentAg ~ DominantCrop)) # yes

# season and permanence
summary(aov(data = invert, Season ~ Permanence)) # yes

# season and detection
summary(aov(data = invert, Season ~ EnvDetection)) # no

# season and percentag
summary(aov(data = invert, Season ~ PercentAg)) # no

# julian and percentag
summary(aov(data = invert, Julian ~ Permanence)) # yes

# julian and envdetection
summary(aov(data = invert, Julian ~ EnvDetection)) # no

# julian and envdetection
summary(aov(data = invert, PercentAg ~ EnvDetection)) # no

# percentveg and envdetection
summary(aov(data = invert, PercentLocalVeg_50m ~ EnvDetection)) # no

# percentveg and envdetection
summary(aov(data = invert, PercentLocalVeg_50m ~ EnvDetection)) # no

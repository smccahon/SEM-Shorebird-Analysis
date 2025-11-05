#-----------------------------------------#
#         Dataset Investigation           #
# Created by Shelby McCahon on 11/04/2025 #
#         Modified on 11/04/2025          #
#-----------------------------------------#

# load packages
library(piecewiseSEM)
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)
library(multcompView)
library(DHARMa)
library(car)
library(AICcmodavg)
library(statmod)
library(MuMIn)

#------------------------------------------------------------------------------#
#                        load data and organize datasets                    ----                        
#------------------------------------------------------------------------------# 

birds <- read.csv("cleaned_data/shorebird_data_cleaned_2025-08-11.csv")
full <- read.csv("cleaned_data/full_data_cleaned_2025-10-14.csv")

# removed a short-billed dowitcher and a greater yellowlegs
full <- full %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# removed four birds
birds <- birds %>% 
  group_by(Species) %>% 
  filter(n() >= 3) %>% 
  ungroup()

# theme for plotting
my_theme <- theme_classic() + theme(
  axis.title.x = element_text(size = 21, margin = margin(t = 12)),
  axis.title.y = element_text(size = 21, margin = margin(r = 12)),
  axis.text.x = element_text(size = 18),
  axis.text.y = element_text(size = 18))

options(tibble.print_max = Inf)
options(digits = 3)

#------------------------------------------------------------------------------#
#                             data patterns                                 ----                        
#------------------------------------------------------------------------------# 

# how much does fattening index and % ag vary by season/event
tapply(full.FI$PercentAg, full.FI$Season, mean)
#    0    1  
# 52.1 63.1 

tapply(full.FI$FatteningIndex, full.FI$Season, mean, na.rm = TRUE)
#   Fall  Spring 
# -0.207  0.390 

tapply(birds.FI$PercentAg, birds.FI$Event, mean)
# Fall 2021   Fall 2023   Spring 2022 Spring 2023 
# 47.5        52.1        47.6        63.1 

tapply(birds.FI$FatteningIndex, birds.FI$Event, mean, na.rm = TRUE)
# Fall 2021   Fall 2023 Spring 2022 Spring 2023 
#   -0.0312     -0.1458     -0.3626      0.4885 

tapply(birds.FI$FatteningIndex, birds.FI$PlasmaDetection, mean, na.rm = TRUE)
# N       Y 
# 0.0965 -0.2453 

tapply(full.FI$FatteningIndex, full.FI$PlasmaDetection, mean, na.rm = TRUE)
# N      Y 
# 0.018 -0.164 

table(birds.FI$Event, birds.FI$PlasmaDetection)
#              N  Y
# Fall 2021    2  1
# Fall 2023   39  8
# Spring 2022  0 14
# Spring 2023 17  5

table(full.FI$Event, full.FI$PlasmaDetection)
#              N  Y
# Fall 2023   39  8
# Spring 2023 17  5


# does the effect of % ag vary by season/event? answer is no
m2 <- glmmTMB(FatteningIndex ~ PercentAg * Event + (1|Species), 
              data = birds)

m2 <- glmmTMB(FatteningIndex ~ PercentAg * Season + (1|Species), 
              data = full)

summary(m2)


# why is there a negative effect of surrounding cropland on fattening index
# only for 2023 birds?
# I think it's because in spring 2023, fattening index was the highest (mean = 0.49)
# and the mean percent ag was also the highest (mean = 70%). When you add in 
# more years, this trend is not as consistent

m1 <- glmmTMB(FatteningIndex ~ PercentAg + (1|Species), 
         data = full)
m2 <- glmmTMB(FatteningIndex ~ PercentAg + (1|Species), 
           data = birds)

summary(m1)
summary(m2)


# no effect
m1 <- glmmTMB(FatteningIndex ~ PercentAg + (1|Species), 
              data = full)

# after accounting for season, there is now an effect of % ag
m1 <- glmmTMB(FatteningIndex ~ PercentAg + Season +
                (1|Species), 
              data = full)

summary(m1)
cor(full$Season, full$PercentAg)


m1 <- glmmTMB(FatteningIndex ~ PercentAg + Event +
                (1|Species), 
              data = birds)
summary(m1)

# view relationships
ggplot(full, aes(x = PercentAg, y = FatteningIndex)) +
  geom_point() + my_theme + 
  geom_hline(yintercept = 0,
             color = "red",
             size = 1,
             linetype = "dashed") +
  labs(x = "% Surrounding Cropland",
       y = "Fattening Index (2023 birds)")

ggplot(birds, aes(x = PercentAg, y = FatteningIndex)) +
  geom_point() + my_theme + 
  geom_hline(yintercept = 0,
             color = "red",
             size = 1,
             linetype = "dashed") +
  labs(x = "% Surrounding Cropland",
       y = "Fattening Index (All birds)")

#-------------------------------------------------------------#
# what's special about the 17 birds from 2021 and 2022?   -----
#-------------------------------------------------------------#

table(birds$FatteningIndex, birds$Event)

#----------------------#
#      FALL 2021       #
#----------------------#

# 3 birds
# two low (Lesser Yellowlegs), one high (Pectoral Sandpiper)
fall.FI <- birds %>% 
  filter(Event == "Fall 2021") %>% 
  filter(!is.na(FatteningIndex))

# mostly all within the same date
# both low values are from a Lesser Yellowlegs
# body condition and fat higher in PESA

fall.FI %>% 
  select(Species, FatteningIndex, Julian,
         Mass, SPEI, PercentAg, PlasmaDetection,
         BCI.NoEvent, Fat)

# Species FatteningIndex Julian  Mass  SPEI PercentAg PlasmaDetection BCI.NoEvent   Fat
# LEYE            -0.359    199  85.2 -1.56      3.75 N                   -0.0538     0
# PESA             0.596    202  91.1 -1.75     65.8  N                    0.293      2
# LEYE            -0.331    207  98.2 -1.51     73.0  Y                    0.0419     0
# 

ggplot(fall.FI, aes(x = Individual, y = FatteningIndex)) +
  geom_point() + geom_hline(yintercept = 0,
                            color = "red",
                            linetype = "dashed") +
  labs(x = "Individual (Fall 2021)",
       y = "Fattening Index") + my_theme

#----------------------#
#      SPRING 2022     #
#----------------------#

# 14 birds
spring.FI <- birds %>% 
  filter(Event == "Spring 2022") %>% 
  filter(!is.na(FatteningIndex))


spring.FI %>% 
  select(Species, FatteningIndex, Julian,
         Mass, SPEI, PercentAg, PlasmaDetection,
         BCI.NoEvent, Fat)

# lots of WIPH
# likely not due to fat, Julian day, drought, or plasma detection
# 10/14 birds have FI lower than 0


# Species     FatteningIndex Julian  Mass  SPEI PercentAg PlasmaDetection BCI.NoEvent   Fat
# LEYE                -0.149    116 110.0  1.00     46.3  Y                     0.114     0
# LEYE                 0.258    117  81.3  1.00     44.6  Y                    -0.154     0
# LEYE                 -1.02    117  88.9  1.00     44.6  Y                   -0.0269     0
# WIPH                  1.09    120  55.7  2.22     60.0  Y                   -0.0194     0
# WIPH                  1.26    120  55.8  2.22     60.0  Y                   -0.0363     0
# WIPH                -0.847    121  49.2  2.22     27.7  Y                    -0.143     0
# WIPH                -0.790    121  69.9  2.22     27.7  Y                     0.152     3
# SESA                -0.241    121  27.9  2.22     27.7  Y                     0.138     0
# SESA                -0.242    121  26.2  2.22     27.7  Y                    0.0882     0
# WIPH                -0.501    122  67.6  2.22     60.0  Y                     0.101     0
# WIPH                 -1.31    123  62.4  2.22     60.0  Y                    0.0571     0
# WIPH                 0.728    123  58.1  2.22     60.0  Y                   -0.0234     0
# WIPH                -0.621    123  60.2  2.22     60.0  Y                   -0.0236     0
# WIPH                 -2.69    123  62.1  2.22     60.0  Y                    0.0615     0

ggplot(spring.FI, aes(x = Individual, y = FatteningIndex,
                    color = PercentAg)) +
  geom_point() + geom_hline(yintercept = 0,
                            color = "red",
                            linetype = "dashed") +
  labs(x = "Individual (Fall 2021)",
       y = "Fattening Index") + my_theme

m <- lm(FatteningIndex ~ PercentAg, 
        data = spring.FI)

summary(m) # no relationship


#------------------------------#
# BIRDS FROM ANALYSIS      -----
#------------------------------#

# 69 birds
full %>% 
  count(!is.na(FatteningIndex)) %>% 
  print()

full.FI <- full %>% 
  filter(!is.na(FatteningIndex))

# only 13 birds with detections
table(full.FI$PlasmaDetection)


# 86 birds
birds %>% 
  count(!is.na(FatteningIndex)) %>% 
  print()

birds.FI <- birds %>% 
  filter(!is.na(FatteningIndex))

# 28 birds with detections
table(birds.FI$PlasmaDetection)

cols <- c(
  "Fall 2023" = "plum",  
  "Fall 2021"   = "darkorange3",  
  "Spring 2022" = "#56B4E9",  
  "Spring 2023"   = "#5DC863"   
)

# 2021-2023 birds dataset
ggplot(full.FI, aes(x = PercentAg, y = FatteningIndex, 
                     colour = Event)) +
  geom_point(size = 1.75) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_colour_manual(values = cols) +
  labs(
    x = "% Surrounding Cropland",
    y = "Fattening Index",
    colour = "Sampling Event"
  ) +
  my_theme +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )



# overlapped dataset
ggplot(full.FI, aes(x = PercentAg, y = FatteningIndex,
                     color = Event)) +
  geom_point() + geom_hline(yintercept = 0,
                            color = "red",
                            linetype = "dashed") +
  labs(x = "% Surrounding Cropland",
       y = "Fattening Index") + my_theme +
  theme(legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

        
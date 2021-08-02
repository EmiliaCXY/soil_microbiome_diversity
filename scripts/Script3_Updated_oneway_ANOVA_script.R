# March 30, 2021
# Updated script for One-way ANOVA
# Used updated metadata file: P2_MICB_421_Soil_Metadata.tsv 

# Load packages
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)

# Load soil metadata
P2_soil_metadata <- read_tsv(file = 'P2_MICB_421_Soil_Metadata.tsv')

# Filter soil metadata to only include 'O horizon' samples that have pH > 0
P2_filtered_soil_metadata <- filter(P2_soil_metadata, Region == "British Columbia", Horizon == "O horizon", pH > 0)

# Rename columns for "Moisture Content" metadata column to Moisture_Content (we don't want any spaces)
colnames(P2_filtered_soil_metadata)[colnames(P2_filtered_soil_metadata) %in% c("Moisture Content")] <- c("Moisture_Content")

# Find mean and standard deviation of moisture content for each precipitation group
P2_filtered_soil_metadata %>%
  group_by(Mean_Annual_Precipitation_mm) %>%
  summarise(mean = mean(Moisture_Content), 
            sd = sd(Moisture_Content))

# Make boxplots to look at moisture content distribution for each precipitation group
boxplot(Moisture_Content ~ Mean_Annual_Precipitation_mm, 
        data = P2_filtered_soil_metadata,
        xlab = "Mean Annual Precipitation (mm)",
        ylab = "Moisture Content",
        col = "steelblue",
        border = "black")

# One-way anova for moisture content and precipitation
one.way.aov_moisture_precipitation <- aov(Moisture_Content ~ Mean_Annual_Precipitation_mm, data = P2_filtered_soil_metadata)
summary(one.way.aov_moisture_precipitation)

# Find mean and standard deviation of moisture content for each compaction treatment
P2_filtered_soil_metadata %>%
  group_by(Compaction_Treatment) %>%
  summarise(mean = mean(Moisture_Content), 
            sd = sd(Moisture_Content))

# Make boxplots to look at moisture content distribution for each compaction treatment
boxplot(Moisture_Content ~ Compaction_Treatment, 
        data =P2_filtered_soil_metadata,
        xlab = "Compaction Treatment",
        ylab = "Moisture Content",
        col = "steelblue",
        border = "black")

# One-way anova for moisture content and compaction
one.way.aov_moisture_compaction <- aov(Moisture_Content ~ Compaction_Treatment, data = P2_filtered_soil_metadata)
summary(one.way.aov_moisture_compaction)

# Find mean and standard deviation of moisture content for each OM removal group
P2_filtered_soil_metadata %>%
  group_by(OM_Removal) %>%
  summarise(mean = mean(Moisture_Content), 
            sd = sd(Moisture_Content))

# Make boxplots to look at moisture content distribution for each OM removal group
boxplot(Moisture_Content ~ OM_Removal, 
        data =P2_filtered_soil_metadata,
        xlab = "LTSP Treatment",
        ylab = "Moisture Content",
        col = "steelblue",
        border = "black")

# One-way anova for moisture content and compaction
one.way.aov_moisture_OM_removal <- aov(Moisture_Content ~ OM_Removal, data = P2_filtered_soil_metadata)
summary(one.way.aov_moisture_OM_removal)


# Find mean and standard deviation of moisture content for each moisture bin
P2_filtered_soil_metadata %>%
  group_by(Moisture_bin) %>%
  summarise(mean = mean(Moisture_Content), 
            sd = sd(Moisture_Content))

# Make boxplots to look at moisture content distribution for each moisture bin
boxplot(Moisture_Content ~ Moisture_bin, 
        data =P2_filtered_soil_metadata,
        xlab = "Moisture Bin",
        ylab = "Moisture Content",
        col = "steelblue",
        border = "black") 

# One-way anova for moisture content and compaction
one.way.aov_moisture_unbinned_binned <- aov(Moisture_Content ~ Moisture_bin, data = P2_filtered_soil_metadata)
summary(one.way.aov_moisture_unbinned_binned)



# MICB421 Project 2 Soil - OM Removal Differential abundance analysis
# April 12, 2021
# Jessica Nguyen, Team 5

# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages
library(phyloseq)
library(DESeq2)

# Calculate relative abundance
calculate_relative_abundance <- function(x) x / sum(x)

# Calculate geometric mean 
calculate_gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Define the set of random numbers
set.seed(711)

# Import qiime files
biom_file <- import_biom("table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("MICB_421_Soil_Metadata.tsv")
tree <- read_tree_greengenes("tree.nwk")

# Convert tree from multichotomous to dichotomous tree
tree <- multi2di(tree)

# Combine information into single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)
physeq

# Look at metadata, head shows only the first 6 lines
head(sample_data(physeq))

# Look at and save data as data frame
metadata_df <- as.data.frame(sample_data(physeq))

# Change taxonomic levels from numbers to actual names
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
head(tax_table(physeq))

# Differential abundance analysis: OM1 vs REF
# Step 1: Remove OM2 treatment so that only OM1 and REF treatments remain
OM1_REF_physeq <- subset_samples(physeq, LTSP.Treatment != "OM2")
OM1_REF_physeq

# Step 2: Determine the total counts and relative abundance for features across all OM1 and REF samples
OM1_REF_total_counts <- taxa_sums(OM1_REF_physeq)
OM1_REF_relative_abundance <- calculate_relative_abundance(OM1_REF_total_counts)

# Remove low abundant features
OM1_REF_abundant <- OM1_REF_relative_abundance > 0.0005
OM1_REF_abundant_taxa <- prune_taxa(OM1_REF_abundant, OM1_REF_physeq)
OM1_REF_abundant_taxa

# Set taxonomic level for analysis (Genus)
OM1_REF_genera <- tax_glom(OM1_REF_abundant_taxa, taxrank = "Genus")
OM1_REF_genera

# Create a DESeq2 object for OM1 and REF
deseq_OM1_REF <- phyloseq_to_deseq2(OM1_REF_genera, ~ LTSP.Treatment)
geo_means_OM1_REF <- apply(counts(deseq_OM1_REF), 1, calculate_gm_mean)
deseq_OM1_REF <- estimateSizeFactors(deseq_OM1_REF, geoMeans = geo_means_OM1_REF)
deseq_OM1_REF <- DESeq(deseq_OM1_REF, fitType = "local")

# Look at DESeq2 results for OM1 vs REF
OM1_REF_diff_abund <- results(deseq_OM1_REF)
OM1_REF_no_alpha_diff_abund <- as.data.frame(OM1_REF_diff_abund)
dim(OM1_REF_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_OM1_REF <- as.data.frame(OM1_REF_diff_abund)
sig_OM1_REF <- filter(sig_OM1_REF, padj < alpha)
dim(sig_OM1_REF)

# sig_OM1_REF does not show any significant differences differential abundance between OM1 and REF
# Output: [1] 0 6
# Use OM1_REF_no_alpha_diff_abund for the next steps instead

# Merge OM1 vs REF DESeq2 results with table of taxonomic info (OM1_REF_genera)
genera_df <- as.data.frame(tax_table(OM1_REF_genera))
OM1_REF_no_alpha_diff_abund <- merge(OM1_REF_no_alpha_diff_abund, genera_df, by = "row.names")
OM1_REF_no_alpha_diff_abund <- arrange(OM1_REF_no_alpha_diff_abund, log2FoldChange)
dim(OM1_REF_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_OM1_REF_no_alpha_diff_abund <- filter(OM1_REF_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_OM1_REF_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_OM1_REF_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()




# Differential abundance analysis: OM2 vs REF
# Remove OM1 treatment so that only OM1 and REF treatments remain
OM2_REF_physeq <- subset_samples(physeq, LTSP.Treatment != "OM1")
OM2_REF_physeq

# Determine the total counts and relative abundance for features across all OM2 and REF samples
OM2_REF_total_counts <- taxa_sums(OM2_REF_physeq)
OM2_REF_relative_abundance <- calculate_relative_abundance(OM2_REF_total_counts)

# Remove low abundant features
OM2_REF_abundant <- OM2_REF_relative_abundance > 0.0005
OM2_REF_abundant_taxa <- prune_taxa(OM2_REF_abundant, OM2_REF_physeq)
OM2_REF_abundant_taxa

# Set taxonomic level for analysis (Genus)
OM2_REF_genera <- tax_glom(OM2_REF_abundant_taxa, taxrank = "Genus")
OM2_REF_genera

# Create a DESeq2 object for OM2 and REF
deseq_OM2_REF <- phyloseq_to_deseq2(OM2_REF_genera, ~ LTSP.Treatment)
geo_means_OM2_REF <- apply(counts(deseq_OM2_REF), 1, calculate_gm_mean)
deseq_OM2_REF <- estimateSizeFactors(deseq_OM2_REF, geoMeans = geo_means_OM2_REF)
deseq_OM2_REF <- DESeq(deseq_OM2_REF, fitType = "local")

# Look at DESeq2 results for OM2 vs REF and save as a data frame
OM2_REF_diff_abund <- results(deseq_OM2_REF)
OM2_REF_no_alpha_diff_abund <- as.data.frame(OM2_REF_diff_abund)
dim(OM2_REF_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_OM2_REF <- as.data.frame(OM2_REF_diff_abund)
sig_OM2_REF <- filter(sig_OM2_REF, padj < alpha)
dim(sig_OM2_REF)

# sig_OM2_REF does not show any significant differences differential abundance between OM2 and REF
# Output: [1] 0 6
# Use OM2_REF_no_alpha_diff_abund for the next steps instead

# Merge OM2 vs REF DESeq2 results with table of taxonomic info (OM2_REF_genera)
genera_OM2_REF_df <- as.data.frame(tax_table(OM2_REF_genera))
OM2_REF_no_alpha_diff_abund <- merge(OM2_REF_no_alpha_diff_abund, genera_OM2_REF_df, by = "row.names")
OM2_REF_no_alpha_diff_abund <- arrange(OM2_REF_no_alpha_diff_abund, log2FoldChange)
dim(OM2_REF_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_OM2_REF_no_alpha_diff_abund <- filter(OM2_REF_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_OM2_REF_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_OM2_REF_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()



# Differential abundance analysis: OM1 vs OM2
# Step 1: Remove REF treatment so that only OM1 and OM2 treatments remain
OM1_2_physeq <- subset_samples(physeq, LTSP.Treatment != "REF")
OM1_2_physeq

# Step 2: Determine the total counts and relative abundance for features across all OM1 and OM2 samples
OM1_2_total_counts <- taxa_sums(OM1_2_physeq)
OM1_2_relative_abundance <- calculate_relative_abundance(OM1_2_total_counts)

# Remove low abundant features
OM1_2_abundant <- OM1_2_relative_abundance > 0.0005
OM1_2_abundant_taxa <- prune_taxa(OM1_2_abundant, OM1_2_physeq)
OM1_2_abundant_taxa

# Set taxonomic level for analysis (Genus)
OM1_2_genera <- tax_glom(OM1_2_abundant_taxa, taxrank = "Genus")
OM1_2_genera

# Create a DESeq2 object for OM1 and OM2
deseq_OM1_2 <- phyloseq_to_deseq2(OM1_2_genera, ~ LTSP.Treatment)
geo_means_OM1_2 <- apply(counts(deseq_OM1_2), 1, calculate_gm_mean)
deseq_OM1_2 <- estimateSizeFactors(deseq_OM1_2, geoMeans = geo_means_OM1_2)
deseq_OM1_2 <- DESeq(deseq_OM1_2, fitType = "local")

# Look at DESeq2 results for OM1 vs OM2 and save as a data frame
OM1_2_diff_abund <- results(deseq_OM1_2)
OM1_2_no_alpha_diff_abund <- as.data.frame(OM1_2_diff_abund)
dim(OM1_2_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_OM1_2 <- as.data.frame(OM1_2_diff_abund)
sig_OM1_2 <- filter(sig_OM1_2, padj < alpha)
dim(sig_OM1_2)

# sig_OM1_OM2 does not show any significant differences differential abundance between OM1 and OM2
# Output: [1] 0 6
# Use OM1_2_no_alpha_diff_abund for the next steps instead

# Merge OM2 vs REF DESeq2 results with table of taxonomic info (OM2_REF_genera)
genera_OM1_2_df <- as.data.frame(tax_table(OM1_2_genera))
OM1_2_no_alpha_diff_abund <- merge(OM1_2_no_alpha_diff_abund, genera_OM1_2_df, by = "row.names")
OM1_2_no_alpha_diff_abund <- arrange(OM1_2_no_alpha_diff_abund, log2FoldChange)
dim(OM1_2_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_OM1_2_no_alpha_diff_abund <- filter(OM1_2_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_OM1_2_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_OM1_2_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()



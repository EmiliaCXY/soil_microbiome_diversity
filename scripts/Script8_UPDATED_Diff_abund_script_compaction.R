# MICB421 Project 2 Soil - Compaction Treatment Differential abundance analysis
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



# Differential abundance analysis: C0 vs REF
# Remove C1 treatment so that only C0, C2, and REF treatments remain
C0_C2_REF_physeq <- subset_samples(physeq, Compaction_Treatment != "C1")
C0_C2_REF_physeq

# Remove C2 so that only C0 and REF treatments remain
C0_REF_physeq <- subset_samples(C0_C2_REF_physeq, Compaction_Treatment != "C2")
C0_REF_physeq

# Determine the total counts and relative abundance for features across all C0 and REF samples
C0_REF_total_counts <- taxa_sums(C0_REF_physeq)
C0_REF_relative_abundance <- calculate_relative_abundance(C0_REF_total_counts)

# Remove low abundant features
C0_REF_abundant <- C0_REF_relative_abundance > 0.0005
C0_REF_abundant_taxa <- prune_taxa(C0_REF_abundant, C0_REF_physeq)
C0_REF_abundant_taxa

# Set taxonomic level for analysis (Genus)
C0_REF_genera <- tax_glom(C0_REF_abundant_taxa, taxrank = "Genus")
C0_REF_genera

# Create a DESeq2 object for C0 and REF
deseq_C0_REF <- phyloseq_to_deseq2(C0_REF_genera, ~ Compaction_Treatment)
geo_means_C0_REF <- apply(counts(deseq_C0_REF), 1, calculate_gm_mean)
deseq_C0_REF <- estimateSizeFactors(deseq_C0_REF, geoMeans = geo_means_C0_REF)
deseq_C0_REF <- DESeq(deseq_C0_REF, fitType = "local")

# Look at DESeq2 results for C0 vs REF
C0_REF_diff_abund <- results(deseq_C0_REF)
C0_REF_no_alpha_diff_abund <- as.data.frame(C0_REF_diff_abund)
dim(C0_REF_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C0_REF <- as.data.frame(C0_REF_diff_abund)
sig_C0_REF <- filter(sig_C0_REF, padj < alpha)
dim(sig_C0_REF)

# sig_C0_REF does not show any significant differences differential abundance between C0 and REF
# Output: [1] 0 6
# Use C0_REF_no_alpha_diff_abund for the next steps instead

# Merge C0 vs REF DESeq2 results with table of taxonomic info (C0_REF_genera)
genera_C0_REF_df <- as.data.frame(tax_table(C0_REF_genera))
C0_REF_no_alpha_diff_abund <- merge(C0_REF_no_alpha_diff_abund, genera_C0_REF_df, by = "row.names")
C0_REF_no_alpha_diff_abund <- arrange(C0_REF_no_alpha_diff_abund, log2FoldChange)
dim(C0_REF_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_C0_REF_no_alpha_diff_abund <- filter(C0_REF_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_C0_REF_no_alpha_diff_abund <- filter(filter_C0_REF_no_alpha_diff_abund, log2FoldChange < 5)
filter_C0_REF_no_alpha_diff_abund <- filter(filter_C0_REF_no_alpha_diff_abund, log2FoldChange > -2.7)
filter_C0_REF_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_C0_REF_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()




# Differential abundance analysis: C1 vs REF
# Remove C0 treatment so that only C1, C2, and REF treatments remain
C1_C2_REF_physeq <- subset_samples(physeq, Compaction_Treatment != "C0")
C1_C2_REF_physeq

# Remove C2 so that only C1 and REF treatments remain
C1_REF_physeq <- subset_samples(C1_C2_REF_physeq, Compaction_Treatment != "C2")
C1_REF_physeq

# Determine the total counts and relative abundance for features across all C1 and REF samples
C1_REF_total_counts <- taxa_sums(C1_REF_physeq)
C1_REF_relative_abundance <- calculate_relative_abundance(C1_REF_total_counts)

# Remove low abundant features
C1_REF_abundant <- C1_REF_relative_abundance > 0.0005
C1_REF_abundant_taxa <- prune_taxa(C1_REF_abundant, C1_REF_physeq)
C1_REF_abundant_taxa

# Set taxonomic level for analysis (Genus)
C1_REF_genera <- tax_glom(C1_REF_abundant_taxa, taxrank = "Genus")
C1_REF_genera

# Create a DESeq2 object for C1 and REF
deseq_C1_REF <- phyloseq_to_deseq2(C1_REF_genera, ~ Compaction_Treatment)
geo_means_C1_REF <- apply(counts(deseq_C1_REF), 1, calculate_gm_mean)
deseq_C1_REF <- estimateSizeFactors(deseq_C1_REF, geoMeans = geo_means_C1_REF)
deseq_C1_REF <- DESeq(deseq_C1_REF, fitType = "local")

# Look at DESeq2 results for C1 vs REF
C1_REF_diff_abund <- results(deseq_C1_REF)
C1_REF_no_alpha_diff_abund <- as.data.frame(C1_REF_diff_abund)
dim(C1_REF_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C1_REF <- as.data.frame(C1_REF_diff_abund)
sig_C1_REF <- filter(sig_C1_REF, padj < alpha)
dim(sig_C1_REF)

# sig_C1_REF does not show any significant differences differential abundance between C1 and REF
# Output: [1] 0 6
# Use C1_REF_no_alpha_diff_abund for the next steps instead

# Merge C0 vs REF DESeq2 results with table of taxonomic info (C0_REF_genera)
genera_C1_REF_df <- as.data.frame(tax_table(C1_REF_genera))
C1_REF_no_alpha_diff_abund <- merge(C1_REF_no_alpha_diff_abund, genera_C1_REF_df, by = "row.names")
C1_REF_no_alpha_diff_abund <- arrange(C1_REF_no_alpha_diff_abund, log2FoldChange)
dim(C1_REF_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_C1_REF_no_alpha_diff_abund <- filter(C1_REF_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_C1_REF_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_C1_REF_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()




# Differential abundance analysis: C2 vs REF
# Remove C1 so that only C2 and REF treatments remain
C2_REF_physeq <- subset_samples(C1_C2_REF_physeq, Compaction_Treatment != "C1")
C2_REF_physeq

# Determine the total counts and relative abundance for features across all C2 and REF samples
C2_REF_total_counts <- taxa_sums(C2_REF_physeq)
C2_REF_relative_abundance <- calculate_relative_abundance(C2_REF_total_counts)

# Remove low abundant features
C2_REF_abundant <- C2_REF_relative_abundance > 0.0005
C2_REF_abundant_taxa <- prune_taxa(C2_REF_abundant, C2_REF_physeq)
C2_REF_abundant_taxa

# Set taxonomic level for analysis (Genus)
C2_REF_genera <- tax_glom(C2_REF_abundant_taxa, taxrank = "Genus")
C2_REF_genera

# Create a DESeq2 object for C2 and REF
deseq_C2_REF <- phyloseq_to_deseq2(C2_REF_genera, ~ Compaction_Treatment)
geo_means_C2_REF <- apply(counts(deseq_C2_REF), 1, calculate_gm_mean)
deseq_C2_REF <- estimateSizeFactors(deseq_C2_REF, geoMeans = geo_means_C2_REF)
deseq_C2_REF <- DESeq(deseq_C2_REF, fitType = "local")

# Look at DESeq2 results for C2 vs REF
C2_REF_diff_abund <- results(deseq_C2_REF)
C2_REF_no_alpha_diff_abund <- as.data.frame(C2_REF_diff_abund)
dim(C2_REF_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C2_REF <- as.data.frame(C2_REF_diff_abund)
sig_C2_REF <- filter(sig_C2_REF, padj < alpha)
dim(sig_C2_REF)

# sig_C2_REF does not show any significant differences differential abundance between C2 and REF
# Output: [1] 0 6
# Use C2_REF_no_alpha_diff_abund for the next steps instead

# Merge C2 vs REF DESeq2 results with table of taxonomic info (C2_REF_genera)
genera_C2_REF_df <- as.data.frame(tax_table(C2_REF_genera))
C2_REF_no_alpha_diff_abund <- merge(C2_REF_no_alpha_diff_abund, genera_C2_REF_df, by = "row.names")
C2_REF_no_alpha_diff_abund <- arrange(C2_REF_no_alpha_diff_abund, log2FoldChange)
dim(C2_REF_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_C2_REF_no_alpha_diff_abund <- filter(C2_REF_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_C2_REF_no_alpha_diff_abund <- filter(filter_C2_REF_no_alpha_diff_abund, log2FoldChange < 1)
filter_C2_REF_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_C2_REF_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()




# Differential abundance analysis: C0 vs C1
# Remove C2 treatment so that only C0, C1, and REF treatments remain
C0_C1_REF_physeq <- subset_samples(physeq, Compaction_Treatment != "C2")
C0_C1_REF_physeq

# Remove REF so that only C0 and C1 treatments remain
C0_C1_physeq <- subset_samples(C0_C1_REF_physeq, Compaction_Treatment != "REF")
C0_C1_physeq

# Determine the total counts and relative abundance for features across all C0 and C1 samples
C0_C1_total_counts <- taxa_sums(C0_C1_physeq)
C0_C1_relative_abundance <- calculate_relative_abundance(C0_C1_total_counts)

# Remove low abundant features
C0_C1_abundant <- C0_C1_relative_abundance > 0.0005
C0_C1_abundant_taxa <- prune_taxa(C0_C1_abundant, C0_C1_physeq)
C0_C1_abundant_taxa

# Set taxonomic level for analysis (Genus)
C0_C1_genera <- tax_glom(C0_C1_abundant_taxa, taxrank = "Genus")
C0_C1_genera

# Create a DESeq2 object for C0 and C1
deseq_C0_C1 <- phyloseq_to_deseq2(C0_C1_genera, ~ Compaction_Treatment)
geo_means_C0_C1 <- apply(counts(deseq_C0_C1), 1, calculate_gm_mean)
deseq_C0_C1 <- estimateSizeFactors(deseq_C0_C1, geoMeans = geo_means_C0_C1)
deseq_C0_C1 <- DESeq(deseq_C0_C1, fitType = "local")

# Look at DESeq2 results for C0 and C1
C0_C1_diff_abund <- results(deseq_C0_C1)
C0_C1_no_alpha_diff_abund <- as.data.frame(C0_C1_diff_abund)
dim(C0_C1_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C0_C1 <- as.data.frame(C0_C1_diff_abund)
sig_C0_C1 <- filter(sig_C0_C1, padj < alpha)
dim(sig_C0_C1)

# sig_C0_C1 shows that there is only 1 significant difference in differential abundance between C0 and C1
# Output: [1] 1 6
# Use sig_C0_C1 for the next steps instead

# Merge C0 vs C1 DESeq2 results with table of taxonomic info (C0_REF_genera)
genera_C0_C1_df <- as.data.frame(tax_table(C0_C1_genera))
sig_C0_C1 <- merge(sig_C0_C1, genera_C0_C1_df, by = "row.names")
sig_C0_C1 <- arrange(sig_C0_C1, log2FoldChange)
dim(sig_C0_C1)


# Visualize differential abundance results in a barplot
# g_ is not assigned, so use order instead
ggplot(sig_C0_C1, aes(x = log2FoldChange, y = Order)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant order",
       x = expression(log[2]~fold~change),
       y = "Order") +
  theme_bw()

# Bar plot for differential abundance for C0 and C1 doesn't look great, so do relative abundance
# Calculate relative abundance of C0 and C1
C0_C1_RA <- transform_sample_counts(C0_C1_physeq, calculate_relative_abundance)
C0_C1_counts <- taxa_sums(C0_C1_physeq)
rela_abund_C0_C1 <- calculate_relative_abundance(C0_C1_counts)
abund_C0_C1 <- rela_abund_C0_C1 > 0.0005
abund_C0_C1_RA_taxa <- prune_taxa(abund_C0_C1, C0_C1_RA)

# Subset order of interest
abund_C0_C1_RA_order <-  tax_glom(abund_C0_C1_RA_taxa, taxrank = "Order")
Solirubrobacterales <- subset_taxa(abund_C0_C1_RA_order, Order == "o__Solirubrobacterales")
otu_table(Solirubrobacterales)

# Transform data frame from wide to long format
Solirubrobacterales_long <- psmelt(Solirubrobacterales)
Solirubrobacterales_long

# Visualize relative abundance of Solirubrobacterales between C0 and C1 in a boxplot 
ggplot(Solirubrobacterales_long, aes(x = Compaction_Treatment, y = Abundance)) +
  geom_boxplot() +
  theme_bw()+
  labs(title = "Relative abundance of Solirubrobacterales",
       x     = "Compaction treatment",
       y     = "Relative abundance")





# Differential abundance analysis: C0 vs C2
# Remove C1 treatment so that only C0, C2, and REF treatments remain
C0_C2_REF_physeq <- subset_samples(physeq, Compaction_Treatment != "C1")
C0_C2_REF_physeq

# Remove  REF so that only C0 and C2 treatments remain
C0_C2_physeq <- subset_samples(C0_C2_REF_physeq, Compaction_Treatment != "REF")
C0_C2_physeq

# Determine the total counts and relative abundance for features across all C0 and C2 samples
C0_C2_total_counts <- taxa_sums(C0_C2_physeq)
C0_C2_relative_abundance <- calculate_relative_abundance(C0_C2_total_counts)

# Remove low abundant features
C0_C2_abundant <- C0_C2_relative_abundance > 0.0005
C0_C2_abundant_taxa <- prune_taxa(C0_C2_abundant, C0_C2_physeq)
C0_C2_abundant_taxa

# Set taxonomic level for analysis (Genus)
C0_C2_genera <- tax_glom(C0_C2_abundant_taxa, taxrank = "Genus")
C0_C2_genera

# Create a DESeq2 object for C0 and C2
deseq_C0_C2 <- phyloseq_to_deseq2(C0_C2_genera, ~ Compaction_Treatment)
geo_means_C0_C2 <- apply(counts(deseq_C0_C2), 1, calculate_gm_mean)
deseq_C0_C2 <- estimateSizeFactors(deseq_C0_C2, geoMeans = geo_means_C0_C2)
deseq_C0_C2 <- DESeq(deseq_C0_C2, fitType = "local")

# Look at DESeq2 results for C0 and C2
C0_C2_diff_abund <- results(deseq_C0_C2)
C0_C2_no_alpha_diff_abund <- as.data.frame(C0_C2_diff_abund)
dim(C0_C2_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C0_C2 <- as.data.frame(C0_C2_diff_abund)
sig_C0_C2 <- filter(sig_C0_C2, padj < alpha)
dim(sig_C0_C2)

# sig_C0_C2 does not show any significant difference in differential abundance between C0 and C2
# Output: [1] 0 6
# Use C0_C2_no_alpha_diff_abund for the next steps instead

# Merge C0 vs C2 DESeq2 results with table of taxonomic info (C0_C2_genera)
genera_C0_C2_df <- as.data.frame(tax_table(C0_C2_genera))
C0_C2_no_alpha_diff_abund <- merge(C0_C2_no_alpha_diff_abund, genera_C0_C2_df, by = "row.names")
C0_C2_no_alpha_diff_abund <- arrange(C0_C2_no_alpha_diff_abund, log2FoldChange)
dim(C0_C2_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_C0_C2_no_alpha_diff_abund <- filter(C0_C2_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_C0_C2_no_alpha_diff_abund


# Visualize differential abundance results in a barplot
ggplot(filter_C0_C2_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()





# Differential abundance analysis: C1 vs C2
# Remove C0 treatment so that only C1, C2, and REF treatments remain
C1_C2_REF_physeq <- subset_samples(physeq, Compaction_Treatment != "C0")
C1_C2_REF_physeq

# Remove REF so that only C1 and C2 treatments remain
C1_C2_physeq <- subset_samples(C1_C2_REF_physeq, Compaction_Treatment != "REF")
C1_C2_physeq

# Determine the total counts and relative abundance for features across all C1 and C2 samples
C1_C2_total_counts <- taxa_sums(C1_C2_physeq)
C1_C2_relative_abundance <- calculate_relative_abundance(C1_C2_total_counts)

# Remove low abundant features
C1_C2_abundant <- C1_C2_relative_abundance > 0.0005
C1_C2_abundant_taxa <- prune_taxa(C1_C2_abundant, C1_C2_physeq)
C1_C2_abundant_taxa

# Set taxonomic level for analysis (Genus)
C1_C2_genera <- tax_glom(C1_C2_abundant_taxa, taxrank = "Genus")
C1_C2_genera

# Create a DESeq2 object for C1 and C2
deseq_C1_C2 <- phyloseq_to_deseq2(C1_C2_genera, ~ Compaction_Treatment)
geo_means_C1_C2 <- apply(counts(deseq_C1_C2), 1, calculate_gm_mean)
deseq_C1_C2 <- estimateSizeFactors(deseq_C1_C2, geoMeans = geo_means_C1_C2)
deseq_C1_C2 <- DESeq(deseq_C1_C2, fitType = "local")

# Look at DESeq2 results for C1 vs C2
C1_C2_diff_abund <- results(deseq_C1_C2)
C1_C2_no_alpha_diff_abund <- as.data.frame(C1_C2_diff_abund)
dim(C1_C2_no_alpha_diff_abund)

# Define the alpha level to determine significant changes
# Convert results from DESeq2 to a data frame for further processing
alpha <- 0.05
sig_C1_C2 <- as.data.frame(C1_C2_diff_abund)
sig_C1_C2 <- filter(sig_C1_C2, padj < alpha)
dim(sig_C1_C2)

# sig_C1_C2 does not show any significant differences differential abundance between C1 and c2
# Output: [1] 0 6
# Use C1_C2_no_alpha_diff_abund for the next steps instead

# Merge C1 vs C2 DESeq2 results with table of taxonomic info (C1_C2_genera)
genera_C1_C2_df <- as.data.frame(tax_table(C1_C2_genera))
C1_C2_no_alpha_diff_abund <- merge(C1_C2_no_alpha_diff_abund, genera_C1_C2_df, by = "row.names")
C1_C2_no_alpha_diff_abund <- arrange(C1_C2_no_alpha_diff_abund, log2FoldChange)
dim(C1_C2_no_alpha_diff_abund)

# Remove samples with log2foldchange = 0.00000000
filter_C1_C2_no_alpha_diff_abund <- filter(C1_C2_no_alpha_diff_abund, log2FoldChange != 0.00000000)
filter_C1_C2_no_alpha_diff_abund

# Visualize differential abundance results in a barplot
ggplot(filter_C1_C2_no_alpha_diff_abund, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()





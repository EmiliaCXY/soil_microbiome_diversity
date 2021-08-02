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

# IMPORTANT!!! Convert tree from multichotomous to dichotomous tree
tree <- multi2di(tree)

# Combine information into single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

# Rename Taxa group name
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")

# Filter by seq depth of 3038, our rarefaction depth
physeq_rar <- rarefy_even_depth(physeq, sample.size = 3038)

# Calculate relative abundance
physeq_RA <- transform_sample_counts(physeq, calculate_relative_abundance)

# Keep only abundant ASVs
# First determine total count of a feature in experiment
total_counts <- taxa_sums(physeq)

# Calculate relative abundance of each feature within experiment
relative_abundance <- calculate_relative_abundance(total_counts)

# Determine which ASVs are more abundant than 0.5% or ratio 0.005
abundant <- relative_abundance > 0.005 
abundant_taxa <- prune_taxa(abundant, physeq)

# select samples of interest
moisture4_7 <- subset_samples(physeq, Moisture_level != "Median")
moisture4_5 <- subset_samples(physeq, Moisture_level != "High")
moisture5_7 <- subset_samples(physeq, Moisture_level != "Low")

# choose the pair we want to compare !!!
moisture_for_comparison <- moisture4_7

# Remove low abundant features
moisture_counts <- taxa_sums(moisture_for_comparison)
relative_abundance_m4_8 <- calculate_relative_abundance(moisture_counts)
abundant_m <- relative_abundance_m4_8 > 0.0001
abundant_m_taxa <- prune_taxa(abundant_m, moisture_for_comparison)

# Set taxonomic level for analysis
abundant_m_taxa_genera <- tax_glom(abundant_m_taxa, taxrank = "Genus")

# Use DESEQ2 to calculate differential abundance
deseq_m <- phyloseq_to_deseq2(abundant_m_taxa_genera, ~ Moisture_level)
geo_means <- apply(counts(deseq_m), 1, calculate_gm_mean)
deseq_m <- estimateSizeFactors(deseq_m, geoMeans = geo_means)
deseq_m <- DESeq(deseq_m, fitType = "local")

# Load results
m_diff_abund <- results(deseq_m)

# Filter for significant results
alpha <- 0.05
significant_m <- as.data.frame(m_diff_abund)
significant_m <- filter(significant_m, padj < alpha)

# Tabulate results
genera_df <- as.data.frame(tax_table(abundant_m_taxa_genera))
significant_m <- merge(significant_m, genera_df, by = "row.names")
significant_m <- arrange(significant_m, log2FoldChange)

# Plotting bar plot of genus that are differentially abundant
no_dup_sig_m <- significant_m[!duplicated(significant_m$Genus), ]
no_dup_sig_m <- mutate(no_dup_sig_m,
                       Genus = factor(Genus, levels = Genus))

ggplot(no_dup_sig_m, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()


#Calculate relative abundance for whole physeq genera
moisture_RA <- transform_sample_counts(moisture_for_comparison, calculate_relative_abundance)
moisture_counts <- taxa_sums(moisture_for_comparison)

# Caculate relative abundacne within just samples of interests
relative_abundance_m <- calculate_relative_abundance(moisture_counts)

# Filter for abundant variants
abundant_m <- relative_abundance_m > 0.0005
abundant_m_taxa <- prune_taxa(abundant_m, moisture4_7_RA)
abundant_m_taxa_genera <- tax_glom(abundant_m_taxa, taxrank = "Genus")
mycobacterium <- subset_taxa(abundant_m_taxa_genera, Genus == "g__Mycobacterium")
mycobacterium_long <- psmelt(mycobacterium)

# Plotting
ggplot(mycobacterium_long, aes(x = Moisture_level, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Differential abundant Mycobacterium",
       x = "moisture level",
       y = "relative abundance") +
  theme_bw()
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
set.seed(711)


# Export .biom file and tree from QIIME2 with metadata file
biom_file <- import_biom("table-with-taxonomy.biom")
metadata  <- import_qiime_sample_data("MICB_421_Soil_Metadata.tsv")
tree      <- read_tree_greengenes("tree.nwk")


# Convert from multichotomous to dichotmous tree
tree <- multi2di(tree)


# Combine all information into a single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class","Order", "Family",
                                 "Genus", "Species")
#Calculate relative abundance
physeq_RA <- transform_sample_counts(physeq, calculate_relative_abundance)


#Remove low-abundance features
total_counts <- taxa_sums(physeq)
relative_abundance <- calculate_relative_abundance(total_counts)
abundant <- relative_abundance > 0.0005 
abundant_taxa <- prune_taxa(abundant, physeq)


#Set taxonomic level for analysis
abundant_soil_genera <- tax_glom(abundant_taxa, taxrank = "Genus")
abundant_soil_genera

#Creating a DeSeq2 object 
deseq_soilgen <- phyloseq_to_deseq2(abundant_soil_genera, ~ Mean_Annual_Precipitation_mm)
geo_means <- apply(counts(deseq_soilgen), 1, calculate_gm_mean)
deseq_soilgen <- estimateSizeFactors(deseq_soilgen, geoMeans = geo_means)
deseq_soilgen <- DESeq(deseq_soilgen, fitType = "local")

soilgen_diff_abund <- results(deseq_soilgen)



#Set alpha to find taxa that are significantly different between precipitation levels
alpha <- 0.05
significant_soilgen <- as.data.frame(soilgen_diff_abund)
significant_soilgen <- filter(significant_soilgen, padj < alpha)


#Sort taxa as a data table
genera_df <- as.data.frame(tax_table(abundant_soil_genera))
significant_soilgen  <- merge(significant_soilgen, genera_df, by = "row.names")
significant_soilgen  <- arrange(significant_soilgen, log2FoldChange)
dim(significant_soilgen)


#Visualizing differential abundance data
significant_soilgen <- significant_soilgen[!duplicated(significant_soilgen$Genus), ]
significant_soilgen <- mutate(significant_soilgen,
                              Genus = factor(Genus, levels = Genus))

ggplot(significant_soilgen, aes(x = log2FoldChange, y = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Differential abundant genera",
       x = expression(log[2]~fold~change),
       y = "Genus") +
  theme_bw()

#Calculate relative abundance for whole physeq genera
soil_RA <- transform_sample_counts(physeq, calculate_relative_abundance)
soil_counts <- taxa_sums(physeq)
relative_abundance_soil <- calculate_relative_abundance(soil_counts)
abundant_soil <- relative_abundance_soil > 0.0005
abundant_soil_RA_taxa <- prune_taxa(abundant_soil, soil_RA)

#Subset to genus of interest (mycobacterium)
abundant_soil_RA_taxa <-  tax_glom(abundant_soil_RA_taxa, taxrank = "Genus")
mycobacterium <- subset_taxa(abundant_soil_RA_taxa, Genus == "g__Mycobacterium")
otu_table(mycobacterium)

mycobacterium_long <- psmelt(mycobacterium)
mycobacterium_long
ggplot(mycobacterium_long, aes(x = Mean_Annual_Precipitation_mm, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Mycobacterium",
       x     = "Mean Annual Precipitation (mm)",
       y     = "Relative abundance")+
  theme_bw()

bradyhizobium <- subset_taxa(abundant_soil_RA_taxa, Genus == "g__Bradyrhizobium")
bradyhizobium_long <- psmelt(bradyhizobium)
ggplot(bradyhizobium_long, aes(x = Mean_Annual_Precipitation_mm, y = Abundance)) +
  geom_boxplot() +
  labs(title = "Bradyrhizobium",
       x     = "Mean Annual Precipitation (mm)",
       y     = "Relative abundance")+
  theme_bw()

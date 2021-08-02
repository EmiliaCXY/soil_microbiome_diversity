# Load CRAN packages
library(tidyverse)
library(vegan)
library(ape)

# Load Bioconductor packages
library(phyloseq)

# Import qiime files
biom_file <- import_biom("table-with-taxonomy.biom")
metadata <- import_qiime_sample_data("MICB_421_Soil_Metadata.tsv")
tree <- read_tree_greengenes("tree.nwk")

# IMPORTANT!!! Convert tree from multichotomous to dichotomous tree
tree <- multi2di(tree)

# Combine information into single phyloseq object
physeq <- merge_phyloseq(biom_file, metadata, tree)

physeq_rar <- rarefy_even_depth(physeq, sample.size = 3038)

#================================================================================
# Plotting weighted unifrac PCA plots
#================================================================================
ord <- ordinate(physeq_rar, method = "PCoA", distance = "wunifrac")

# plotting weighted unifrac PCA on annual precipitation
plot_ordination(physeq_rar,
                ord,
                color = "Mean_Annual_Precipitation_mm") +
  # Define title of plot
  labs(title = "PCoA (weighted UniFrac)") +
  theme_bw()+
  # Postion a ellipses around each group of data, `type` determines how centre
  # is determined, `size` the thickness of the line
  stat_ellipse(type = "norm", size = 0.5) + 
  geom_point(size = 2) +
  scale_colour_manual(values = c("#00cc77", "#ff6030"),
                      labels = c("146-193", "300")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Mean annual\nprecipitation (mm)"))

# plotting weighted unifrac PCA on moisture bins
plot_ordination(physeq_rar,
                ord,
                color = "Moisture_bin") +
  # Define title of plot
  labs(title = "PCoA (weighted UniFrac)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(labels = c("40-50", "50-60", "60-70", "70-80", "80-90","90-100"),
                      values = c("#d73027","#fdae61","#fee090","#7dcbff","#69a0ff","#0000ff")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Moisture content"))

# plotting weighted unifrac PCA on OM
plot_ordination(physeq_rar,
                ord,
                color = "LTSP.Treatment") +
  # Define title of plot
  labs(title = "PCoA (weighted UniFrac)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(labels = c("OM1", "OM2", "Reference"),
                      values = c("#ff5ca2","#69a0ff","#898b8a")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Organic matter removal"))

# plotting weighted unifrac PCA on soil compaction
plot_ordination(physeq_rar,
                ord,
                color = "Compaction.Treatment") +
  # Define title of plot
  labs(title = "PCoA (weighted UniFrac)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(labels = c("C0", "C1", "C2", "Reference"),
                      values = c("#db6400","#61b15a","#892cdc", "#898b8a")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Soil compaction level"))


#================================================================================
# Plotting Unweighted Unifrac distance
#================================================================================ 
ord <- ordinate(physeq_rar, method = "PCoA", distance = "uunifrac")
plot_ordination(physeq_rar,
                ord,
                color = "LTSP.Treatment") +
  # Define title of plot
  labs(title = "PCoA (unweighted UniFrac)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(labels = c("OM1", "OM2", "Reference"),
                      values = c("#ff5ca2","#69a0ff","#898b8a")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Organic matter removal"))

#================================================================================
# Plotting Jaccard index
#================================================================================ 
ord <- ordinate(physeq_rar, method = "PCoA", distance = "jaccard")
plot_ordination(physeq_rar,
                ord,
                color = "Mean_Annual_Precipitation_mm") +
  # Define title of plot
  labs(title = "PCoA (Jaccard)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(values = c("#00cc77", "#ff6030"),
                      labels = c("146-193", "300")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Mean annual\nprecipitation (mm)"))

#================================================================================
# Plotting Bray Curtis index
#================================================================================ 
ord <- ordinate(physeq_rar, method = "PCoA", distance = "bray")
plot_ordination(physeq_rar,
                ord,
                color = "Mean_Annual_Precipitation_mm") +
  # Define title of plot
  labs(title = "PCoA (Bray)") +
  theme_bw()+
  geom_point(size = 2) +
  scale_colour_manual(values = c("#00cc77", "#ff6030"),
                      labels = c("146-193", "300")) +
  # Themes change the basic look of your plot
  # Rename the title of your legend
  guides(color = guide_legend("Mean annual\nprecipitation (mm)"))
#================================================================================
# Plotting Shannon diversity
#================================================================================

# inputing table
shannon_precipitation <- read.table(file='shannon_precipitation.tsv', sep = '\t', header = TRUE)

# ploting for precipitation
ggplot(shannon_precipitation, aes(x = Mean_Annual_Precipitation_mm, y = shannon_entropy, fill=Mean_Annual_Precipitation_mm))+
  geom_boxplot() +
  labs(title = "Shannon diveristy analysis",
       x     = "Mean annual precipitation (mm)",
       y     = "Shannon diveristy index",
       fill = "Mean annual\nprecipitation (mm)") +
  scale_fill_manual(labels = c("146-193", "300"),
                    values = c("#00cc77", "#ff6030"))+
  theme_bw() +
  guides(color = guide_legend("Mean annual precipitation level (mm)"))

# shannon for moisture
shannon_moisture <- read.table(file='shannon_soil_moisture.tsv', sep = '\t', header = TRUE)

# performing spearman correlation test on moisture
res <-cor.test(shannon_moisture$Moisture_bin, shannon_moisture$shannon_entropy,  method = "spearman")

shannon_precipitation$Moisture_bin <- as.character(shannon_precipitation$Moisture_bin)
ggplot(shannon_precipitation, aes(x = Moisture_bin, y = shannon_entropy, fill=Moisture_bin))+
  geom_boxplot() +
  labs(title = "Shannon diveristy analysis",
       x     = "Moisture content",
       y     = "Shannon diveristy index",
       fill = "Moisture content") +
  theme_bw() +
  scale_fill_manual(labels = c("40-50", "50-60", "60-70", "70-80", "80-90","90-100"),
                      values = c("#d73027","#fdae61","#fee090","#7dcbff","#69a0ff","#0000ff"))

# shannon for OM
ggplot(shannon_precipitation, aes(x = OM_Removal, y = shannon_entropy, fill=OM_Removal))+
  geom_boxplot() +
  labs(title = "Shannon diveristy analysis",
       x     = "Organic matter removal",
       y     = "Shannon diveristy index",
       fill = "Organic matter removal") +
  theme_bw() +
  scale_fill_manual(labels = c("OM1", "OM2", "Reference"),
                                  values = c("#ff5ca2","#69a0ff","#898b8a"))

# shannon for compaction
ggplot(shannon_precipitation, aes(x = Compaction_Treatment, y = shannon_entropy, fill=Compaction_Treatment))+
  geom_boxplot() +
  labs(title = "Shannon diveristy analysis",
       x     = "Soil compaction",
       y     = "Shannon diveristy index",
       fill = "Soil compaction") +
  theme_bw() +
  scale_fill_manual(labels = c("C0", "C1", "C2", "Reference"),
                      values = c("#db6400","#61b15a","#892cdc", "#898b8a"))



#================================================================================
# Plotting Faith diversity
#================================================================================

# inputting data
faith <- read.table(file='faith.tsv', sep = '\t', header = TRUE)

# Plot faith for precipitation
ggplot(faith, aes(x = Mean_Annual_Precipitation_mm, y = faith_pd, fill=Mean_Annual_Precipitation_mm))+
  geom_boxplot() +
  labs(title = "Faith diveristy analysis",
       x     = "Mean annual precipitation (mm)",
       y     = "Faith index",
       fill = "Mean annual\nprecipitation (mm)") +
  scale_fill_manual(labels = c("146-193", "300"),
                    values = c("#00cc77", "#ff6030"))+
  theme_bw() +
  guides(color = guide_legend("Mean annual precipitation level (mm)"))

# Faith for moisture

# performing spearman correlation test on moisture
res2 <-cor.test(faith$Moisture.Content, faith$faith_pd,  method = "spearman")

faith$Moisture_bin <- faith$Moisture_bin <- faith$Moisture_bin <- as.character(faith$Moisture_bin)
ggplot(faith, aes(x = Moisture_bin, y = faith_pd, fill=Moisture_bin))+
  geom_boxplot() +
  labs(title = "Faith diveristy analysis",
       x     = "Moisture content",
       y     = "Faith index",
       fill = "Moisture content") +
  theme_bw() +
  scale_fill_manual(labels = c("40-50", "50-60", "60-70", "70-80", "80-90","90-100"),
                    values = c("#d73027","#fdae61","#fee090","#7dcbff","#69a0ff","#0000ff"))

# Faith for OM
ggplot(faith, aes(x = LTSP.Treatment, y = faith_pd, fill=LTSP.Treatment))+
  geom_boxplot() +
  labs(title = "Faith diveristy analysis",
       x     = "Organic matter removal",
       y     = "Faith index",
       fill = "Organic matter removal") +
  theme_bw() +
  scale_fill_manual(labels = c("OM1", "OM2", "Reference"),
                    values = c("#ff5ca2","#69a0ff","#898b8a"))

# Faith for compaction
ggplot(faith, aes(x = Compaction.Treatment, y = faith_pd, fill=Compaction.Treatment))+
  geom_boxplot() +
  labs(title = "Faith diveristy analysis",
       x     = "Soil compaction",
       y     = "Faith diveristy index",
       fill = "Soil compaction") +
  theme_bw() +
  scale_fill_manual(labels = c("C0", "C1", "C2", "Reference"),
                    values = c("#db6400","#61b15a","#892cdc", "#898b8a"))




#================================================================================
# Plotting Pielou evenness
#================================================================================

# input table
evenness <- read.table(file='evenness.tsv', sep = '\t', header = TRUE)

# Plot evenness for precipitation
ggplot(evenness, aes(x = Mean_Annual_Precipitation_mm, y = pielou_evenness, fill=Mean_Annual_Precipitation_mm))+
  geom_boxplot() +
  labs(title = "Evenness analysis",
       x     = "Mean annual precipitation (mm)",
       y     = "Pieloue evenness",
       fill = "Mean annual\nprecipitation (mm)") +
  scale_fill_manual(labels = c("146-193", "300"),
                    values = c("#00cc77", "#ff6030"))+
  theme_bw() +
  guides(color = guide_legend("Mean annual precipitation level (mm)"))

# evenness for moisture
evenness$Moisture_bin <- as.character(evenness$Moisture_bin)
ggplot(evenness, aes(x = Moisture_bin, y = pielou_evenness, fill=Moisture_bin))+
  geom_boxplot() +
  labs(title = "Evenness analysis",
       x     = "Moisture content",
       y     = "Pieloue evenness",
       fill = "Moisture content") +
  theme_bw() +
  scale_fill_manual(labels = c("40-50", "50-60", "60-70", "70-80", "80-90","90-100"),
                    values = c("#d73027","#fdae61","#fee090","#7dcbff","#69a0ff","#0000ff"))

# evenness for OM
ggplot(evenness, aes(x = LTSP.Treatment, y = pielou_evenness, fill=LTSP.Treatment))+
  geom_boxplot() +
  labs(title = "Evenness analysis",
       x     = "Organic matter removal",
       y     = "Pieloue evenness",
       fill = "Organic matter removal") +
  theme_bw() +
  scale_fill_manual(labels = c("OM1", "OM2", "Reference"),
                    values = c("#ff5ca2","#69a0ff","#898b8a"))

# evenness for compaction
ggplot(evenness, aes(x = Compaction.Treatment, y = pielou_evenness, fill=Compaction.Treatment))+
  geom_boxplot() +
  labs(title = "Evenness analysis",
       x     = "Soil compaction",
       y     = "Pielou evenness",
       fill = "Soil compaction") +
  theme_bw() +
  scale_fill_manual(labels = c("C0", "C1", "C2", "Reference"),
                    values = c("#db6400","#61b15a","#892cdc", "#898b8a"))


#================================================================================
# Plotting PERMANOVA restuls on weighted unifrac
#================================================================================

# pairwise comparison on precipitation
# input table
pair_preci <- read.table(file='weighted_unifrac_precipitation_PERMANOVA.tsv', sep = '\t', header = TRUE)

pair_preci <- pair_preci %>% filter(Group1 == "146-193")
ggplot(pair_preci, aes(x = Group2, y = Distance))+
  geom_boxplot() +
  labs(title = "Weighted unifrac distance to 146-193",
       x     = "Mean annual precipitation (mm)",
       y     = "Distance") +
  theme_bw()

# pairwise comparison on moisture
# input table
pair_moisture <- read.table(file='moisture_weighted_unifrac_PERMANOVA.tsv', sep='\t', header= TRUE)

pair_preci$Group1 = as.character(pair_preci$Group1)
pair_preci$Group2 = as.character(pair_preci$Group2)
pair_preci <- pair_preci %>% filter(Group1 == "4")
ggplot(pair_preci, aes(x = Group2, y = Distance))+
  geom_boxplot() +
  labs(title = "Weighted unifrac distance to 40-50",
       x     = "Moisture content",
       y     = "Distance") +
  theme_bw()


# pairwise comparison on organic matter removal
# input table
pair_OM <- read.table(file='OM_pair_PERMANOVA.tsv', sep='\t', header= TRUE)
pair_OM$Group1 = as.character(pair_OM$Group1)
pair_OM$Group2 = as.character(pair_OM$Group2)
pair_OM <- pair_OM %>% filter(Group1 == "OM1")
ggplot(pair_OM, aes(x = Group2, y = Distance))+
  geom_boxplot() +
  labs(title = "Weighted unifrac distance to OM1",
       x     = "Organic Matter Removal",
       y     = "Distance") +
  theme_bw()

# pairwise comparison on soil compaction
# input table
pair_compaction <- read.table(file='compaction_pair_PERMANOVA.tsv', sep='\t', header= TRUE)
pair_compaction$Group1 = as.character(pair_compaction$Group1)
pair_compaction$Group2 = as.character(pair_compaction$Group2)
pair_compaction <- pair_compaction %>% filter(Group1 == "C1")
ggplot(pair_compaction, aes(x = Group2, y = Distance))+
  geom_boxplot() +
  labs(title = "Weighted unifrac distance to C1",
       x     = "Soil compaction",
       y     = "Distance") +
  theme_bw()
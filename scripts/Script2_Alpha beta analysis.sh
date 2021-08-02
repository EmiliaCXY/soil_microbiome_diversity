#======================================================================================
# Alpha Diversity Analysis: Shannon and Faith index analysis
#======================================================================================

# Visualizing alpha diversity results
qiime diversity alpha-group-significance \
   --i-alpha-diversity ../core-metrics-results/faith_pd_vector.qza \
   --m-metadata-file ../metadata_related/MICB_421_Soil_Metadata.tsv \
   --o-visualization ../core-metrics-results/faiths_pd_statistics.qzv


qiime diversity alpha-group-significance \
   --i-alpha-diversity ../core-metrics-results/faith_pd_vector.qza \
   --m-metadata-file ../metadata_related/MICB_421_Soil_Metadata.tsv \
   --o-visualization ../core-metrics-results/faiths_pd_statistics.qzv


#======================================================================================
# Beta Diversity Analysis: Weighted and Unweighted UniFrac
#======================================================================================

# Perform PERMANOVA to test whether mean annual precipitation is associated with significant 
# differences in unweighted and weighted UniFrac distances

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Mean Annual Precipitation (mm)' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-precipitation-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Mean Annual Precipitation (mm)' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-precipitation-significance.qzv

# Pairwise PERMANOVA for mean annual precipitation
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Mean Annual Precipitation (mm)' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-pairwise-precipitation-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Mean Annual Precipitation (mm)' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-pairwise-precipitation-significance.qzv \
  --p-pairwise

# Perform PERMANOVA to test whether compaction treatment is associated with significant 
# differences in unweighted and weighted UniFrac distances
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Compaction Treatment' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-compaction-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Compaction Treatment' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-compaction-significance.qzv

# Pairwise PERMANOVA for compaction treatment
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Compaction Treatment' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-pairwise-compaction-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Compaction Treatment' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-pairwise-compaction-significance.qzv \
  --p-pairwise

# Perform PERMANOVA to test whether moisture content (binned) is associated with significant
# differences in unweighted and weighted UniFrac distances
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Moisture bin' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-moisture-bin-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Moisture bin' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-moisture-bin-significance.qzv

# Pairwise PERMANOVA for moisture (binned)
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Moisture bin' \
  --o-visualization ./beta-diversity-results/unweighted_unifrac/unweighted-unifrac-pairwise-moisture-bin-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Moisture bin' \
  --o-visualization ./beta-diversity-results/weighted_unifrac/weighted-unifrac-pairwise-moisture-bin-significance.qzv \
  --p-pairwise


# Perform PERMANOVA to test whether organic matter (OM) removal is associated with significant
# differences in unweighted and weighted UniFrac distances
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'LTSP Treatment' \
  --o-visualization ./beta-diversity-results/unweighted-unifrac-OM-removal-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'LTSP Treatment' \
  --o-visualization ./beta-diversity-results/weighted-unifrac-OM-removal-significance.qzv

# Pairwise PERMANOVA for OM removal
qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'LTSP Treatment' \
  --o-visualization ./beta-diversity-results/unweighted_unifrac/unweighted-unifrac-pairwise-OM-removal-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'LTSP Treatment' \
  --o-visualization ./beta-diversity-results/weighted_unifrac/weighted-unifrac-pairwise-OM-removal-significance.qzv \
  --p-pairwise


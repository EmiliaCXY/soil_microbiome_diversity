# Activate QIIME2
conda activate qiime2-2020.8


# Tabulated metadata
qiime metadata tabulate \
--m-input-file MICB_421_Soil_Metadata.tsv \
--o-visualization soil_metadata.qzv


# Reimporting sequences
# Import and visualize demultiplexed data
qiime tools import \
 --type "SampleData[SequencesWithQuality]" \
 --input-format SingleEndFastqManifestPhred33V2 \
 --input-path ./'bc_o_ph_soil_manifest.txt'\
 --output-path ./bc_o_ph_demux_seqs.qza

# Visualize demultiplexed data
qiime demux summarize \
  --i-data ./bc_o_ph_demux_seqs.qza \
  --o-visualization ./bc_o_ph_demux_seqs.qzv

# Denoise data with DADA2 and trim reads to 340bp
qiime dada2 denoise-single \
--i-demultiplexed-seqs ./bc_o_ph_demux_seqs.qza  \
--p-trunc-len 340 \
--o-table ./bc_o_ph_dada2_table.qza \
--o-representative-sequences ./bc_o_ph_dada2_rep_set.qza \
--o-denoising-stats ./bc_o_ph_dada2_stats.qza

# Change tmp folder to ensure we have enough space on the server
Export $TMPDIR = ./data

# Building a phylogenetic tree for all the samples
qiime fragment-insertion sepp \
  --i-representative-sequences ./try_filter_manifest/ bc_o_ph_dada2_rep_set.qza \
  --i-reference-database sepp-refs-gg-13-8.qza \
  --o-tree ./bc_tree.qza \
  --o-placements ./bc_tree_placements.qza \
  --p-threads 1

# Generating diversity metrics 
qiime diversity core-metrics-phylogenetic \
  --i-table ./try_filter_manifest/bc_o_ph_data2_table.qza \
  --i-phylogeny ./bc_tree.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --p-sampling-depth 3038 \
  --output-dir ./core-metrics-results2

qiime diversity beta-group-significance \
  --i-distance-matrix ./core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ./metadata_related/MICB_421_Soil_Metadata.tsv \
  --m-metadata-column 'Moisture_level' \
  --o-visualization ./beta-diversity-results/weighted_unifrac/weighted-unifrac-pairwise-moisture-bin-significance.qzv \
  --p-pairwise


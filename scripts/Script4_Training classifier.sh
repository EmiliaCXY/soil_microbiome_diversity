# Building a reference database for classifier

mkdir training-feature-classifiers
cd training-feature-classifiers

# building reference OTUs
qiime tools import \   
--type 'FeatureData[Sequence]' \  
--input-path gg_13_8_otus/rep_set/99_otus.fasta \  
--output-path ref-otus.qza

# building reference taxonomy
qiime tools import \
 --type 'FeatureData[Taxonomy]' \
 --input-format HeaderlessTSVTaxonomyFormat \
 --input-path gg_13_8_otus/taxonomy/97_otu_taxonomy.txt \
 --output-path ref-taxonomy.qza

# extract V1-V3 from 16rRNA reads based on the primers used in the paper
qiime feature-classifier extract-reads \
 --i-sequences ref-otus.qza \
 --p-f-primer AGAGTTTGATCMTGGCTCAG \
 --p-r-primer GWATTACCGCGGCKGCTG\
 --p-trunc-len 340 \
 --o-reads ref-seqs.qza

# Training a naive Baysian model 
qiime feature-classifier fit-classifier-naive-bayes \
 --i-reference-reads ref-seqs.qza \
 --i-reference-taxonomy ref-taxonomy.qza \
 --o-classifier classifier.qza

# apply our data to the trained model
qiime feature-classifier classify-sklearn \
 --i-classifier classifier.qza \
 --i-reads ../try_filtering_manifest/bc_o_ph_rep-seqs.qza \
 --o-classification taxonomy.qza

# tabulate tax present in our samples
qiime metadata tabulate \
 --m-input-file taxonomy.qza \
 --o-visualization taxonomy.qzv

# produce taxa bar plot
qiime taxa barplot \
 --i-table ../try_filtering_manifest/bc_o_ph_dada2_table.qza \
 --i-taxonomy taxonomy.qza \
 --m-metadata-file ../metadata_related/MICB_421_Soil_Metadata.txt \
 --o-visualization taxa-bar-plots.qzv
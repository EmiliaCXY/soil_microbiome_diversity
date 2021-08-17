# Annual precipitation and soil moisture level strongly associate with the bacterial community structure in Interior Douglas-fir and Sub-Boreal Spruce ecozones in British Columbia

A bioinformatic pipeline using QIIME2 and R to analyze 16S rRNA sequencing data from forest soil metagenome to reveal bacterial community diversity and composition of two ecozones in British Columbia.

## Dataset
Sequencing data obtained from Wilhelm et al (2017) that examined the metagenome of forest soil microbial communities across 16 ecozones from sites in the Long-Term Soil Productivity (LTSP) in North America.

## Data preprocessing
The data were subset to 104 sites from the Interior Douglas-fir (IDF) and Sub-Boreal Spruce (SBS) ecozones in British Columbia (BC). 

Demultiplexing of the 16S rRNA data was acheived using QIIME2. DADA2 was used to denoise the data and classify amplicon sequence variants (ASV). 

## Phylogenetic tree construction
QIIME2 fragment insertion approach was used to generate a phylogenetic tree for the ASVs detected in the samples. `Greengenes (sepp-refs-gg-13-8.qza)` was used as a reference tree.

## Diversity anlaysis
Alpha and beta diversity analyses were performed in QIIME2 to correlate mean annual precipitation (MAP) and soil moisture levels with bacterial diversity.

## Taxanomic analysis
We trained a Naive Bayes classifier using the `q2-feature-classifier` plugin and the `Greengenes (release 13_8) 97% database` as training data.

## Differential abundance analysis
DESEQ2 package in R was adapted to analyze differental DNA sequence abundance.

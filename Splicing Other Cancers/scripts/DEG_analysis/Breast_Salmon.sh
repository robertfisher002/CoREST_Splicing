#!/bin/bash -l
#$ -pe omp 24
#$ -P epigen

module load salmon

# Set the path to the Salmon index 
SALMON_INDEX="/projectnb/epigen/Robert/ATRT/salmon/hg38_index"

# Set the input directory
INPUT_DIR="/projectnb/epigen/Robert/Breast/fastq"

# Set the output directory for Salmon quantification results
OUTPUT_DIR="/projectnb/epigen/Robert/Breast/salmon"

# Run Salmon to quantify expression
salmon quant -i $SALMON_INDEX \
             -l A \
             -r /projectnb/epigen/Robert/Breast/fastq/SRR13926994_cut.fastq.gz \
             -p 8 \
             -o /projectnb/epigen/Robert/Breast/salmon/SRR13926994_cut

# Run Salmon to quantify expression
salmon quant -i $SALMON_INDEX \
             -l A \
             -r /projectnb/epigen/Robert/Breast/fastq/SRR13926995_cut.fastq.gz \
             -p 8 \
             -o /projectnb/epigen/Robert/Breast/salmon/SRR13926995_cut

# Run Salmon to quantify expression
salmon quant -i $SALMON_INDEX \
             -l A \
             -r /projectnb/epigen/Robert/Breast/fastq/SRR13926998_cut.fastq.gz \
             -p 8 \
             -o /projectnb/epigen/Robert/Breast/salmon/SRR13926998_cut

# Run Salmon to quantify expression
salmon quant -i $SALMON_INDEX \
             -l A \
             -r /projectnb/epigen/Robert/Breast/fastq/SRR13926999_cut.fastq.gz \
             -p 8 \
             -o /projectnb/epigen/Robert/Breast/salmon/SRR13926999_cut
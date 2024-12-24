#!/bin/bash -l
#$ -pe omp 4
#$ -P epigen

module load star
module load python3
module load R
module load samtools
module load regtools

 # Run STAR
      STAR --genomeDir /projectnb/epigen/Robert/hg38 \
           --twopassMode Basic \
           --readFilesIn ../fastq/SRR13926994_cut.fastq.gz \
           --readFilesCommand zcat \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --outFileNamePrefix "DMSO_1_"

# Run STAR
      STAR --genomeDir /projectnb/epigen/Robert/hg38 \
           --twopassMode Basic \
           --readFilesIn ../fastq/SRR13926995_cut.fastq.gz \
           --readFilesCommand zcat \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --outFileNamePrefix "DMSO_2_"

for bamfile in `ls /projectnb/epigen/Robert/Breast/BAM/*.bam`; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -M 500000 -s unstranded $bamfile -o $bamfile.junc
    echo $bamfile.junc >> breast_juncfiles.txt
done

python /projectnb/epigen/Robert/melanoma_splice/leafcutter/clustering/leafcutter_cluster_regtools.py -j breast_juncfiles.txt –checkchrom flag -m 50 -o intron_cluster -l 500000

/projectnb/epigen/Robert/melanoma_splice/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -i 2 -g 2 /projectnb/epigen/Robert/Breast/BAM/intron_cluster/intron_cluster_perind_numers.counts.gz leaf_legend.txt  




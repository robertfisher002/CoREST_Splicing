#!/bin/bash -l
#$ -pe omp 8
#$ -P epigen

module load python3
module load R
module load samtools
module load regtools

for bamfile in `ls junctions/junc_WM793/*.bam`; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -M 500000 -s unstranded $bamfile -o $bamfile.junc
    echo $bamfile.junc >> WM793_juncfiles.txt
done

#cluster introns
python /projectnb/epigen/Robert/melanoma_splice/leafcutter/clustering/leafcutter_cluster_regtools.py -j WM793_juncfiles.txt –checkchrom flag -m 50 -o WM793_clusters/WM793_cluster -l 500000

#differential analysis DMSO v Corin
/projectnb/epigen/Robert/melanoma_splice/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -i 2 -g 2 /projectnb/epigen/Robert/ATRT/leafcutter/BAM/WM793_clusters/WM793_cluster_perind_numers.counts.gz leaf_legend_WM793.txt 
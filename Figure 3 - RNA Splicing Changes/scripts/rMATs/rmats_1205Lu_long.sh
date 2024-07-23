#!/bin/bash -l
#$ -pe omp 4
#$ -P epigen

module load miniconda/4.9.2
module load rmats/4.1.1

rmats.py --b1 txt_bam/1205Lu_DMSO.txt --b2 txt_bam/1205Lu_Corin.txt --gtf gtf/gencode.v42.annotation.gtf -t paired --readLength 150 --allow-clipping --nthread 4 --od rmats/1205_new --tmp rmats/tmp

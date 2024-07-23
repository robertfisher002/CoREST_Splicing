#!/bin/bash -l
#$ -pe omp 4
#$ -P epigen

module load python3
module load samtools

python src/rmats2sashimiplot/rmats2sashimiplot.py --b1 ../D1_chr.bam,../D2_chr.bam --b2 ../J1_chr.bam,../J2_chr.bam -c chr2:-:215379129:215382325:../gencode.v42.gff3 --l1 DMSO --l2 Corin --exon_s 1 --intron_s 20 -o FN1_output_long --color '#b2b2b2,#b2b2b2,#D2042D,#D2042D' --fig-width 5 --fig-height 6 --font-size 7

python src/rmats2sashimiplot/rmats2sashimiplot.py --b1 ../D1_chr.bam,../D2_chr.bam --b2 ../J1_chr.bam,../J2_chr.bam -c chr15:-:29718265:29720708:../gencode.v42.gff3 --l1 DMSO --l2 Corin --exon_s 1 --intron_s 20 -o TJP1_output_long --color '#b2b2b2,#b2b2b2,#D2042D,#D2042D' --fig-width 5 --fig-height 6 --font-size 7

python src/rmats2sashimiplot/rmats2sashimiplot.py --b1 ../D1_chr.bam,../D2_chr.bam --b2 ../J1_chr.bam,../J2_chr.bam -c chr2:+:191400381:191402718:../gencode.v42.gff3 --l1 DMSO --l2 Corin --exon_s 1 --intron_s 20 -o MYO1B_output_long --color '#b2b2b2,#b2b2b2,#D2042D,#D2042D' --fig-width 5 --fig-height 6 --font-size 7
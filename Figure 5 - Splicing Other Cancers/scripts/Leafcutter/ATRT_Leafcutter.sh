#!/bin/bash -l
#$ -pe omp 4
#$ -P epigen

module load star
module load python3
module load R
module load samtools
module load regtools

INPUT_DIR="/projectnb/epigen/ATRT_RNAseq_FASTQ"
OUTPUT_DIR="/projectnb/epigen/Robert/ATRT/leafcutter/BAM"  # Specify the path to your output directory

# Loop through subdirectories (sample replicates)
for sample_dir in ${INPUT_DIR}/*; do
  if [ -d "$sample_dir" ]; then
    # Get the sample replicate name (directory name)
    sample_name=$(basename "$sample_dir")
    output_bam="${OUTPUT_DIR}/${sample_name}_Aligned.sortedByCoord.out.bam"

    # Check if the output BAM file already exists
    if [ -f "$output_bam" ]; then
      echo "Output BAM for ${sample_name} already exists. Skipping STAR alignment."
      continue
    fi

    # Find read files in the sample replicate subdirectory
    read1_file="$(find "$sample_dir" -type f -name "*_1.fq.gz" | head -n 1)"
    read2_file="$(find "$sample_dir" -type f -name "*_2.fq.gz" | head -n 1)"

    if [ -n "$read1_file" ] && [ -n "$read2_file" ]; then
      # Run STAR to convert FASTQ to BAM
      STAR --genomeDir /projectnb/epigen/Robert/hg38 \
           --twopassMode Basic \
           --readFilesIn "$read1_file" "$read2_file" \
           --readFilesCommand zcat \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --outFileNamePrefix "${OUTPUT_DIR}/${sample_name}_"

      echo "Converted ${sample_name} to BAM."
    else
      echo "ERROR: Read files not found in ${sample_dir}" >&2
    fi
  fi
done

echo "Conversion complete."


for bamfile in `ls /projectnb/epigen/Robert/ATRT/leafcutter/BAM/BT37_BAM/*.bam`; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -M 500000 -s unstranded $bamfile -o $bamfile.junc
    echo $bamfile.junc >> BT37_juncfiles.txt
done


for bamfile in `ls /projectnb/epigen/Robert/ATRT/leafcutter/BAM/CHO6_BAM/*.bam`; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -M 500000 -s unstranded $bamfile -o $bamfile.junc
    echo $bamfile.junc >> CHO6_juncfiles.txt
done


python /projectnb/epigen/Robert/melanoma_splice/leafcutter/clustering/leafcutter_cluster_regtools.py -j CHO6_juncfiles.txt –checkchrom flag -m 50 -o CHO6_clusters/CHO6_cluster -l 500000

python /projectnb/epigen/Robert/melanoma_splice/leafcutter/clustering/leafcutter_cluster_regtools.py -j BT37_juncfiles.txt –checkchrom flag -m 50 -o BT37_clusters/BT37_cluster -l 500000

/projectnb/epigen/Robert/melanoma_splice/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -i 2 -g 2 /projectnb/epigen/Robert/ATRT/leafcutter/BAM/BT37_clusters/BT37_cluster_perind_numers.counts.gz leaf_legend_BT37.txt

/projectnb/epigen/Robert/melanoma_splice/leafcutter/scripts/leafcutter_ds.R --num_threads 4 -i 2 -g 2 /projectnb/epigen/Robert/ATRT/leafcutter/BAM/CHO6_clusters/CHO6_cluster_perind_numers.counts.gz leaf_legend_CHO6.txt 




#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Please provide the following arguments: REF_GENOME READS1 READS2 PREFIX OUTPUT_DIR"
    exit 1
fi

# Assign input arguments to variables
REF_GENOME=$1
READS1=$2
READS2=$3
PREFIX=$4
OUTPUT_DIR=$5

# Step 1: Align reads to reference genome
bwa mem -t 16 $REF_GENOME $READS1 $READS2 > $OUTPUT_DIR/$PREFIX"_aligned.sam"

# Step 2: Convert SAM to BAM
samtools view -bS $OUTPUT_DIR/$PREFIX"_aligned.sam" > $OUTPUT_DIR/$PREFIX"_aligned.bam"

# Step 3: Sort BAM file
samtools sort $OUTPUT_DIR/$PREFIX"_aligned.bam" -o $OUTPUT_DIR/$PREFIX"_aligned_sorted.bam"

# Step 4: Index sorted BAM file
samtools index $OUTPUT_DIR/$PREFIX"_aligned_sorted.bam"

# Step 5: Extract specific genomic region
samtools view -b $OUTPUT_DIR/$PREFIX"_aligned_sorted.bam" chrX:20000000-40000000 > $OUTPUT_DIR/$PREFIX"_part.bam"

# Step 6: Sort extracted region BAM file
samtools sort $OUTPUT_DIR/$PREFIX"_part.bam" -o $OUTPUT_DIR/$PREFIX"_part_sorted.bam"

# Step 7: Index sorted extracted region BAM file
samtools index $OUTPUT_DIR/$PREFIX"_part_sorted.bam"
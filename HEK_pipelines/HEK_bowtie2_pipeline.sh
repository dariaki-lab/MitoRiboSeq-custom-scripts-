#!/bin/bash

set -e

THREADS=4

INPUT_DIR="processed"
OUTPUT_DIR="processed"
LOG_DIR="logs"

BOWTIE2_INDEX="/Users/annkrg/Library/CloudStorage/OneDrive-KI.SE/Mac/Documents/ribo-profiling/refseq/Human_mtDNA_polyA2"

mkdir -p "$LOG_DIR"

echo "Starting Bowtie2..."

for file in "$INPUT_DIR"/*_cut.fastq.gz
do
    basename=$(basename "$file" _cut.fastq.gz)

    echo "Processing $basename..."

    bowtie2 \
        --local \
        -x "$BOWTIE2_INDEX" \
        -U "$file" \
        -p "$THREADS" \
        2> "$LOG_DIR/${basename}_bowtie2.log" \
    | samtools view -bS - \
    | tee "$OUTPUT_DIR/${basename}.bam" \
    | samtools view -b -F 4 - \
    > "$OUTPUT_DIR/${basename}_mapped.bam"

    samtools sort \
        "$OUTPUT_DIR/${basename}_mapped.bam" \
        -o "$OUTPUT_DIR/${basename}_mapped.sorted.bam"

    samtools index \
        "$OUTPUT_DIR/${basename}_mapped.sorted.bam"

    echo "$basename done"
done

echo "Bowtie2 finished!"

echo -e "Sample\tTotal\tUnaligned\tAligned_1x\tAligned_multi" > bowtie2_summary.tsv

for log in "$LOG_DIR"/*_bowtie2.log
do
    sample=$(basename "$log" _bowtie2.log)

    total=$(grep "were unpaired" "$log" | awk '{print $1}')
    unaligned=$(grep "aligned 0 times" "$log" | awk '{print $1}')
    aligned1=$(grep "aligned exactly 1 time" "$log" | awk '{print $1}')
    aligned_multi=$(grep "aligned >1 times" "$log" | awk '{print $1}')

    echo -e "$sample\t$total\t$unaligned\t$aligned1\t$aligned_multi" >> bowtie2_summary.tsv
done

echo "Bowtie2 summary written!"
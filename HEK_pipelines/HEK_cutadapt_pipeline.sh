#!/bin/bash

set -e

THREADS=4
ADAPTER="TGGAATTCTCGGGTGCCAAGG"
MINLEN=25
MAXLEN=45

INPUT_DIR="fastq"
OUTPUT_DIR="processed"
LOG_DIR="logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

echo "Starting cutadapt..."

for file in "$INPUT_DIR"/*.fastq.gz
do
    basename=$(basename "$file" .fastq.gz)

    echo "Processing $basename..."

    cutadapt \
        -a "$ADAPTER" \
        -o "$OUTPUT_DIR/${basename}_cut.fastq.gz" \
        "$file" \
        --max-n 0.2 \
        -m "$MINLEN" \
        -M "$MAXLEN" \
        --cores "$THREADS" \
        > "$LOG_DIR/${basename}_cutadapt.log"

    echo "$basename done"
done

echo "Cutadapt finished!"

echo -e "Sample\tTotal\tWithAdapters\tWritten" > cutadapt_summary.tsv

for log in "$LOG_DIR"/*_cutadapt.log
do
    sample=$(basename "$log" _cutadapt.log)

    total=$(grep "Total reads processed" "$log" | awk '{print $4}')
    adapted=$(grep "Reads with adapters" "$log" | awk '{print $4}')
    written=$(grep "Reads written (passing filters)" "$log" | awk '{print $5}')

    echo -e "$sample\t$total\t$adapted\t$written" >> cutadapt_summary.tsv
done

echo "Cutadapt summary written!"
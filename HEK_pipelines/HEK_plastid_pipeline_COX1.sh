#!/bin/bash

source /Users/annkrg/Applications/miniconda3/etc/profile.d/conda.sh
conda activate plastid

set -e

# ------------ PATHS ------------
DATA_DIR="processed"
PLASTID_DIR="plastid"

OFFSET_FILE="${PLASTID_DIR}/offsets/final_offsets_WT_mod.txt"
ROI_START="${PLASTID_DIR}/mtDNA_cds_start_COX1.txt"

ANNOT_GTF="${PLASTID_DIR}/Human_mtDNA_polyA2_mod2.gtf2"
ANNOT_BED="${PLASTID_DIR}/Human_mtDNA_polyA2_mod3_extended.bed"

# Create output folders
mkdir -p ${PLASTID_DIR}/phasing_COX1

echo "Starting plastid pipeline..."

# ------------ LOOP OVER BAM FILES ------------
for bam in ${DATA_DIR}/*_mapped.sorted.bam
do
    sample=$(basename ${bam} _mapped.sorted.bam)

    echo "Processing ${sample}..."

    # -------- PHASING --------
    phase_by_size \
        ${ROI_START} \
        ${PLASTID_DIR}/phasing_COX1/${sample} \
        --count_files ${bam} \
        --fiveprime_variable \
        --offset ${OFFSET_FILE} \
        --codon_buffer 5 \
        --min_length 25 --max_length 45

    echo "${sample} done"
done

echo "Plastid pipeline finished ✅"
``
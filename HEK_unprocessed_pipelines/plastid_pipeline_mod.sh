#!/bin/bash

conda activate plastid

set -e

# ------------ PATHS ------------
DATA_DIR="processed"
PLASTID_DIR="plastid"

OFFSET_FILE="${PLASTID_DIR}/offsets/final_offsets_WT_mod.txt"
ROI_START="${PLASTID_DIR}/mtDNA1_cds_start_rois_mod.txt"

ANNOT_GTF="${PLASTID_DIR}/mtDNA1_UTR.gtf2"
ANNOT_BED="${PLASTID_DIR}/mtDNA1_UTR_extended.bed"

# Create output folders
mkdir -p "${PLASTID_DIR}/phasing_mod"
mkdir -p "${PLASTID_DIR}/counts_mod"
mkdir -p "${PLASTID_DIR}/vectors_mod"

echo "Starting plastid pipeline..."

# ------------ LOOP OVER BAM FILES ------------
for bam in "${DATA_DIR}"/*_mapped.sorted.bam
do
    sample=$(basename "${bam}" _mapped.sorted.bam)

    echo "Processing ${sample}..."

    # -------- PHASING --------
    phase_by_size \
        "${ROI_START}" \
        "${PLASTID_DIR}/phasing_mod/${sample}" \
        --count_files "${bam}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --codon_buffer 5 \
        --min_length 28 \
        --max_length 37

    # -------- COUNTS --------
    counts_in_region \
        "${PLASTID_DIR}/counts_mod/Counts_${sample}.txt" \
        --count_files "${bam}" \
        --annotation_files "${ANNOT_GTF}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --min_length 28 \
        --max_length 37

    # -------- COUNT VECTORS --------
    get_count_vectors \
        --annotation_files "${ANNOT_BED}" \
        --annotation_format BED \
        --count_files "${bam}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --min_length 28 \
        --max_length 37 \
        "${PLASTID_DIR}/vectors_mod/${sample}"

    echo "${sample} done"
done

echo "Plastid pipeline finished!"
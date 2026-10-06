#!/bin/bash

source /Users/annkrg/Applications/miniconda3/etc/profile.d/conda.sh
conda activate plastid

# Continue even if one sample fails
set +e

# ------------ PATHS ------------

DATA_DIR="processed"
PLASTID_DIR="plastid"

OFFSET_FILE="${PLASTID_DIR}/offsets/final_offsets_WT_mod.txt"
ROI_START="${PLASTID_DIR}/mtDNA_cds_start_rois_mod.txt"

ANNOT_GTF="${PLASTID_DIR}/Human_mtDNA_polyA2_mod2.gtf2"
ANNOT_BED="${PLASTID_DIR}/Human_mtDNA_polyA2_mod3_extended.bed"

# ------------ OUTPUT FOLDERS ------------

mkdir -p "${PLASTID_DIR}/phasing_mod"
mkdir -p "${PLASTID_DIR}/counts_mod"
mkdir -p "${PLASTID_DIR}/vectors_mod"

ERROR_LOG="${PLASTID_DIR}/plastid_recovery_errors.log"

echo "Starting plastid pipeline..."
echo "Started: $(date)" >> "${ERROR_LOG}"

# ------------ LOOP OVER BAM FILES ------------

for bam in "${DATA_DIR}"/*_mapped.sorted.bam
do

    sample=$(basename "${bam}" _mapped.sorted.bam)

    phasing_out="${PLASTID_DIR}/phasing_mod/${sample}"
    counts_out="${PLASTID_DIR}/counts_mod/Counts_${sample}.txt"
    vectors_out="${PLASTID_DIR}/vectors_mod/${sample}"

    # ------------ SKIP COMPLETED ------------
    # Only use get_count_vectors output as completion check

    if [ -d "${vectors_out}" ]; then
        echo "Skipping ${sample} (vectors already exist)"
        continue
    fi

    echo ""
    echo "======================================"
    echo "Processing ${sample}"
    echo "======================================"

    # ------------ PHASING ------------

    phase_by_size \
        "${ROI_START}" \
        "${phasing_out}" \
        --count_files "${bam}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --codon_buffer 5 \
        --min_length 28 \
        --max_length 37

    if [ $? -ne 0 ]; then
        echo "$(date) | phase_by_size failed | ${sample}" >> "${ERROR_LOG}"
        continue
    fi

    # ------------ COUNTS ------------

    counts_in_region \
        "${counts_out}" \
        --count_files "${bam}" \
        --annotation_files "${ANNOT_GTF}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --min_length 28 \
        --max_length 37

    if [ $? -ne 0 ]; then
        echo "$(date) | counts_in_region failed | ${sample}" >> "${ERROR_LOG}"
        continue
    fi

    # ------------ COUNT VECTORS ------------

    get_count_vectors \
        --annotation_files "${ANNOT_BED}" \
        --annotation_format BED \
        --count_files "${bam}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --min_length 28 \
        --max_length 37 \
        "${vectors_out}"

    if [ $? -ne 0 ]; then
        echo "$(date) | get_count_vectors failed | ${sample}" >> "${ERROR_LOG}"
        continue
    fi

    echo "${sample} done"

done

echo ""
echo "======================================"
echo "Plastid pipeline finished"
echo "======================================"

echo ""
echo "Errors logged to:"
echo "${ERROR_LOG}"
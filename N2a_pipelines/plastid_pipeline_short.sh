#!/bin/bash

source /Users/annkrg/Applications/miniconda3/etc/profile.d/conda.sh
conda activate plastid

# Continue even if one sample fails
set +e

# ------------ PATHS ------------

DATA_DIR="processed"
PLASTID_DIR="Alignment_polyA"

ROI_START="${PLASTID_DIR}/mtDNA_cds_start_rois.txt"

ANNOT_GTF="${PLASTID_DIR}/Mouse_mtDNA_polyA.gtf2"
ANNOT_BED="${PLASTID_DIR}/Mouse_polyA_UTR_extended.bed"

# ------------ OUTPUT FOLDERS ------------

mkdir -p "${PLASTID_DIR}/phasing_short"
mkdir -p "${PLASTID_DIR}/counts_short"
mkdir -p "${PLASTID_DIR}/vectors_short"

ERROR_LOG="${PLASTID_DIR}/plastid_recovery_errors.log"

echo "Starting plastid pipeline..."
echo "Started: $(date)" >> "${ERROR_LOG}"

# ------------ LOOP OVER BAM FILES ------------

for bam in "${DATA_DIR}"/*_mapped.sorted.bam
do

    sample=$(basename "${bam}" _mapped.sorted.bam)

    # ----------------------------------------
    # Determine sample-specific offset file
    # ----------------------------------------

    offset_sample=$(echo "${sample}" | sed 's/_cut_aligned$//')

    OFFSET_FILE="${PLASTID_DIR}/offsets/${offset_sample}_p-site_stop_p_offsets.txt"

    if [ ! -f "${OFFSET_FILE}" ]; then
        echo "$(date) | missing offset file | ${sample}" >> "${ERROR_LOG}"
        echo "Missing offset file: ${OFFSET_FILE}"
        continue
    fi

    phasing_out="${PLASTID_DIR}/phasing_short/${sample}"
    counts_out="${PLASTID_DIR}/counts_short/Counts_${sample}.txt"
    vectors_out="${PLASTID_DIR}/vectors_short/${sample}"

    # ------------ SKIP COMPLETED ------------

    if [ -d "${vectors_out}" ]; then
        echo "Skipping ${sample} (vectors already exist)"
        continue
    fi

    echo ""
    echo "======================================"
    echo "Processing ${sample}"
    echo "Offset file: ${OFFSET_FILE}"
    echo "======================================"

    # ------------ PHASING ------------

    phase_by_size \
        "${ROI_START}" \
        "${phasing_out}" \
        --count_files "${bam}" \
        --fiveprime_variable \
        --offset "${OFFSET_FILE}" \
        --codon_buffer 5 \
        --min_length 30 \
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
        --min_length 30 \
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
        --min_length 30 \
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
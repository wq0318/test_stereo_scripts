#!/bin/bash
set -euo pipefail

# Fill these variables on the cloud platform before running.
MASK_FILE="/path/to/barcodeToPos.h5"
RAW_FQ1="/path/to/sample_R1.fq.gz"
RAW_FQ2="/path/to/sample_R2.fq.gz"
REF_DIR="/path/to/reference_dir"
SAMPLE_NAME="sample_lane_or_test_id"
CHR_M="chrM"
OUTPUT_ROOT="/path/to/test_output"

# Step 1
/test_stereo_scripts/tests/test_barcode_mapping.sh \
    "${MASK_FILE}" "${RAW_FQ1}" "${RAW_FQ2}" "${SAMPLE_NAME}" "${OUTPUT_ROOT}"

# Step 2
/test_stereo_scripts/tests/test_mapping.sh \
    "${OUTPUT_ROOT}/01.barcode_mapping/fq/${SAMPLE_NAME}.processed.fq1.gz" \
    "${OUTPUT_ROOT}/01.barcode_mapping/fq/${SAMPLE_NAME}.processed.fq2.gz" \
    "${REF_DIR}" \
    "${SAMPLE_NAME}" \
    "${OUTPUT_ROOT}"

# Step 3
/test_stereo_scripts/tests/test_sam2frag.sh \
    "${OUTPUT_ROOT}/02.mapping/sam/${SAMPLE_NAME}.aligned.sam" \
    "${SAMPLE_NAME}" \
    "${CHR_M}" \
    "${OUTPUT_ROOT}"

# Step 4
/test_stereo_scripts/tests/test_merge_frag.sh \
    "${SAMPLE_NAME}" \
    "${OUTPUT_ROOT}" \
    "${OUTPUT_ROOT}/03.sam2frag/fragments/${SAMPLE_NAME}.fragment.tsv.gz"

# Optional: prepare bin50/bin100 inputs with the real WDL task or your platform binning step first.
# Then run downstream separately for each bin:
# /test_stereo_scripts/tests/test_downstream_analysis.sh "/path/to/${SAMPLE_NAME}_bin50.sorted.fragments.tsv.gz" "${SAMPLE_NAME}_bin50" "50" "${OUTPUT_ROOT}"
# /test_stereo_scripts/tests/test_downstream_analysis.sh "/path/to/${SAMPLE_NAME}_bin100.sorted.fragments.tsv.gz" "${SAMPLE_NAME}_bin100" "100" "${OUTPUT_ROOT}"

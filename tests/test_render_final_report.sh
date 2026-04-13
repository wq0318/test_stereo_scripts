#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage: test_render_final_report.sh <sample_name> <bin_size> <output_root> <barcode_metrics_csv[,more]>"
    exit 1
fi

SAMPLE_NAME=$1
BIN_SIZE=$2
OUTPUT_ROOT=$3
BARCODE_METRICS=$4

REPORT_DIR="${OUTPUT_ROOT}/07.report/bin${BIN_SIZE}"
ANALYSIS_DIR="${OUTPUT_ROOT}/06.analysis/bin${BIN_SIZE}"
HTML_FILE="${REPORT_DIR}/${SAMPLE_NAME}.stereoATAC.report.html"

mkdir -p "${REPORT_DIR}"

python /test_stereo_scripts/reporting/render_final_report.py \
    "${SAMPLE_NAME}" \
    "${BIN_SIZE}" \
    "${HTML_FILE}" \
    "${BARCODE_METRICS}" \
    "${OUTPUT_ROOT}/02.mapping/report/${SAMPLE_NAME}.mapping.metrics.csv" \
    "${OUTPUT_ROOT}/03.sam2frag/report/${SAMPLE_NAME}.sam2frag.metrics.csv" \
    "${OUTPUT_ROOT}/04.merge/report/${SAMPLE_NAME}.merge.metrics.csv" \
    "${ANALYSIS_DIR}/04_metrics/${SAMPLE_NAME}.downstream.metrics.csv" \
    --saturation-plot "${ANALYSIS_DIR}/01_qc/${SAMPLE_NAME}_Saturation.png" \
    --qc-violin "${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_qc_violin.png" \
    --tss-scatter "${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_tss_scatter.png" \
    --fragment-size "${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_fragment_size.png" \
    --spatial-qc "${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_spatial_qc.png" \
    --cluster-plots "${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_clusters_res0_4.png,${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_clusters_res0_6.png,${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_clusters_res0_8.png,${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_clusters_res1.png,${ANALYSIS_DIR}/03_signac/${SAMPLE_NAME}_clusters_res1_2.png"

echo "Final HTML report: ${HTML_FILE}"

#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage: test_downstream_analysis.sh <fragment_tsv_gz> <sample_name> <bin_size> <output_root>"
    exit 1
fi

FRAGMENT_TSV=$1
SAMPLE_NAME=$2
BIN_SIZE=$3
OUTPUT_ROOT=$4

OUTPUT_DIR="${OUTPUT_ROOT}/06.analysis/bin${BIN_SIZE}"
SUMMARY_FILE="${OUTPUT_DIR}/04_metrics/${SAMPLE_NAME}.downstream_summary.tsv"
DOWNSTREAM_METRICS="${OUTPUT_DIR}/04_metrics/${SAMPLE_NAME}.downstream.metrics.csv"

/test_stereo_scripts/analysis/downstream_analysis.sh "${FRAGMENT_TSV}" "${SAMPLE_NAME}" "${OUTPUT_DIR}"

python /test_stereo_scripts/reporting/parse_downstream_analysis_outputs.py \
    "${SUMMARY_FILE}" \
    "${OUTPUT_DIR}/04_metrics/${SAMPLE_NAME}.signac_stats.tsv" \
    "${SAMPLE_NAME}" \
    "${BIN_SIZE}" \
    "${DOWNSTREAM_METRICS}"

echo "downstreamAnalysis done: ${OUTPUT_DIR}"
echo "Metrics: ${DOWNSTREAM_METRICS}"

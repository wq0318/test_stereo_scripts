#!/bin/bash
set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage: test_barcode_mapping.sh <mask_file> <input_fq1> <input_fq2> <sample_name> <output_root>"
    exit 1
fi

MASK_FILE=$1
INPUT_FQ1=$2
INPUT_FQ2=$3
SAMPLE_NAME=$4
OUTPUT_ROOT=$5

OUTPUT_DIR="${OUTPUT_ROOT}/01.barcode_mapping"
OUTPUT_FQ1="${OUTPUT_DIR}/fq/${SAMPLE_NAME}.processed.fq1.gz"
OUTPUT_FQ2="${OUTPUT_DIR}/fq/${SAMPLE_NAME}.processed.fq2.gz"
LOG_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.barcodemap.log"
METRICS_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.barcode_mapping.metrics.csv"

mkdir -p "${OUTPUT_DIR}/fq" "${OUTPUT_DIR}/report"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting barcodeMapping for ${SAMPLE_NAME}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input R1: ${INPUT_FQ1}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input R2: ${INPUT_FQ2}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mask: ${MASK_FILE}"

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/opt/conda/lib:/usr/local/lib"

/ST_BarcodeMap/ST_BarcodeMap-0.0.1 \
    --in "${MASK_FILE}" \
    --in1 "${INPUT_FQ1}" \
    --in2 "${INPUT_FQ2}" \
    --umiStart -1 \
    --out "${OUTPUT_FQ1}" \
    --out2 "${OUTPUT_FQ2}" \
    --mismatch 1 \
    --thread 2 \
    --PEout 2>&1 | tee "${LOG_FILE}"

python /test_stereo_scripts/reporting/parse_barcode_mapping_log.py \
    "${LOG_FILE}" \
    "${SAMPLE_NAME}" \
    "${SAMPLE_NAME}" \
    "${METRICS_FILE}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] barcodeMapping done"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Log: ${LOG_FILE}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Metrics: ${METRICS_FILE}"

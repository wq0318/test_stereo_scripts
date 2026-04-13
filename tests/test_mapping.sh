#!/bin/bash
set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage: test_mapping.sh <processed_fq1> <processed_fq2> <ref_dir> <sample_name> <output_root>"
    exit 1
fi

INPUT_FQ1=$1
INPUT_FQ2=$2
REF_DIR=$3
SAMPLE_NAME=$4
OUTPUT_ROOT=$5

OUTPUT_DIR="${OUTPUT_ROOT}/02.mapping"
SAM_FILE="${OUTPUT_DIR}/sam/${SAMPLE_NAME}.aligned.sam"
SUMMARY_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.chromap_summary.csv"
STDERR_LOG="${OUTPUT_DIR}/report/${SAMPLE_NAME}.chromap.stderr.log"
METRICS_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.mapping.metrics.csv"

mkdir -p "${OUTPUT_DIR}/sam" "${OUTPUT_DIR}/report"

export PATH="/usr/local/bin:${PATH}"
export LD_LIBRARY_PATH="/usr/local/lib:/usr/lib64/:${LD_LIBRARY_PATH:-}"

chromap \
    -l 2000 \
    --low-mem \
    --trim-adapters \
    -x "${REF_DIR}/genome.index" \
    -r "${REF_DIR}/genome.fasta" \
    -1 "${INPUT_FQ1}" \
    -2 "${INPUT_FQ2}" \
    -o "${SAM_FILE}" \
    -t 2 \
    --SAM \
    --summary "${SUMMARY_FILE}" \
    2> "${STDERR_LOG}"

python /test_stereo_scripts/reporting/parse_mapping_summary.py \
    "${SUMMARY_FILE}" \
    "${SAMPLE_NAME}" \
    "${SAMPLE_NAME}" \
    "${METRICS_FILE}"

echo "Mapping done: ${SAM_FILE}"
echo "Summary: ${SUMMARY_FILE}"
echo "Metrics: ${METRICS_FILE}"

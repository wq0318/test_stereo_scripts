#!/bin/bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: test_merge_frag.sh <sample_name> <output_root> <fragment1> [fragment2 ...]"
    exit 1
fi

SAMPLE_NAME=$1
OUTPUT_ROOT=$2
shift 2
FRAGMENTS=("$@")

OUTPUT_DIR="${OUTPUT_ROOT}/04.merge"
MERGED_TSV="${OUTPUT_DIR}/fragments/${SAMPLE_NAME}.merged.fragments.tsv.gz"
REPORT_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.merge_stats.tsv"
METRICS_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.merge.metrics.csv"

mkdir -p "${OUTPUT_DIR}/fragments" "${OUTPUT_DIR}/report"

python /bin/StereATAC/script/merge_fragment.py -i "${FRAGMENTS[@]}" -o "${OUTPUT_DIR}" -s "${SAMPLE_NAME}"
sort -k 1,1V -k 2,2n -k 3,3n "${OUTPUT_DIR}/${SAMPLE_NAME}_fragment.tsv" > "${OUTPUT_DIR}/fragments/${SAMPLE_NAME}.merged.fragments.tsv"
bgzip "${OUTPUT_DIR}/fragments/${SAMPLE_NAME}.merged.fragments.tsv"
tabix -p bed "${MERGED_TSV}"
mv "${OUTPUT_DIR}/${SAMPLE_NAME}_frag_stat.tsv" "${REPORT_FILE}"
rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}_fragment.tsv"

python /test_stereo_scripts/reporting/parse_merge_report.py \
    "${REPORT_FILE}" \
    "${SAMPLE_NAME}" \
    "${METRICS_FILE}"

echo "mergeFrag done: ${MERGED_TSV}"
echo "Metrics: ${METRICS_FILE}"

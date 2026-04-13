#!/bin/bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
    echo "Usage: test_sam2frag.sh <input_sam> <sample_name> <chrM_name> <output_root>"
    exit 1
fi

INPUT_SAM=$1
SAMPLE_NAME=$2
CHR_M=$3
OUTPUT_ROOT=$4

OUTPUT_DIR="${OUTPUT_ROOT}/03.sam2frag"
ALIGNED_BAM="${OUTPUT_DIR}/bam/${SAMPLE_NAME}.aligned.bam"
SORTED_BAM="${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam"
FRAGMENT_TSV="${OUTPUT_DIR}/fragments/${SAMPLE_NAME}.fragment.tsv.gz"
REPORT_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.sam2bam_report.txt"
METRICS_FILE="${OUTPUT_DIR}/report/${SAMPLE_NAME}.sam2frag.metrics.csv"

mkdir -p "${OUTPUT_DIR}/bam" "${OUTPUT_DIR}/fragments" "${OUTPUT_DIR}/report"

/scATAC-seq_pipeline_v2.4/bin/PISA sam2bam \
    -t 4 \
    "${INPUT_SAM}" \
    -o "${ALIGNED_BAM}" \
    -report "${REPORT_FILE}" \
    -mito "${CHR_M}"

samtools sort --threads 4 -o "${SORTED_BAM}" "${ALIGNED_BAM}"
/scATAC-seq_pipeline_v2.4/bin/PISA bam2frag -tag CB -isize 2000 -@ 4 -o "${FRAGMENT_TSV}" "${SORTED_BAM}"

python /test_stereo_scripts/reporting/parse_sam2frag_report.py \
    "${REPORT_FILE}" \
    "${SAMPLE_NAME}" \
    "${SAMPLE_NAME}" \
    "${METRICS_FILE}"

echo "sam2frag done: ${FRAGMENT_TSV}"
echo "Report: ${REPORT_FILE}"
echo "Metrics: ${METRICS_FILE}"

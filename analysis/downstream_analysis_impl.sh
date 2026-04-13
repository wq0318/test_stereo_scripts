#!/bin/bash
set -euo pipefail

SCRIPT_DIR="${SCRIPT_DIR:-/test_stereo_scripts/analysis}"

FRAGMENT=${1:-}
SAMPLE=${2:-}
OUTDIR=${3:-}

if [ -z "$FRAGMENT" ] || [ -z "$SAMPLE" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: downstream_analysis.sh <fragment_file> <sample_name> <output_dir>"
    exit 1
fi

mkdir -p "$OUTDIR/01_qc" "$OUTDIR/02_filtered_data" "$OUTDIR/03_signac" "$OUTDIR/04_metrics"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting downstream analysis for $SAMPLE"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input fragment: $FRAGMENT"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output directory: $OUTDIR"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: saturation analysis"
perl "${SCRIPT_DIR}/atac_saturation.pl" \
    -i "$FRAGMENT" \
    -o "$OUTDIR/01_qc/${SAMPLE}"

mv "$OUTDIR/01_qc/${SAMPLE}_result.txt" "$OUTDIR/01_qc/${SAMPLE}_saturation_stats.tsv" 2>/dev/null || true
mv "$OUTDIR/01_qc/${SAMPLE}.cutoff.txt" "$OUTDIR/01_qc/${SAMPLE}_cutoff_value.txt" 2>/dev/null || true
mv "$OUTDIR/01_qc/${SAMPLE}.barcode_dist.tmp" "$OUTDIR/01_qc/${SAMPLE}_barcode_distribution.tsv" 2>/dev/null || true
mv "$OUTDIR/01_qc/${SAMPLE}.barcode_pass_signal.txt" "$OUTDIR/01_qc/${SAMPLE}_barcodes_pass_signal.txt" 2>/dev/null || true

cutoff_file="$OUTDIR/01_qc/${SAMPLE}_cutoff_value.txt"
if [ -f "$cutoff_file" ]; then
    cutoff=$(cat "$cutoff_file")
else
    cutoff=1000
    echo "$cutoff" > "$cutoff_file"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: spatial filter"
Rscript "${SCRIPT_DIR}/draw_atac_qc.R" \
    "$OUTDIR/01_qc/${SAMPLE}_saturation_stats.tsv" \
    "$OUTDIR/01_qc/${SAMPLE}_barcode_distribution.tsv" \
    "$OUTDIR/01_qc/${SAMPLE}" \
    SPATIAL_FILTER \
    "$cutoff"

if [ -f "$OUTDIR/01_qc/${SAMPLE}.valid_barcodes.txt" ]; then
    mv "$OUTDIR/01_qc/${SAMPLE}.valid_barcodes.txt" "$OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt"
else
    cp "$OUTDIR/01_qc/${SAMPLE}_barcodes_pass_signal.txt" "$OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: filter fragments"
VALID_BC="$OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt"
zcat "$FRAGMENT" | \
    awk -v bc_file="$VALID_BC" 'BEGIN{while((getline k < bc_file)>0)bc[k]=1} $1!="#" && $1!="chrM" && bc[$4]==1 {print}' | \
    bgzip -c > "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz"

zcat "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz" | \
    sort -t$'\t' -k1,1 -k2,2n -k3,3n | \
    bgzip -c > "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.sorted.tsv.gz"
mv "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.sorted.tsv.gz" "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz"
tabix -p bed "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 4: Signac analysis"
Rscript "${SCRIPT_DIR}/signac_analysis.R" \
    "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz" \
    "$SAMPLE" \
    "$OUTDIR/03_signac"

if [ -f "$OUTDIR/03_signac/${SAMPLE}_stats.txt" ]; then
    cp "$OUTDIR/03_signac/${SAMPLE}_stats.txt" "$OUTDIR/04_metrics/${SAMPLE}.signac_stats.tsv"
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 5: summarize downstream outputs"
total_barcodes=$(tail -n +2 "$OUTDIR/01_qc/${SAMPLE}_barcode_distribution.tsv" 2>/dev/null | wc -l | awk '{print $1}')
valid_barcodes=$(wc -l < "$OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt")
filtered_fragments=$(zcat "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz" | wc -l | awk '{print $1}')
raw_fragment_size=$(ls -lh "$FRAGMENT" | awk '{print $5}')
filtered_fragment_size=$(ls -lh "$OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz" | awk '{print $5}')
chrM_count=$(pigz -dc "$FRAGMENT" | awk '$1=="chrM"{count++} END{print count+0}')
total_count=$(pigz -dc "$FRAGMENT" | wc -l | awk '{print $1}')
if [ "$total_count" -gt 0 ]; then
    chrM_ratio=$(awk -v chrM="$chrM_count" -v total="$total_count" 'BEGIN { printf "%.2f", chrM * 100 / total }')
else
    chrM_ratio="0.00"
fi
saturation_value=$(tail -n 1 "$OUTDIR/01_qc/${SAMPLE}_saturation_stats.tsv" | awk '{print $4}')

cat > "$OUTDIR/04_metrics/${SAMPLE}.downstream_summary.tsv" <<EOF
metric	value	unit
total_barcodes	${total_barcodes}	barcodes
valid_barcodes	${valid_barcodes}	barcodes
filtered_fragments	${filtered_fragments}	fragments
raw_fragment_size	${raw_fragment_size}	file_size
filtered_fragment_size	${filtered_fragment_size}	file_size
chrM_count	${chrM_count}	fragments
chrM_ratio	${chrM_ratio}	%
saturation	${saturation_value}	%
EOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downstream analysis complete for $SAMPLE"

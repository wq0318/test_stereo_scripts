#!/bin/bash
set +e  # 测试模式：允许部分步骤失败

# 参数设置 - WDL 环境中通过命令行传入
FRAGMENT=${1}
SAMPLE=${2}
OUTDIR=${3}

# 检查必需参数
if [ -z "$FRAGMENT" ] || [ -z "$SAMPLE" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: run_pipeline.sh <fragment_file> <sample_name> <output_dir>"
    exit 1
fi

# 脚本目录（在 Docker 中的固定路径）
SCRIPT_DIR="/test_stereo_scripts"

# 创建输出目录结构
mkdir -p $OUTDIR/{01_qc,02_filtered_data,03_signac,report}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting pipeline for $SAMPLE"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Input: $FRAGMENT"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output: $OUTDIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Directory structure created:"
echo "  - 01_qc/       : Quality control plots and statistics"
echo "  - 02_filtered_data/ : Filtered fragments and barcodes"
echo "  - 03_signac/   : Signac analysis results"
echo "  - report/      : HTML report"

# Step 1: Perl 脚本 - 信号过滤 + 动态 cutoff
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: Running saturation analysis..."
perl ${SCRIPT_DIR}/ATAC_saturation_final.pl \
    -i $FRAGMENT \
    -o $OUTDIR/01_qc/${SAMPLE}

# 移动文件到对应目录，保留sample_name前缀以便区分多个样本
mv $OUTDIR/01_qc/${SAMPLE}_result.txt $OUTDIR/01_qc/${SAMPLE}_saturation_stats.txt 2>/dev/null || true
mv $OUTDIR/01_qc/${SAMPLE}.cutoff.txt $OUTDIR/01_qc/${SAMPLE}_cutoff_value.txt 2>/dev/null || true
mv $OUTDIR/01_qc/${SAMPLE}.barcode_dist.tmp $OUTDIR/01_qc/${SAMPLE}_barcode_distribution.txt 2>/dev/null || true
mv $OUTDIR/01_qc/${SAMPLE}.barcode_pass_signal.txt $OUTDIR/01_qc/${SAMPLE}_barcodes_pass_signal.txt 2>/dev/null || true
# Saturation 图保留 sample_name 在文件名中

# 计算chrM比例（基于原始输入文件，非过滤后）
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Calculating chrM ratio..."
chrM_count=$(pigz -dc $FRAGMENT | grep -c "^chrM" 2>/dev/null || echo 0)
total_count=$(pigz -dc $FRAGMENT | wc -l 2>/dev/null || echo 1)
if [ "$total_count" -gt 0 ]; then
    chrM_ratio=$(echo "scale=2; $chrM_count * 100 / $total_count" | bc)
else
    chrM_ratio="N/A"
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] chrM ratio: ${chrM_ratio}% (${chrM_count}/${total_count})"

# 获取 cutoff 值
cutoff_file="$OUTDIR/01_qc/${SAMPLE}_saturation_stats.txt"
if [ ! -f "$cutoff_file" ]; then
    echo "Warning: Saturation analysis failed, using default cutoff 1000"
    cutoff_file="/dev/null"
    echo -e "1000\t1000\t1000" > $OUTDIR/01_qc/${SAMPLE}_saturation_stats.txt
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: Extracting valid barcodes..."

# Step 2: 读取 cutoff 值
if [ -f "$OUTDIR/01_qc/${SAMPLE}_cutoff_value.txt" ]; then
    cutoff=$(cat "$OUTDIR/01_qc/${SAMPLE}_cutoff_value.txt")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cutoff value: $cutoff"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Warning: Cutoff file not found, using default 1000"
    cutoff=1000
fi

# 过滤 barcode，先进行信号过滤
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2a: Signal filter (cutoff: $cutoff)..."
# 注意：signal filter 结果已由 Perl 脚本生成到 .barcode_pass_signal.txt

# Step 3: R 脚本 - 空间过滤 (网格过滤) + 绘制QC图
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: Running spatial filter and generating QC plots..."
Rscript ${SCRIPT_DIR}/draw_ATAC_refined.R \
    $cutoff_file \
    $OUTDIR/01_qc/${SAMPLE}_barcode_distribution.txt \
    $OUTDIR/01_qc/${SAMPLE} \
    SPATIAL_FILTER \
    $cutoff

# 检查空间过滤结果
if [ ! -f "$OUTDIR/01_qc/${SAMPLE}.valid_barcodes.txt" ]; then
    echo "Warning: Spatial filter did not produce barcode list, using signal filter result"
    mv $OUTDIR/01_qc/${SAMPLE}_barcodes_pass_signal.txt $OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt
else
    mv $OUTDIR/01_qc/${SAMPLE}.valid_barcodes.txt $OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt
    rm $OUTDIR/01_qc/${SAMPLE}_barcodes_pass_signal.txt 2>/dev/null || true
fi

# Step 4: 过滤 fragment 文件 (去除 chrM，只保留 valid barcode)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 4: Filtering fragment file..."
 
VALID_BC="$OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt"

zcat $FRAGMENT | \
    awk -v bc_file="$VALID_BC" 'BEGIN{while((getline k < bc_file)>0)bc[k]=1} \
    $1!="#" && $1!="chrM" && bc[$4]==1 {print}' | \
    bgzip -c > $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtered fragments saved (BGZF format)"

# 为 Signac 建立索引
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing fragment file for Signac..."
# 先解压、排序、再用 bgzip 压缩（BGZF格式，tabix 需要）
zcat $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz | \
    sort -t$'\t' -k1,1 -k2,2n -k3,3n | \
    bgzip -c > $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.sorted.tsv.gz

# 使用排序后的文件
mv $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.sorted.tsv.gz $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz

# 建立 tabix 索引（Signac 需要）
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating tabix index..."
tabix -p bed $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fragment file ready with index"

# Step 5: Signac 下游分析
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 5: Running Signac analysis..."
Rscript ${SCRIPT_DIR}/signac.R \
    $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz \
    $SAMPLE \
    $OUTDIR/03_signac

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Signac analysis complete"

# Step 6: 生成 HTML 报告
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 6: Generating HTML report..."

# 获取统计信息
total_barcodes=$(wc -l < $OUTDIR/01_qc/${SAMPLE}_barcode_distribution.txt)
valid_barcodes=$(wc -l < $OUTDIR/02_filtered_data/${SAMPLE}_valid_barcodes.txt)
filtered_frags=$(zcat $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz | wc -l)

# 从饱和度文件获取饱和度 (最后一行的 Percent_Duplicates)
saturation=$(tail -1 $OUTDIR/01_qc/${SAMPLE}_saturation_stats.txt | awk '{printf "%.1f", $4}')

# fragment 文件大小 (原始)
raw_frag_size=$(ls -lh $FRAGMENT 2>/dev/null | awk '{print $5}')
[ -z "$raw_frag_size" ] && raw_frag_size="N/A"

# fragment 文件大小 (过滤后)
filtered_frag_size=$(ls -lh $OUTDIR/02_filtered_data/${SAMPLE}_filtered_fragments.tsv.gz 2>/dev/null | awk '{print $5}')
[ -z "$filtered_frag_size" ] && filtered_frag_size="N/A"

# 从 Signac stats 读取统计值
# 注意：signac.R 已经将所有值转为百分数，这里直接读取即可
signac_stats_file="$OUTDIR/03_signac/${SAMPLE}_stats.txt"
if [ -f "$signac_stats_file" ]; then
    median_nfrags=$(grep "median_nFrags" $signac_stats_file | awk '{print $2}')
    median_tss=$(grep "median_TSS" $signac_stats_file | awk '{print $2}')
    # FRiP和Blacklist在signac.R中已经乘以100转为百分数，直接读取即可
    median_frip=$(grep "median_FRiP" $signac_stats_file | awk '{print $2}')
    median_blacklist=$(grep "median_blacklist" $signac_stats_file | awk '{print $2}')
else
    median_nfrags="N/A"
    median_tss="N/A"
    median_frip="N/A"
    median_blacklist="N/A"
fi

# 生成 HTML 报告 - 图片转为base64内嵌
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting images to base64..."

# 转换函数：图片转为base64
img_to_base64() {
    local img_path=$1
    if [ -f "$img_path" ]; then
        echo "data:image/png;base64,$(base64 -w 0 $img_path)"
    else
        echo "NOT_FOUND"
    fi
}

# 预转换所有图片
SATURATION_IMG=$(img_to_base64 "$OUTDIR/01_qc/${SAMPLE}_Saturation.png")
QC_VIOLIN_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_qc_violin.png")
TSS_SCATTER_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_tss_scatter.png")
FRAG_SIZE_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_fragment_size.png")
SPATIAL_QC_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_spatial_qc.png")

# 多分辨率聚类图
CLUSTER_04_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_clusters_res0_4.png")
CLUSTER_06_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_clusters_res0_6.png")
CLUSTER_08_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_clusters_res0_8.png")
CLUSTER_10_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_clusters_res1.png")
CLUSTER_12_IMG=$(img_to_base64 "$OUTDIR/03_signac/${SAMPLE}_clusters_res1_2.png")

# 处理图片不存在的情况：将NOT_FOUND替换为空字符串
# 同时验证base64数据是否有效（必须以data:image开头）
[[ ! "$SATURATION_IMG" =~ ^data:image ]] && SATURATION_IMG=""
[[ ! "$QC_VIOLIN_IMG" =~ ^data:image ]] && QC_VIOLIN_IMG=""
[[ ! "$TSS_SCATTER_IMG" =~ ^data:image ]] && TSS_SCATTER_IMG=""
[[ ! "$FRAG_SIZE_IMG" =~ ^data:image ]] && FRAG_SIZE_IMG=""
[[ ! "$SPATIAL_QC_IMG" =~ ^data:image ]] && SPATIAL_QC_IMG=""
[[ ! "$CLUSTER_04_IMG" =~ ^data:image ]] && CLUSTER_04_IMG=""
[[ ! "$CLUSTER_06_IMG" =~ ^data:image ]] && CLUSTER_06_IMG=""
[[ ! "$CLUSTER_08_IMG" =~ ^data:image ]] && CLUSTER_08_IMG=""
[[ ! "$CLUSTER_10_IMG" =~ ^data:image ]] && CLUSTER_10_IMG=""
[[ ! "$CLUSTER_12_IMG" =~ ^data:image ]] && CLUSTER_12_IMG=""

# 生成 HTML 报告 - 新设计，base64数据存储在JS中
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating HTML report..."

# 从 report.tsv 读取测序和比对报告数据
# 这个文件来自WDL的report task，包含合并后的统计信息
# OUTDIR结构: ./06.analysis/bin100，需要 ../../../05.report/report.tsv
REPORT_TSV="${OUTDIR}/../../../05.report/report.tsv"
# 如果找不到，尝试其他可能的路径
if [ ! -f "$REPORT_TSV" ]; then
    REPORT_TSV="${OUTDIR}/../../05.report/report.tsv"
fi
if [ ! -f "$REPORT_TSV" ]; then
    REPORT_TSV="./05.report/report.tsv"
fi

# 报告文件名包含sample_name
REPORT_FILE="$OUTDIR/report/${SAMPLE}_stereoATAC_report.html"

# 写入HTML头部
cat > "$REPORT_FILE" << 'HTMLEOF'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
HTMLEOF

echo "    <title>${SAMPLE} - Stereo-ATAC Report</title>" >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'HTMLEOF'
    <style>
        :root {
            --primary: #3b82f6;
            --primary-dark: #2563eb;
            --secondary: #64748b;
            --success: #10b981;
            --warning: #f59e0b;
            --danger: #ef4444;
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --text: #334155;
            --text-light: #64748b;
            --border: #e2e8f0;
            --shadow: 0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.04);
            --shadow-lg: 0 4px 6px -1px rgba(0,0,0,0.08), 0 2px 4px -1px rgba(0,0,0,0.04);
        }
        
        * { box-sizing: border-box; margin: 0; padding: 0; }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1000px;
            margin: 0 auto;
            padding: 1.5rem;
        }
        
        /* Header */
        .header {
            background: linear-gradient(135deg, #1e40af 0%, #3b82f6 100%);
            color: white;
            padding: 2rem 1.5rem;
            border-radius: 12px;
            margin-bottom: 1.5rem;
            box-shadow: var(--shadow-lg);
        }
        .header h1 {
            font-size: 1.75rem;
            font-weight: 600;
            margin-bottom: 0.4rem;
        }
        .header .subtitle {
            opacity: 0.85;
            font-size: 0.9rem;
        }
        
        /* Cards */
        .card {
            background: var(--card-bg);
            border-radius: 10px;
            padding: 1.25rem;
            margin-bottom: 1.25rem;
            box-shadow: var(--shadow);
            border: 1px solid var(--border);
        }
        .card-title {
            font-size: 1.1rem;
            font-weight: 600;
            color: var(--text);
            margin-bottom: 1rem;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border);
        }
        
        /* Stats Table */
        .stats-table {
            width: 100%;
            max-width: 480px;
            margin: 0 auto;
        }
        .stats-table table {
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9rem;
        }
        .stats-table th {
            background: var(--bg);
            padding: 10px 12px;
            text-align: left;
            font-weight: 600;
            color: var(--text);
            border-bottom: 2px solid var(--border);
        }
        .stats-table td {
            padding: 8px 12px;
            border-bottom: 1px solid var(--border);
        }
        .stats-table tr:hover {
            background: var(--bg);
        }
        .stats-table .metric-name {
            color: var(--text-light);
            font-weight: 500;
            width: 60%;
        }
        .stats-table .metric-value {
            font-weight: 500;
            color: var(--text);
            text-align: right;
            width: 45%;
        }
        
        /* Section Headers */
        .section-header {
            font-size: 1.25rem;
            font-weight: 600;
            color: var(--text);
            margin: 1.5rem 0 0.75rem;
            padding-left: 0.75rem;
            border-left: 3px solid var(--primary);
        }
        
        /* Image Containers */
        .image-container {
            background: var(--card-bg);
            border-radius: 10px;
            padding: 1.25rem;
            margin-bottom: 1.25rem;
            box-shadow: var(--shadow);
            border: 1px solid var(--border);
            text-align: center;
            min-height: 150px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
        }
        .image-container h3 {
            font-size: 0.9rem;
            color: var(--text-light);
            margin-bottom: 0.75rem;
            font-weight: 500;
            width: 100%;
        }
        .image-container img {
            max-width: 100%;
            height: auto;
            border-radius: 6px;
            display: block;
            margin: 0 auto;
        }
        .image-placeholder {
            padding: 3rem 1.5rem;
            background: var(--bg);
            border-radius: 6px;
            color: var(--text-light);
            font-style: italic;
            display: flex;
            align-items: center;
            justify-content: center;
            min-height: 120px;
        }
        
        /* Metrics Grid - Two Column Layout */
        .metrics-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1.5rem;
        }
        @media (max-width: 768px) {
            .metrics-grid { grid-template-columns: 1fr; }
        }
        .metrics-column {
            min-width: 0;
        }
        .metrics-column h4 {
            font-size: 0.95rem;
            font-weight: 600;
            color: var(--primary);
            margin-bottom: 0.75rem;
            padding-bottom: 0.4rem;
            border-bottom: 2px solid var(--border);
        }
        .metrics-column .stats-table {
            max-width: 100%;
        }
        
        /* Grid Layouts */
        .grid-2 {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 1.25rem;
        }
        @media (max-width: 768px) {
            .grid-2 { grid-template-columns: 1fr; }
            .container { padding: 1rem; max-width: 100%; }
            .header { padding: 1.5rem; }
            .header h1 { font-size: 1.5rem; }
        }
        
        /* Slider */
        .slider-container {
            background: var(--card-bg);
            border-radius: 10px;
            padding: 1rem 1.25rem;
            margin: 0 0 1rem auto;
            max-width: 320px;
            box-shadow: var(--shadow);
            border: 1px solid var(--border);
        }
        .slider-label {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 0.75rem;
            font-weight: 500;
        }
        .res-badge {
            background: var(--primary);
            color: white;
            padding: 0.2rem 0.6rem;
            border-radius: 12px;
            font-size: 0.85rem;
            font-weight: 600;
        }
        .slider-wrapper {
            position: relative;
            margin: 0.75rem 0;
        }
        input[type="range"] {
            width: 100%;
            height: 6px;
            border-radius: 3px;
            background: linear-gradient(to right, var(--primary) 0%, var(--primary) 50%, var(--border) 50%, var(--border) 100%);
            outline: none;
            -webkit-appearance: none;
            cursor: pointer;
        }
        input[type="range"]::-webkit-slider-thumb {
            -webkit-appearance: none;
            width: 18px;
            height: 18px;
            border-radius: 50%;
            background: white;
            border: 2px solid var(--primary);
            cursor: pointer;
            box-shadow: var(--shadow);
        }
        .slider-labels {
            display: flex;
            justify-content: space-between;
            margin-top: 0.5rem;
            font-size: 0.75rem;
            color: var(--text-light);
        }
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
HTMLEOF

echo "            <h1>${SAMPLE}</h1>" >> "$REPORT_FILE"
echo "            <div class=\"subtitle\">Stereo-ATAC Analysis Report • Generated on $(date '+%Y-%m-%d %H:%M:%S')</div>" >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'HTMLEOF'
        </div>
        
        <!-- Metrics -->
        <div class="card">
            <div class="card-title">📊 Quality Control Metrics</div>
            <div class="metrics-grid">
                <!-- Left Column: Sequencing & Alignment Report -->
                <div class="metrics-column">
                    <h4>🔬 Sequencing & Alignment</h4>
                    <div class="stats-table">
                        <table>
                            <tr>
                                <th>Metric</th>
                                <th style="text-align:right">Value</th>
                            </tr>
                            <tr><td class="metric-name">Total Reads</td><td class="metric-value">${total_reads_report:-N/A}</td></tr>
                            <tr><td class="metric-name">Valid CID Reads</td><td class="metric-value">${valid_cid_report:-N/A}</td></tr>
                            <tr><td class="metric-name">Total Mapping Ratio</td><td class="metric-value">${total_mapping_ratio:-N/A}%</td></tr>
                            <tr><td class="metric-name">High Quality Mapping</td><td class="metric-value">${hq_mapping_ratio:-N/A}%</td></tr>
                            <tr><td class="metric-name">chrM Ratio</td><td class="metric-value">${chrm_ratio:-N/A}%</td></tr>
                        </table>
                    </div>
                </div>
                <!-- Right Column: Sample QC Metrics -->
                <div class="metrics-column">
                    <h4>📈 Sample QC Statistics</h4>
                    <div class="stats-table">
                        <table>
                            <tr>
                                <th>Metric</th>
                                <th style="text-align:right">Value</th>
                            </tr>
HTMLEOF

# 写入统计数据
echo "                    <tr><td class=\"metric-name\">Saturation</td><td class=\"metric-value\">${saturation:-N/A}%</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">chrM Ratio</td><td class=\"metric-value\">${chrM_ratio:-N/A}%</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">Median nFrags</td><td class=\"metric-value\">${median_nfrags:-N/A}</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">TSS Enrichment</td><td class=\"metric-value\">${median_tss:-N/A}</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">FRiP</td><td class=\"metric-value\">${median_frip:-N/A}%</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">Blacklist Ratio</td><td class=\"metric-value\">${median_blacklist:-N/A}%</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">Valid Barcodes</td><td class=\"metric-value\">${valid_barcodes:-N/A}</td></tr>" >> "$REPORT_FILE"
echo "                    <tr><td class=\"metric-name\">Filtered Fragments</td><td class=\"metric-value\">${filtered_frag_size:-N/A}</td></tr>" >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'HTMLEOF'
                        </table>
                    </div>
                </div>
            </div>
        </div>
        
        <!-- Saturation -->
        <div class="section-header">Saturation Analysis</div>
        <div class="image-container" id="saturation-container">
            <h3>Library Saturation Curve</h3>
            <div class="image-placeholder">Loading...</div>
        </div>
        
        <!-- QC Violin -->
        <div class="section-header">QC Distribution</div>
        <div class="image-container" id="qc-container">
            <h3>QC Violin Plots</h3>
            <div class="image-placeholder">Loading...</div>
        </div>
        
        <!-- TSS & Fragment Size -->
        <div class="grid-2">
            <div class="image-container" id="tss-container">
                <h3>TSS Enrichment vs Fragments</h3>
                <div class="image-placeholder">Loading...</div>
            </div>
            <div class="image-container" id="frag-container">
                <h3>Fragment Size Distribution</h3>
                <div class="image-placeholder">Loading...</div>
            </div>
        </div>
        
        <!-- Spatial -->
        <div class="section-header">Spatial Visualization</div>
        <div class="image-container" id="spatial-container">
            <h3>Spatial QC</h3>
            <div class="image-placeholder">Loading...</div>
        </div>
        
        <!-- Clustering -->
        <div class="section-header">Clustering Results</div>
        <div class="slider-container">
            <div class="slider-label">
                <span>Select Resolution</span>
                <span class="res-badge" id="resDisplay">0.8</span>
            </div>
            <div class="slider-wrapper">
                <input type="range" id="resSlider" min="0" max="4" value="2" step="1" oninput="updateCluster(this.value)">
            </div>
            <div class="slider-labels">
                <span>0.4</span><span>0.6</span><span>0.8</span><span>1.0</span><span>1.2</span>
            </div>
        </div>
        <div class="image-container" id="cluster-container">
            <div class="image-placeholder">Loading cluster plot...</div>
        </div>
    </div>
    
    <script>
HTMLEOF

# 写入JS图片数据 - 关键：base64数据只存在于JS中，不在HTML标签中
echo "        // 图片数据存储在JS变量中，页面加载后动态插入" >> "$REPORT_FILE"
echo "        const imageData = {" >> "$REPORT_FILE"
echo "            saturation: '${SATURATION_IMG}'," >> "$REPORT_FILE"
echo "            qc: '${QC_VIOLIN_IMG}'," >> "$REPORT_FILE"
echo "            tss: '${TSS_SCATTER_IMG}'," >> "$REPORT_FILE"
echo "            frag: '${FRAG_SIZE_IMG}'," >> "$REPORT_FILE"
echo "            spatial: '${SPATIAL_QC_IMG}'," >> "$REPORT_FILE"
echo "            clusters: ['${CLUSTER_04_IMG}', '${CLUSTER_06_IMG}', '${CLUSTER_08_IMG}', '${CLUSTER_10_IMG}', '${CLUSTER_12_IMG}']" >> "$REPORT_FILE"
echo "        };" >> "$REPORT_FILE"

cat >> "$REPORT_FILE" << 'HTMLEOF'
        const resLabels = ['0.4', '0.6', '0.8', '1.0', '1.2'];
        
        // 加载图片函数 - 动态创建img元素
        function loadImage(containerId, imgData, altText) {
            const container = document.getElementById(containerId);
            if (!container) return;
            
            // 检查是否有有效的base64数据
            if (imgData && imgData.length > 100 && imgData.startsWith('data:image')) {
                // 创建img元素
                const img = document.createElement('img');
                img.src = imgData;
                img.alt = altText;
                img.style.maxWidth = '100%';
                img.style.height = 'auto';
                img.style.borderRadius = '8px';
                img.onerror = function() {
                    container.innerHTML = '<div class="image-placeholder">' + altText + ' Failed to Load</div>';
                };
                // 清空容器并添加图片
                container.innerHTML = '';
                container.appendChild(img);
            } else {
                // 无有效数据，显示占位符
                container.innerHTML = '<div class="image-placeholder">' + altText + ' Not Available</div>';
            }
        }
        
        // 初始化加载所有图片
        window.onload = function() {
            loadImage('saturation-container', imageData.saturation, 'Saturation Plot');
            loadImage('qc-container', imageData.qc, 'QC Violin Plots');
            loadImage('tss-container', imageData.tss, 'TSS Scatter Plot');
            loadImage('frag-container', imageData.frag, 'Fragment Size Distribution');
            
            loadImage('spatial-container', imageData.spatial, 'Spatial QC Plot');
            
            // 加载默认聚类图 (res=0.8, index=2)
            updateCluster(2);
        };
        
        // 更新聚类图
        function updateCluster(index) {
            document.getElementById('resDisplay').textContent = resLabels[index];
            const clusterData = imageData.clusters[index];
            loadImage('cluster-container', clusterData, 'Cluster Plot (res=' + resLabels[index] + ')');
        }
        
        // 更新滑块背景
        document.getElementById('resSlider').addEventListener('input', function(e) {
            const percent = (e.target.value / 4) * 100;
            e.target.style.background = 'linear-gradient(to right, var(--primary) 0%, var(--primary) ' + percent + '%, var(--border) ' + percent + '%, var(--border) 100%)';
        });
    </script>
</body>
</html>
HTMLEOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline completed successfully!"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Results in: $OUTDIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] View report: $REPORT_FILE"

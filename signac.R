#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(Signac)
    library(Seurat)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
    library(Matrix)
    library(tidyr)
    library(dplyr)
    library(SeuratObject)
    library(future)
    library(data.table)
    library(ArchR)
})

# 命令行参数
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
    stop("Usage: Rscript signac.R <filtered_fragments.tsv.gz> <sample_name> <output_dir>")
}

fragpath <- args[1]
sample_name <- args[2]
outdir <- args[3]

# 参数化配置
WORKERS <- 8
MAX_MEM <- 64 * 1024^3  # 64GB
GENOME <- "mm10"  # 可根据物种修改

plan("multisession", workers = WORKERS) 
options(future.globals.maxSize = MAX_MEM)

# 定义颜色向量
pal <- c('#D51F26', '#272E6A', '#208A42', '#89288F', '#F47D2B', 
            '#FEE500', '#8A9FD1', '#C06CAB', '#E6C2DC', '#90D5E4', 
            '#89C75F', '#F37B7D', '#9983BD', '#D24B27', '#3BBCA8', 
            '#6E4B9E', '#0C727C', '#7E1416', '#D8A767', '#3D3D3D', 
            '#371377', '#7700FF', '#9E0142', '#FF0080', '#DC494C', 
            '#F88D51', '#FAD510', '#FFFF5F', '#88CFA4', '#238B45', 
            '#02401B', '#0AD7D3', '#046C9A', '#A2A475', 'grey35')

message(sprintf("[Signac] Processing: %s", sample_name))
message(sprintf("[Signac] Fragment file: %s", fragpath))

# 检查输入文件
if(!file.exists(fragpath)) {
    stop(sprintf("Fragment file not found: %s", fragpath))
}

# 读取 Fragment 文件
total_counts <- CountFragments(fragpath)
cutoff <- 0  # 输入的 fragment 已经过滤过
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB

message(sprintf("[Signac] Barcodes: %d", length(barcodes)))

# 检查并创建索引（如果没有的话）
index_file <- paste0(fragpath, ".tbi")
if (!file.exists(index_file)) {
    message("[Signac] Fragment index not found, creating...")
    IndexFragments(fragment.path = fragpath, verbose = TRUE)
    message("[Signac] Index created")
}

frags <- CreateFragmentObject(path = fragpath, cells = barcodes)
peaks <- CallPeaks(frags)
# 保存 peaks 用于后续分析
saveRDS(peaks, file.path(outdir, paste0(sample_name, "_peaks.rds")))
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)

chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragpath,
    min.cells = 0,
    min.features = 0
)
pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    project = sample_name
)

# 基因组注释 - 根据物种自动选择
if(GENOME == "mm10") {
    library(EnsDb.Mmusculus.v79)
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
} else if(GENOME == "hg38") {
    library(EnsDb.Hsapiens.v86)
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
} else {
    stop(sprintf("Unsupported genome: %s", GENOME))
}
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- GENOME
Annotation(pbmc) <- annotations

other_ids <- rownames(pbmc@meta.data)
pbmc$coor_x <- as.numeric(sapply(other_ids,function(x){strsplit(x,'_')[[1]][1]}))
pbmc$coor_y <- as.numeric(sapply(other_ids,function(x){strsplit(x,'_')[[1]][2]}))


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 20)

# 动态调整 LSI 维度，避免细胞数太少导致错误
n_cells <- ncol(pbmc)
n_features <- nrow(pbmc)
max_svd_dims <- min(n_cells, n_features) - 2  # 必须严格小于 min(nrow, ncol)
if (max_svd_dims < 8) {
    message(sprintf("[Signac] Warning: Low dimensionality (cells=%d, features=%d), adjusting LSI dims", n_cells, n_features))
    n_svd <- max(5, max_svd_dims)  # 至少5个，但不超过限制
    dims_to_use <- 2:max(2, max_svd_dims)
} else {
    n_svd <- 30
    dims_to_use <- 2:10
}
message(sprintf("[Signac] SVD dims: %d, Using dims: %s", n_svd, paste(dims_to_use, collapse=":")))

pbmc <- RunSVD(pbmc, n = n_svd)
pbmc <- RunUMAP(pbmc, dims = dims_to_use, reduction = 'lsi')
pbmc <- FindNeighbors(object = pbmc, dims = dims_to_use, reduction = 'lsi')

# 多分辨率聚类
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)
message(sprintf("[Signac] Clustering with resolutions: %s", paste(resolutions, collapse=", ")))

# 获取默认assay名称，用于构建正确的聚类列名
default_assay <- DefaultAssay(pbmc)
message(sprintf("[Signac] Default assay: %s", default_assay))

for (res in resolutions) {
    res_str <- gsub("\\.", "_", as.character(res))
    # 动态构建聚类列名，如 ATAC_snn_res.0.4
    cluster_col <- paste0(default_assay, "_snn_res.", res)
    
    # 容错：单个分辨率失败不影响其他
    tryCatch({
        # 聚类
        pbmc <- FindClusters(object = pbmc, algorithm = 3, resolution = res, verbose = FALSE)
        
        # 保存聚类结果到新列（避免覆盖）
        pbmc[[paste0("clusters_res_", res_str)]] <- Idents(pbmc)
        
        # 检查列名是否存在
        if (!cluster_col %in% colnames(pbmc@meta.data)) {
            # 尝试备选列名格式
            alt_col <- paste0("seurat_clusters")
            if (alt_col %in% colnames(pbmc@meta.data)) {
                cluster_col <- alt_col
                message(sprintf("[Signac] Warning: Using fallback column '%s' for res=%s", cluster_col, res))
            } else {
                stop(sprintf("Cluster column '%s' not found in metadata", cluster_col))
            }
        }
        
        # 设置当前ident为当前分辨率，用于画图
        Idents(pbmc) <- cluster_col
        
        # 画图 (UMAP 和 空间图都使用当前分辨率的聚类结果)
        p_umap <- DimPlot(pbmc, label = TRUE, cols = pal, pt.size = 0.1, reduction = "umap")
        # 使用 .data[[ 安全引用含点的列名
        p_spatial <- ggplot()+
            geom_tile(data=pbmc@meta.data,aes(x=coor_x,y=coor_y,fill=.data[[cluster_col]]))+
            scale_fill_manual(values = pal)+
            theme_void()+
            coord_fixed()+
            theme(text=element_text(size=16))
        
        cluster_plot <- p_umap + p_spatial
        plot_file <- file.path(outdir, paste0(sample_name, "_clusters_res", res_str, ".png"))
        ggsave(plot_file, cluster_plot, width = 12, height = 6, dpi = 300)
        message(sprintf("[Signac] Cluster plot (res=%s): %s", res, plot_file))
    }, error = function(e) {
        message(sprintf("[Signac] Warning: Failed to generate cluster plot for res=%s: %s", res, e$message))
        # 生成空白占位图
        blank_plot <- ggplot() + 
            annotate("text", x=0.5, y=0.5, label=sprintf("Plot failed\nres=%s", res)) +
            theme_void()
        plot_file <- file.path(outdir, paste0(sample_name, "_clusters_res", res_str, ".png"))
        ggsave(plot_file, blank_plot, width = 6, height = 6, dpi = 300)
    })
}

# 默认使用 res=0.8 作为当前 ident
pbmc <- FindClusters(object = pbmc, algorithm = 3, resolution = 0.8, verbose = FALSE)


# 计算 QC 指标，添加容错处理（测试数据可能太稀疏）
ns_result <- tryCatch({
    NucleosomeSignal(object = pbmc)
}, error = function(e) {
    message("[Signac] Warning: NucleosomeSignal failed: ", e$message)
    pbmc$nucleosome_signal <- NA
    pbmc
})
if (!is.null(ns_result)) pbmc <- ns_result

tss_result <- tryCatch({
    TSSEnrichment(object = pbmc)
}, error = function(e) {
    message("[Signac] Warning: TSSEnrichment failed (possibly due to sparse test data): ", e$message)
    pbmc$TSS.enrichment <- 0
    pbmc
})
if (!is.null(tss_result)) pbmc <- tss_result

bl_result <- tryCatch({
    pbmc$blacklist_ratio <- FractionCountsInRegion(object = pbmc,regions = blacklist_mm10)
    pbmc
}, error = function(e) {
    message("[Signac] Warning: Blacklist ratio calculation failed: ", e$message)
    pbmc$blacklist_ratio <- 0
    pbmc
})
if (!is.null(bl_result)) pbmc <- bl_result

counts_vector <- total_counts$frequency_count
names(counts_vector) <- total_counts$CB
pbmc$nFrags <- counts_vector[Cells(pbmc)]
pbmc <- FRiP(object = pbmc, assay = 'ATAC', total.fragments = 'nFrags')


##指控数据nFrags，TSS，FRiP，blacklist_ratio
p1 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=nFrags),color='black',fill='#5A8100',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=nFrags),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(median(pbmc@meta.data$nFrags))

p2 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=TSS.enrichment),color='black',fill='#178CA4',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=TSS.enrichment),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(round(median(pbmc@meta.data$TSS.enrichment), 2))

p3 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=FRiP*100),color='black',fill='#FFB400',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=FRiP*100),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(paste0(round(median(pbmc@meta.data$FRiP)*100, 1), "%"))

p4 <- ggplot()+
geom_violin(data=pbmc@meta.data,aes(x=orig.ident,y=blacklist_ratio*100),color='black',fill='#FF6C02',width=1,lwd=0.5,show.legend = F)+
geom_boxplot(data=pbmc@meta.data,aes(x=orig.ident,y=blacklist_ratio*100),color='black',fill='white',lwd=0.5,width=0.3,outlier.shape = NA,show.legend = F)+
theme_classic()+
theme(text=element_text(size=16))+
theme(axis.title.x=element_blank())+
ggtitle(paste0(round(median(pbmc@meta.data$blacklist_ratio)*100, 2), "%"))

qc_violin_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
tryCatch({
    ggsave(file.path(outdir, paste0(sample_name, "_qc_violin.png")), qc_violin_plot, width = 16, height = 4, dpi = 300)
    message("[Signac] QC violin plot saved")
}, error = function(e) {
    message("[Signac] Warning: QC violin plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="QC plot failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_qc_violin.png")), blank, width = 16, height = 4, dpi = 300)
})

##片段TSS散点图
tryCatch({
    # 避免ArchR ggPoint尝试创建Rplots.pdf
    pdf(NULL)
    p5 <- ggPoint(
        x = log10(pbmc@meta.data$nFrags), 
        y = pbmc@meta.data$TSS.enrichment, 
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(log10(pbmc@meta.data$nFrags), probs = 0.99, na.rm = TRUE)),
        ylim = c(0, quantile(pbmc@meta.data$TSS.enrichment, probs = 0.99, na.rm = TRUE))
    ) + 
    theme(
        legend.position = "bottom",
        legend.justification = "center"
    )
    dev.off()
    ggsave(file.path(outdir, paste0(sample_name, "_tss_scatter.png")), p5, width = 6, height = 6, dpi = 300)
    message("[Signac] TSS scatter plot saved")
}, error = function(e) {
    message("[Signac] Warning: TSS scatter plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="TSS scatter failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_tss_scatter.png")), blank, width = 6, height = 6, dpi = 300)
})

##片段分布
threads <- 8 
cmd <- paste0("pigz -dc -p ", threads, " ", fragpath, 
  " | awk '$1 != \"chrM\" && $3-$2 <= 750 && $3-$2 > 0 { counts[$3-$2]++ } END { for (i=1; i<=750; i++) print i, counts[i]+0 }'")
size_df <- fread(cmd = cmd, col.names = c("size", "count"))
size_df[, percentage := (count / sum(count)) * 100]
write.csv(size_df, file.path(outdir, paste0(sample_name, "_fragment_size.csv")))

p6 <- ggplot(size_df, aes(x = size, y = percentage)) +
  geom_line(color = "#1f77b4", linewidth = 1) + 
  geom_area(fill = "#1f77b4", alpha = 0.1) + 
  scale_x_continuous(breaks = seq(0, 750, 150), limits = c(0, 750)) +
  labs(
    x = "Fragment Size (bp)",
    y = "Percentage of Fragments (%)"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth =1), 
    panel.grid = element_blank(), 
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  )
p6_formatted <- p6+theme(text=element_text(size=16))
tryCatch({
    ggsave(file.path(outdir, paste0(sample_name, "_fragment_size.png")), p6_formatted, width = 8, height = 5, dpi = 300)
    message("[Signac] Fragment size plot saved")
}, error = function(e) {
    message("[Signac] Warning: Fragment size plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="Fragment size plot failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_fragment_size.png")), blank, width = 8, height = 5, dpi = 300)
})

##空间可视化nFrags，TSS
tryCatch({
    p7 <- ggplot()+
    geom_tile(data=pbmc@meta.data,aes(x=coor_x,y=coor_y,fill=nFrags))+
    scale_fill_gradientn(colours = rainbow(6))+
    theme_void()+
    coord_fixed()+
    theme(text=element_text(size=16))
    p8 <- ggplot()+
    geom_tile(data=pbmc@meta.data,aes(x=coor_x,y=coor_y,fill=TSS.enrichment))+
    scale_fill_gradientn(colours = rainbow(6))+
    theme_void()+
    coord_fixed()+
    theme(text=element_text(size=16))
    
    # 保存 QC 图 (2个正方形子图并排)
    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc.png")), p7 + p8, width = 12, height = 6, dpi = 300)
    message("[Signac] Spatial QC plot saved")
}, error = function(e) {
    message("[Signac] Warning: Spatial QC plot failed: ", e$message)
    blank <- ggplot() + annotate("text", x=0.5, y=0.5, label="Spatial QC failed") + theme_void()
    ggsave(file.path(outdir, paste0(sample_name, "_spatial_qc.png")), blank, width = 12, height = 6, dpi = 300)
})

# 保存 Seurat 对象
saveRDS(pbmc, file.path(outdir, paste0(sample_name, "_filtered_seurat_object.rds")))

# 输出统计值供 HTML 报告使用
stats_file <- file.path(outdir, paste0(sample_name, "_stats.txt"))
stats_df <- data.frame(
    metric = c("median_nFrags", "median_TSS", "median_FRiP", "median_blacklist", "n_cells"),
    value = c(
        round(median(pbmc@meta.data$nFrags, na.rm = TRUE), 0),
        round(median(pbmc@meta.data$TSS.enrichment, na.rm = TRUE), 2),
        round(median(pbmc@meta.data$FRiP, na.rm = TRUE) * 100, 2),
        round(median(pbmc@meta.data$blacklist_ratio, na.rm = TRUE) * 100, 2),
        ncol(pbmc)
    )
)
write.table(stats_df, stats_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
message(sprintf("[Signac] Stats file: %s", stats_file))

message(sprintf("[Signac] Analysis complete! Results saved to: %s", outdir))
message(sprintf("[Signac] - Peaks: %s", file.path(outdir, paste0(sample_name, "_peaks.rds"))))
message(sprintf("[Signac] - QC violin: %s", file.path(outdir, paste0(sample_name, "_qc_violin.png"))))
message(sprintf("[Signac] - TSS scatter: %s", file.path(outdir, paste0(sample_name, "_tss_scatter.png"))))
message(sprintf("[Signac] - Fragment size: %s", file.path(outdir, paste0(sample_name, "_fragment_size.png"))))
message(sprintf("[Signac] - Spatial QC: %s", file.path(outdir, paste0(sample_name, "_spatial_qc.png"))))
message(sprintf("[Signac] - Cluster plots (multi-res): %s_clusters_res*.png", file.path(outdir, sample_name)))
message(sprintf("[Signac] - Seurat object: %s", file.path(outdir, paste0(sample_name, "_filtered_seurat_object.rds"))))

#!/usr/bin/env Rscript
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)

# --- 参数接收 ---
args <- commandArgs(trailingOnly = TRUE)
sat_file  <- args[1] 
dist_file <- args[2] 
prefix    <- args[3]
mode      <- args[4] # FIND_CUTOFF, SPATIAL_FILTER, DRAW

display_name <- basename(prefix)

# --- 模式 1: 寻找拐点 (输出给 Perl) ---
if (mode == "FIND_CUTOFF") {
    d <- fread(dist_file, header=T)[order(-nFrags)][, Rank := 1:.N]
    
    # 限制搜索范围，确保在大约 100 nFrags 以上的区域寻找
    idx_limit <- max(which(d$nFrags >= 100)) 
    if (is.infinite(idx_limit) || idx_limit < 10) idx_limit <- nrow(d)
    
    # 计算 log 后的斜率变化 (dy/dx)，只在限定范围内计算
    dy_dx <- diff(log10(d$nFrags[1:idx_limit])) / diff(log10(d$Rank[1:idx_limit]))
    
    # 寻找斜率最陡的地方 (拐点)
    knee_idx <- which(dy_dx < -10)[1]
    
    final_cutoff <- d$nFrags[knee_idx]
    
    # 如果没找着，给个保底值
    if(is.na(final_cutoff)) final_cutoff <- 200 
    
    # 唯一输出给 Perl 的标准流
    cat(final_cutoff)

} else if (mode == "SPATIAL_FILTER") {
    # --- 模式 2: 空间过滤 + 生成 QC 图 ---
    cutoff_input <- as.numeric(args[5])
    
    # 读取数据
    dist_all <- fread(dist_file, header=T)[order(-nFrags)]
    dist_all[, c("x", "y") := tstrsplit(Barcode, "_", fixed=TRUE, type.convert=TRUE)]
    
    # 1. 拐点过滤
    dist_candidate <- dist_all[nFrags >= cutoff_input]
    
    # 2. 空间过滤 (Grid Filtering) - N >= 7
    dist_candidate[, grid_id := paste0(round(x / 500), "_", round(y / 500))]
    keep_grids <- dist_candidate[, .N, by = grid_id][N >= 7, grid_id]
    dist_filtered <- dist_candidate[grid_id %in% keep_grids]
    
    # 输出 valid barcode list
    output_file <- paste0(prefix, ".valid_barcodes.txt")
    fwrite(dist_filtered[, .(Barcode)], output_file, col.names=FALSE)
    
    message(sprintf("[SPATIAL_FILTER] Input: %d barcodes", nrow(dist_all)))
    message(sprintf("[SPATIAL_FILTER] After signal filter (>= %d): %d barcodes", 
                    cutoff_input, nrow(dist_all[nFrags >= cutoff_input])))
    message(sprintf("[SPATIAL_FILTER] After spatial filter: %d barcodes", 
                    nrow(dist_filtered)))
    message(sprintf("[SPATIAL_FILTER] Output: %s", output_file))
    
    # 生成简单的空间 QC 图（Raw vs Filtered）
    message("[SPATIAL_FILTER] Generating spatial QC plot...")
    
    # 使用 geom_tile + rainbow 配色
    p_raw <- ggplot(dist_all, aes(x=x, y=y, fill=nFrags)) +
        geom_tile(width=100, height=100) + 
        scale_fill_gradientn(colours = rainbow(6),
                             breaks = c(min(dist_all$nFrags, na.rm=TRUE), max(dist_all$nFrags, na.rm=TRUE)),
                             labels = c("low", "high")) +
        coord_fixed() + 
        theme_void() +
        labs(title = paste0("Raw: ", nrow(dist_all), " spots"), fill = NULL)
    
    median_val <- round(median(dist_filtered$nFrags), 0)
    p_filtered <- ggplot(dist_filtered, aes(x=x, y=y, fill=nFrags)) +
        geom_tile(width=100, height=100) + 
        scale_fill_gradientn(colours = rainbow(6),
                             breaks = c(min(dist_filtered$nFrags, na.rm=TRUE), max(dist_filtered$nFrags, na.rm=TRUE)),
                             labels = c("low", "high")) +
        coord_fixed() + 
        theme_void() +
        labs(title = paste0("Filtered: ", nrow(dist_filtered), " spots"), fill = NULL)
    
} else if (mode == "DRAW") {
    # --- 模式 3: 正式绘图 (由 Perl 传入最终 cutoff) ---
    cutoff_input <- as.numeric(args[5])
    
    # 读取所有数据 (Raw)
    dist_all <- fread(dist_file, header=T)[order(-nFrags)]
    dist_all[, c("x", "y") := tstrsplit(Barcode, "_", fixed=TRUE, type.convert=TRUE)]
    
    # 执行双重过滤得到 Cleaned 数据
    dist_candidate <- dist_all[nFrags >= cutoff_input]
    dist_candidate[, grid_id := paste0(round(x / 500), "_", round(y / 500))]
    keep_grids <- dist_candidate[, .N, by = grid_id][N >= 7, grid_id]
    dist_filtered <- dist_candidate[grid_id %in% keep_grids]
    
    # 计算关键统计值
    real_cutoff <- min(dist_filtered$nFrags)  # 过滤后的最小值 (Min Unique Frags)
    median_val  <- round(median(dist_filtered$nFrags), 0)
    
    message(sprintf("[DRAW] Raw spots: %d", nrow(dist_all)))
    message(sprintf("[DRAW] Filtered spots: %d", nrow(dist_filtered)))
    message(sprintf("[DRAW] Real cutoff (min): %d", real_cutoff))
    message(sprintf("[DRAW] Median nFrags: %d", median_val))
    
    # --- 绘图 ---
    
    # p_raw: 所有检测到的点 (Raw Data) - legend只显示high/low，title只放spot数
    p_raw <- ggplot(dist_all, aes(x=x, y=y, fill=nFrags)) +
        geom_tile(width=100, height=100) + 
        scale_fill_gradientn(colours = rainbow(6),
                             breaks = c(min(dist_all$nFrags, na.rm=TRUE), max(dist_all$nFrags, na.rm=TRUE)),
                             labels = c("low", "high")) +
        coord_fixed() + 
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right") +
        labs(title = paste0("Raw Data: ", nrow(dist_all), " spots"), fill = NULL)
    
    # p_cln: 双重过滤后的点 (Cleaned Data) - legend只显示high/low，title只放spot数
    # nFrags计算方式与signac一致：使用过滤后数据的frequency_count
    p_cln <- ggplot(dist_filtered, aes(x=x, y=y, fill=nFrags)) +
        geom_tile(width=100, height=100) + 
        scale_fill_gradientn(colours = rainbow(6),
                             breaks = c(min(dist_filtered$nFrags, na.rm=TRUE), max(dist_filtered$nFrags, na.rm=TRUE)),
                             labels = c("low", "high")) +
        coord_fixed() + 
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right") +
        labs(title = paste0("Cleaned: ", nrow(dist_filtered), " spots"), fill = NULL)
    
    # p_rank: Rank 图，标注 Min Unique Frags (不是 median！)
    p_rank <- ggplot(dist_all, aes(x=1:nrow(dist_all), y=nFrags)) +
        geom_line(color="dodgerblue") + 
        scale_x_log10() + 
        scale_y_log10() +
        geom_vline(xintercept = nrow(dist_filtered), linetype="dashed", color="red") +
        annotate("label", 
                 x = nrow(dist_filtered), 
                 y = real_cutoff, 
                 label = paste0("Min Unique Frags:", real_cutoff, "\nSpots:", nrow(dist_filtered)), 
                 fill = "red", 
                 color = "white", 
                 size = 3, 
                 vjust = -0.5) +
        theme_bw() + 
        labs(x="Rank", y="Unique Fragments")
    
    # 读取饱和度数据
    s <- read.table(sat_file, header=T)
    s$M_Reads <- round(s$Total_Fragments/1e6, 1)
    
    # p_sens: 灵敏度曲线 (Median_nFrags来自Perl计算，与signac一致)
    p_sens <- ggplot(s, aes(x=M_Reads, y=Median_nFrags)) +
        geom_line(color="darkgreen") + 
        geom_point() +
        geom_text(aes(label=round(Median_nFrags, 0)), 
                  size=2.8, vjust=-1.5, color="darkgreen", check_overlap=T) +
        theme_bw() + 
        labs(x="Reads (M)", y="Median Unique Frags")
    
    # p_sat: 饱和度曲线 (max_sat信息放入subtitle或注释)
    max_sat <- round(max(s$Percent_Duplicates), 1)
    p_sat <- ggplot(s, aes(x=M_Reads, y=Percent_Duplicates)) +
        geom_line(color="firebrick") + 
        geom_point() +
        geom_text(aes(label=paste0(round(Percent_Duplicates, 1), "%")), 
                  size=2.8, vjust=-1.5, color="firebrick", check_overlap=T) +
        theme_bw() + 
        labs(x="Reads (M)", y="Saturation (%)")
    
    # 布局合并 (所有子图接近正方形)
    # 左列: 2个空间图(p_raw, p_cln)，右列: 3个统计图(p_rank, p_sens, p_sat)
    # 宽度比 1.2:1，整体比例让每个子图接近正方形
    final <- (p_raw / p_cln | p_rank / p_sens / p_sat) + 
             plot_layout(widths = c(1.2, 1), heights = c(1, 1, 1))
    
    # 保存 (DPI = 500, 宽高比让子图接近正方形)
    ggsave(paste0(prefix, "_Saturation.png"), final, width=14, height=10, dpi=500, type = "cairo")
    message(sprintf("[DRAW] Final plot saved: %s", paste0(prefix, "_Saturation.png")))
    
} else {
    stop(sprintf("Unknown mode: %s. Use FIND_CUTOFF, SPATIAL_FILTER, or DRAW", mode))
}

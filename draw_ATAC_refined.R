#!/usr/bin/env Rscript
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
sat_file  <- args[1]
dist_file <- args[2]
prefix    <- args[3]
mode      <- args[4]

display_name <- basename(prefix)
my_colors <- c("#0E458F","#0F5298","#0E6BA8","#0C86B8","#3399A1","#3B9C9C","#B2C061","#F2CE38","#F2AB38","#F2AB38","#EB7232","#E65B2E","#E14428","#DC2E22","#DB2921","#CC2623")


if (mode == "FIND_CUTOFF") {
    d <- fread(dist_file, header=T)[order(-nFrags)][, Rank := 1:.N]
    
    # 核心修改：将搜索范围放宽到 100 (或者你认为的最低底线)
    # 这样算法就能看到 200-300 附近的斜率变化了
    idx_limit <- max(which(d$nFrags >= 100)) 
    
    # 计算 log 后的斜率变化
    dy_dx <- diff(log10(d$nFrags)) / diff(log10(d$Rank))
    
    # 寻找斜率最陡的地方（拐点）
    # -2.8 或 -3 是一个经验值，可以根据实际效果微调
    knee_idx <- which(dy_dx < -10)[1]
    
    # 输出算出来的 Cutoff，如果没有算出来，输出 200 作为最低保障
    final_cutoff <- d$nFrags[knee_idx]
    if(is.na(final_cutoff)) final_cutoff <- 200 
    
    cat(final_cutoff)
} else {
    cutoff_val <- as.numeric(args[5])
    dist_all <- fread(dist_file, header=T)[order(-nFrags)]
    dist_all[, c("x", "y") := tstrsplit(Barcode, "_", fixed=TRUE, type.convert=TRUE)]
    
    # 过滤逻辑
    dist_tmp <- dist_all[nFrags >= cutoff_val]
    dist_tmp[, grid_id := paste0(round(x / 500), "_", round(y / 500))]
    keep_grids <- dist_tmp[, .N, by = grid_id][N > 7, grid_id]
    dist_filtered <- dist_tmp[grid_id %in% keep_grids]

    
    ##############################
    # --- 1. 过滤逻辑 (支持 200/300 低捕获样本) ---
    dist_tmp <- dist_all[nFrags >= cutoff_val]
    dist_tmp[, grid_id := paste0(round(x / 500), "_", round(y / 500))]
    keep_grids <- dist_tmp[, .N, by = grid_id][N > 7, grid_id]
    dist_filtered <- dist_tmp[grid_id %in% keep_grids]

    # --- 2. 计算展示用的统计值 ---
    # 计算过滤后的 Median nFrags
    median_val <- round(median(dist_filtered$nFrags), 0)
    
    # --- 3. 空间绘图函数 (完全独立 Legend) ---
    plot_spatial <- function(data, title_name, subtitle_name) {
        ggplot(data, aes(x=x, y=y, color=nFrags)) +
            geom_point(size=0.3) + 
            # 使用各自数据的范围，不设 limits 确保 Legend 独立且适配
            scale_color_gradientn(colors = my_colors, trans = "log10") +
            coord_fixed() + 
            scale_y_reverse() + 
            theme_void() +
            theme(
                panel.background = element_rect(fill = 'black', color = 'black'),
                plot.background = element_rect(fill = 'white', color = NA),
                plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5, size = 11),
                legend.position = "right",
                # 强制显示图例标题
                legend.title = element_text(size = 8)
            ) +
            labs(title = title_name, subtitle = subtitle_name, color = "nFrags")
    }

    # 调用绘图
    p_raw <- plot_spatial(dist_all, "Raw Data", 
                          paste("Spots:", nrow(dist_all)))
    
    # 修改这里：subtitle 显示 Median nFrags
    p_cln <- plot_spatial(dist_filtered, "Cleaned Data", 
                          paste0("Spots: ", nrow(dist_filtered), "\nMedian nFrags: ", median_val))

    #########################
    # 统计图
    s <- read.table(sat_file, header=T)
    s$M_Reads <- round(s$Total_Fragments/1e6, 1)
    
    # Rank图带标注
    p_rank <- ggplot(dist_all, aes(x=1:nrow(dist_all), y=nFrags)) +
        geom_line(color="dodgerblue") + scale_x_log10() + scale_y_log10() +
        geom_vline(xintercept = nrow(dist_filtered), linetype="dashed", color="red") +
        annotate("label", x=nrow(dist_filtered), y=cutoff_val, label=paste0("Cutoff:", cutoff_val, "\nSpots:", nrow(dist_filtered)), fill="red", color="white", size=3, hjust=-0.1) +
        theme_bw() + labs(title="Barcode Rank", x="Rank", y="nFrags")

    # 灵敏度/饱和度 (加 expand 防止标签溢出)
    p_sens <- ggplot(s, aes(x=M_Reads, y=Median_nFrags)) +
        geom_line(color="darkgreen") + geom_point() +
        geom_text(aes(label=paste0("(", M_Reads, ",", round(Median_nFrags,0), ")")), size=2.8, vjust=-1.5, color="darkgreen", check_overlap=T) +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + theme_bw() + labs(title="Sensitivity", x="Reads (M)", y="Median nFrags")

    p_sat <- ggplot(s, aes(x=M_Reads, y=Percent_Duplicates)) +
        geom_line(color="firebrick") + geom_point() +
        geom_text(aes(label=paste0("(", M_Reads, ",", round(Percent_Duplicates,1), "%)")), size=2.8, vjust=-1.5, color="firebrick", check_overlap=T) +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + theme_bw() + labs(title="Saturation", x="Reads (M)", y="Saturation (%)")

    # 合并 (压缩中间间距)
    left_panel <- (p_raw / p_cln) + 
                  plot_layout(guides = "keep")  + theme(plot.margin = margin(5, -15, 5, 5))
    right_panel <- (p_rank / p_sens / p_sat) + theme(plot.margin = margin(5, 5, 5, -15))

    final <- (left_panel | right_panel) + plot_layout(widths = c(1.5, 1)) +
             plot_annotation(title = paste("Stereo-ATAC QC Report:", display_name), theme = theme(plot.title = element_text(size=18, face="bold", hjust=0.5)))

    ggsave(paste0(prefix, "_Saturation.pdf"), final, width=15, height=12)
}

# ==============================================================================
# Figure 4: Volcano Shared - Model 2 (WT & KO at 6h)
# ==============================================================================
# Script: fig4_volcano_shared.R
# Output: (none — provides shared DESeq2 model + volcano plotting utilities)
# Description: Pre-calculates WT/KO 6h DESeq2 model for all fig4 volcano scripts.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
})

# --- STEP 2: Load Core Model ---
dds_path <- file.path(DESEQ_DIR, "core_dds.rds")
if (!file.exists(dds_path)) {
    stop("CRITICAL ERROR: core_dds.rds not found. Run 01_core_model.R first.")
}
dds_6h <- readRDS(dds_path)
message("  Loaded pre-calculated Core DESeq2 Model.")

# --- STEP 3: Extract Primary Contrasts ---
message("  Extracting core 6h contrasts for downstream use...")

save_contrast <- function(dds, contrast_vec, filename) {
    res <- results(dds, contrast = contrast_vec, cooksCutoff = DESEQ_COOKS_CUTOFF)
    if (exists("USE_LFC_SHRINK") && USE_LFC_SHRINK) {
        message("  Applying LFC Shrinkage (ashr) to ", filename, "...")
        res <- lfcShrink(dds, contrast = contrast_vec, res = res, type = "ashr")
    }
    res_df <- as.data.frame(res) %>%
        mutate(gene = rownames(res), gene_name = sapply(gene, get_gene_name)) %>%
        arrange(padj)
    write.csv(res_df, file.path(DESEQ_DIR, filename), row.names = FALSE)
    return(res)
}

res_wt_6h        <- save_contrast(dds_6h, c("group", "WT_Tm_6", "WT_YPD_6"), "deg_WT_6h_TmVsYPD.csv")
res_ko_6h        <- save_contrast(dds_6h, c("group", "KO_Tm_6", "KO_YPD_6"), "deg_KO_6h_TmVsYPD.csv")
res_wt_vs_ko_Tm6h <- save_contrast(dds_6h, c("group", "WT_Tm_6", "KO_Tm_6"), "deg_WT_vs_KO_Tm6h.csv")
message("  Contrasts extracted and saved to: ", DESEQ_DIR)


# Create Labeled Volcano Plot (Publication Style)
create_labeled_volcano <- function(res, title, label_genes = NULL,
                                   show_labels = TRUE, show_legend = TRUE) {
    if (is.data.frame(res) && "gene" %in% colnames(res)) {
        df <- res %>%
            mutate(
                gene_name = sapply(gene, get_gene_name),
                significance = case_when(
                    padj < PADJ_CUTOFF & log2FoldChange > LFC_CUTOFF ~ "Up",
                    padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
                    TRUE ~ "NS"
                )
            ) %>%
            filter(!is.na(padj))
    } else {
        df <- as.data.frame(res) %>%
            mutate(
                gene = rownames(res),
                gene_name = sapply(gene, get_gene_name),
                significance = case_when(
                    padj < PADJ_CUTOFF & log2FoldChange > LFC_CUTOFF ~ "Up",
                    padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
                    TRUE ~ "NS"
                )
            ) %>%
            filter(!is.na(padj))
    }

    n_up   <- sum(df$significance == "Up")
    n_down <- sum(df$significance == "Down")


    df$label <- NA
    if (show_labels && !is.null(label_genes)) {
        idx <- df$gene %in% label_genes | df$gene_name %in% label_genes
        df$label[idx] <- df$gene_name[idx]
    }

    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(data = subset(df, significance == "NS"), color = COLORS$UP_DOWN[["NS"]], alpha = VOLCANO$POINT_ALPHA_NS, size = VOLCANO$POINT_SIZE_NS) +
        geom_point(data = subset(df, significance != "NS"), aes(color = significance), alpha = VOLCANO$POINT_ALPHA_SIG, size = VOLCANO$POINT_SIZE_SIG) +
        scale_color_manual(
            values = COLORS$UP_DOWN,
            labels = c(Down = paste0("Down (", n_down, ")"), NS = "NS", Up = paste0("Up (", n_up, ")")),
            breaks = c("Down", "NS", "Up"),
            name = NULL
        )

    if (show_labels && !is.null(label_genes) && any(!is.na(df$label))) {
        p <- p +
            geom_text_repel(
                data = subset(df, !is.na(label)), aes(label = label),
                color = VOLCANO$LABEL_COLOR, fontface = VOLCANO$LABEL_FONTFACE,
                size = VOLCANO$LABEL_SIZE, segment.color = VOLCANO$LABEL_SEGMENT_COLOR,
                segment.linewidth = VOLCANO$LABEL_SEGMENT_SIZE, min.segment.length = 0,
                box.padding = VOLCANO$LABEL_BOX_PADDING, family = PLOT_FONT_FAMILY, max.overlaps = 20
            )
    }

    p <- p +
        geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = VOLCANO$THRESHOLD_LINE_TYPE, color = VOLCANO$THRESHOLD_LINE_COLOR) +
        geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = VOLCANO$THRESHOLD_LINE_TYPE, color = VOLCANO$THRESHOLD_LINE_COLOR) +
        labs(
            title = title,
            x     = expression(log[2] ~ "Fold Change"),
            y     = expression(-log[10] ~ "(padj)")
        ) +
        theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
        theme(
            plot.background  = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position  = if (show_legend) "top" else "none",
            legend.text      = element_text(size = 10),
            plot.title       = element_text(size = PLOT_FONT_SIZE_TITLE, hjust = 0.5, face = "bold"),
            panel.border     = element_rect(color = "black", fill = NA),
            aspect.ratio     = 1
        ) +
        guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))


    if (!is.null(VOLCANO$XLIM)) {
        p <- p + xlim(VOLCANO$XLIM)
    } else {
        x_max <- max(abs(df$log2FoldChange), na.rm = TRUE) * 1.05
        p <- p + xlim(c(-x_max, x_max))
    }
    y_max <- max(-log10(df$padj[!is.na(df$padj) & df$padj > 0])) * 1.1
    if (!is.null(VOLCANO$YLIM)) p <- p + ylim(VOLCANO$YLIM) else p <- p + ylim(c(0, y_max))

    return(p)
}

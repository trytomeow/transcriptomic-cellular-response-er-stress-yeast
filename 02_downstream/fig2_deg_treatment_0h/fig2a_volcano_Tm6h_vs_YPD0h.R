# ==============================================================================
# Figure 2A: Volcano - WT Tm (6h) vs YPD (0h)
# ==============================================================================
# Script: fig2a_volcano_Tm6h_vs_YPD0h.R
# Output: fig2a_volcano_Tm6h_vs_YPD0h.png/pdf,
#         fig2a_volcano_Tm6h_vs_YPD0h_data.csv
# Description: DEG volcano plot comparing Tm treatment at 6h to baseline.
# ==============================================================================

# --- LOCAL CONFIG ---
SHOW_LABELS <- FALSE
SHOW_LEGEND <- TRUE

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

MARKER_GENES <- MARKERS$VALIDATION

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
})

message("Running: fig2a_volcano_Tm6h_vs_YPD0h")

# --- STEP 2: Load Core DESeq2 Model ---
dds_path <- file.path(DESEQ_DIR, "core_dds.rds")
if (!file.exists(dds_path)) {
    stop("CRITICAL ERROR: core_dds.rds not found. Run 01_core_model.R first.")
}
dds <- readRDS(dds_path)
message("  Loaded core DESeq2 model.")

# --- STEP 3: Extract Results ---
res <- results(dds, contrast = c("group", "WT_Tm_6", "WT_YPD_0"), cooksCutoff = DESEQ_COOKS_CUTOFF)
if (exists("USE_LFC_SHRINK") && USE_LFC_SHRINK) {
    message("  Applying LFC shrinkage (ashr)...")
    res <- lfcShrink(dds, contrast = c("group", "WT_Tm_6", "WT_YPD_0"), res = res, type = "ashr")
}

res_df <- as.data.frame(res) %>%
    mutate(
        gene         = rownames(res),
        gene_name    = sapply(gene, get_gene_name),
        significance = case_when(
            padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up",
            padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
            TRUE ~ "NS"
        )
    ) %>%
    filter(!is.na(padj))

write.csv(res_df, file.path(DESEQ_DIR, "fig2a_volcano_Tm6h_vs_YPD0h_data.csv"), row.names = FALSE)
message("  Saved CSV results.")

# --- STEP 4: Prepare Plot Data ---
n_up   <- sum(res_df$significance == "Up")
n_down <- sum(res_df$significance == "Down")

res_df$label <- NA
if (SHOW_LABELS) {
    idx <- res_df$gene %in% MARKER_GENES | res_df$gene_name %in% MARKER_GENES
    res_df$label[idx] <- res_df$gene_name[idx]
}

plot_title <- if (!is.null(PLOT_TITLES$fig2_volcano_Tm6h_vs_YPD0h$main)) {
    PLOT_TITLES$fig2_volcano_Tm6h_vs_YPD0h$main
} else {
    "Volcano Plot: WT Tm (6h) vs WT YPD (0h)"
}

# --- STEP 5: Generate Volcano Plot ---
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(
        data = subset(res_df, significance == "NS"),
        color = COLORS$UP_DOWN[["NS"]], alpha = VOLCANO$POINT_ALPHA_NS, size = VOLCANO$POINT_SIZE_NS
    ) +
    geom_point(
        data = subset(res_df, significance != "NS"),
        aes(color = significance), alpha = VOLCANO$POINT_ALPHA_SIG, size = VOLCANO$POINT_SIZE_SIG
    ) +
    scale_color_manual(
        values = COLORS$UP_DOWN,
        labels = c(
            Down = paste0("Down (", n_down, ")"),
            NS   = "NS",
            Up   = paste0("Up (", n_up, ")")
        ),
        breaks = c("Down", "NS", "Up"),
        name = NULL
    )

if (SHOW_LABELS && any(!is.na(res_df$label))) {
    p <- p +
        geom_text_repel(
            data               = subset(res_df, !is.na(label)),
            aes(label          = label),
            color              = VOLCANO$LABEL_COLOR,
            fontface           = VOLCANO$LABEL_FONTFACE,
            size               = VOLCANO$LABEL_SIZE,
            segment.color      = VOLCANO$LABEL_SEGMENT_COLOR,
            segment.linewidth  = VOLCANO$LABEL_SEGMENT_SIZE,
            min.segment.length = 0,
            box.padding        = VOLCANO$LABEL_BOX_PADDING,
            family             = PLOT_FONT_FAMILY,
            max.overlaps       = 20
        )
}

p <- p +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = VOLCANO$THRESHOLD_LINE_TYPE, color = VOLCANO$THRESHOLD_LINE_COLOR) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = VOLCANO$THRESHOLD_LINE_TYPE, color = VOLCANO$THRESHOLD_LINE_COLOR) +
    labs(
        title = plot_title,
        x     = expression(log[2] ~ "Fold Change"),
        y     = expression(-log[10] ~ "(padj)")
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title       = element_text(size = PLOT_FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
        legend.position  = if (SHOW_LEGEND) "top" else "none",
        legend.text      = element_text(size = 10),
        panel.border     = element_rect(color = "black", fill = NA),
        aspect.ratio     = 1
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

if (!is.null(VOLCANO$XLIM)) {
    p <- p + xlim(VOLCANO$XLIM)
} else {
    x_max <- max(abs(res_df$log2FoldChange), na.rm = TRUE) * 1.05
    p <- p + xlim(c(-x_max, x_max))
}
if (!is.null(VOLCANO$YLIM)) p <- p + ylim(VOLCANO$YLIM)

# --- STEP 6: Save Figure ---
save_figure(p, "fig2a_volcano_Tm6h_vs_YPD0h", width = VOLCANO$WIDTH, height = VOLCANO$HEIGHT, dir = FIG2_DIR)

message("Done: fig2a_volcano_Tm6h_vs_YPD0h")

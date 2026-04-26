# ==============================================================================
# Figure 3G: Mfuzz Elbow Analysis
# ==============================================================================
# Script: fig3g_mfuzz_elbow.R
# Output: fig3g_mfuzz_elbow.png/pdf
# Description: Determines optimal cluster count via Minimum Centroid Distance
#              and Davies-Bouldin Index curves.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(DESeq2)
    library(Mfuzz)
    library(Biobase)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    if (requireNamespace("clValid", quietly = TRUE)) {
        library(clValid)
    }
})

message("Running: fig3g_mfuzz_elbow")

# --- LOCAL CONFIG ---
C_MIN  <- 2
C_MAX  <- 12
M_VALUE  <- NULL
N_STARTS <- 3
USE_TIMEPOINTS <- c(0, 6, 12, 24)

# --- STEP 2: Load VST Data ---
vsd_path <- file.path(DESEQ_DIR, "core_vsd.rds")
if (!file.exists(vsd_path)) {
    stop("CRITICAL ERROR: core_vsd.rds not found. Run 01_core_model.R first.")
}
vsd <- readRDS(vsd_path)
mat <- assay(vsd)

meta_wt_path <- file.path(DESEQ_DIR, "core_meta.rds")
meta_wt      <- readRDS(meta_wt_path)

sel_wt_tm <- with(meta_wt, genotype == "WT" & treatment == "Tm" & time %in% USE_TIMEPOINTS)
sel_0h    <- with(meta_wt, genotype == "WT" & time == 0 & !sel_wt_tm)
sel       <- sel_wt_tm | sel_0h

meta_use <- meta_wt[sel, ]
meta_use$timepoint <- factor(meta_use$time, levels = USE_TIMEPOINTS)

message("  Loaded pre-calculated VST data.")
message("Samples for clustering: ", nrow(meta_use))

# --- STEP 3: Compute Group Means per Timepoint ---
tp_means <- sapply(USE_TIMEPOINTS, function(tp) {
    sel_tp <- rownames(meta_use)[meta_use$time == tp]
    sel_tp <- sel_tp[sel_tp %in% colnames(mat)]
    if (length(sel_tp) == 0) rep(NA_real_, nrow(mat)) else rowMeans(mat[, sel_tp, drop = FALSE])
})
colnames(tp_means) <- paste0(USE_TIMEPOINTS, "h")
tp_means <- tp_means[complete.cases(tp_means), ]
message("Genes for Mfuzz clustering: ", nrow(tp_means))

# --- STEP 4: Select DEGs and Build ExpressionSet ---
deg_files <- list(
    "6h"  = file.path(DESEQ_DIR, "fig3a_volcano_6h_vs_0h_data.csv"),
    "12h" = file.path(DESEQ_DIR, "fig3b_volcano_12h_vs_6h_data.csv"),
    "24h" = file.path(DESEQ_DIR, "fig3c_volcano_24h_vs_12h_data.csv")
)

deg_genes <- character(0)
for (label in names(deg_files)) {
    fpath <- deg_files[[label]]
    if (!file.exists(fpath)) {
        warning("Missing: ", fpath)
        next
    }
    df <- read.csv(fpath)
    if (!"gene" %in% colnames(df)) df$gene <- df[, 1]
    if (!"padj" %in% colnames(df)) df$padj <- 1
    if (!"log2FoldChange" %in% colnames(df)) df$log2FoldChange <- 0
    sig <- df$gene[!is.na(df$padj) & df$padj < PADJ_CUTOFF &
        abs(df$log2FoldChange) > LFC_CUTOFF]
    message(
        "  ", label, ": ", length(sig), " DEGs (padj < ", PADJ_CUTOFF,
        ", |LFC| > ", LFC_CUTOFF, ")"
    )
    deg_genes <- union(deg_genes, sig)
}

available_genes <- intersect(deg_genes, rownames(tp_means))
message("Union DEGs: ", length(deg_genes), " | Complete data: ", length(available_genes))
if (length(available_genes) < 10) stop("Too few DEGs. Run fig3a/b/c first.")

tp_deg <- tp_means[available_genes, ]

eset <- new("ExpressionSet", exprs = tp_deg)
eset <- filter.NA(eset, thres = 0.25)
eset <- fill.NA(eset, mode = "mean")
eset <- standardise(eset)
message("Genes after standardisation: ", nrow(exprs(eset)))

if (is.null(M_VALUE)) {
    M_VALUE <- mestimate(eset)
    message("Auto-estimated fuzzifier: m = ", round(M_VALUE, 3))
}

# --- STEP 5: Trial Clustering — Range of c values ---
message("\nRunning Dmin() for c = ", C_MIN, " to ", C_MAX, " (repeats = ", N_STARTS, ") ...")
set.seed(42)
mcd_values_vec <- Dmin(eset, m = M_VALUE, crange = C_MIN:C_MAX, repeats = N_STARTS, visu = FALSE)

c_range <- C_MIN:C_MAX
for (i in seq_along(c_range)) {
    message("  c = ", c_range[i], " | Min Centroid Distance = ", round(mcd_values_vec[i], 4))
}

# --- STEP 6: Compute Rate of Change ---
mcd_df <- data.frame(
    c   = c_range,
    mcd = mcd_values_vec,
    roc = c(NA, diff(mcd_values_vec)) # rate of change
)

# --- STEP 7: Plot — MCD Elbow and Rate of Change ---
# Suggest optimal c = where |roc| drops below 5% of total range
thresh_roc <- 0.05 * (max(mcd_df$mcd) - min(mcd_df$mcd))
suggested_c <- with(
    mcd_df[!is.na(mcd_df$roc), ],
    c[which(abs(roc) < thresh_roc)[1]]
)
if (is.na(suggested_c)) suggested_c <- ceiling((C_MIN + C_MAX) / 2)

message("Suggested optimal c = ", suggested_c, " (first k where |delta MCD| < 5% of total range)")

p_mcd <- ggplot(mcd_df, aes(x = c, y = mcd)) +
    geom_line(color = COLORS$TREATMENT[["Tm"]], linewidth = 1.2) +
    geom_point(color = COLORS$TREATMENT[["Tm"]], size = 3) +
    geom_vline(
        xintercept = suggested_c,
        linetype = "dashed", color = "gray40",
        linewidth = 0.8
    ) +
    annotate("label",
        x = suggested_c + 0.3,
        y = max(mcd_df$mcd) * 0.95,
        label = paste0("Suggested c = ", suggested_c),
        hjust = 0,
        size = 3.5,
        color = "gray20",
        fill = "white"
    ) +
    scale_x_continuous(breaks = c_range) +
    labs(
        title    = "Mfuzz: Minimum Centroid Distance vs. Number of Clusters",
        subtitle = paste0("WT–Tm time course (0, 6, 12, 24 h) | Union DEGs | m = ", round(M_VALUE, 3)),
        x        = "Number of Clusters (c)",
        y        = "Minimum Centroid Distance"
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title       = element_text(size = PLOT_FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "gray40"),
        axis.title       = element_text(face = "bold"),
        panel.border     = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank()
    )

# Rate of change subplot
p_roc <- ggplot(
    mcd_df[!is.na(mcd_df$roc), ],
    aes(x = c, y = roc)
) +
    geom_col(fill = COLORS$TREATMENT[["YPD"]], color = "white", width = 0.7) +
    geom_hline(
        yintercept = -thresh_roc,
        linetype = "dashed", color = "gray50", linewidth = 0.7
    ) +
    geom_vline(
        xintercept = suggested_c,
        linetype = "dashed", color = "gray40", linewidth = 0.8
    ) +
    scale_x_continuous(breaks = c_range) +
    labs(
        x = "Number of Clusters (c)",
        y = "Delta Min Centroid Distance\n(rate of change)"
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.title       = element_text(face = "bold"),
        panel.border     = element_rect(color = "black", fill = NA),
        panel.grid.minor = element_blank()
    )

use_patchwork <- requireNamespace("patchwork", quietly = TRUE)
use_cowplot   <- requireNamespace("cowplot",   quietly = TRUE)

if (use_patchwork) {
    library(patchwork)
    p_final <- (p_mcd / p_roc) +
        plot_layout(heights = c(2, 1)) +
        plot_annotation(
            title = "Supplementary S4: Mfuzz Cluster Number Optimisation",
            subtitle = paste0(
                "Optimal c = ", suggested_c,
                " based on minimum centroid distance (MCD) elbow"
            ),
            theme = theme(
                plot.title = element_text(
                    face = "bold", hjust = 0.5,
                    size = PLOT_FONT_SIZE_TITLE
                ),
                plot.subtitle = element_text(
                    hjust = 0.5, size = 10,
                    color = "gray40"
                )
            )
        )
} else if (use_cowplot) {
    library(cowplot)
    p_final <- plot_grid(p_mcd, p_roc, nrow = 2, rel_heights = c(2, 1))
} else {
    # Fallback: save individually
    p_final <- p_mcd
    save_figure(p_roc, "fig3g_mfuzz_elbow_roc", width = 8, height = 3, dir = FIG3_DIR)
}

save_figure(p_final, "fig3g_mfuzz_elbow", width = 8, height = 8, dir = FIG3_DIR)

# --- STEP 8: Save MCD Table ---
write.csv(mcd_df, file.path(FIG3_DIR, "fig3g_mfuzz_elbow_data.csv"), row.names = FALSE)
message("Saved MCD table: fig3g_mfuzz_elbow_data.csv")

message("Done: fig3g_mfuzz_elbow | Suggested c = ", suggested_c)

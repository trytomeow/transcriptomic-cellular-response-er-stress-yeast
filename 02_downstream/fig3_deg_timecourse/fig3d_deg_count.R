# ==============================================================================
# Figure 3: DEG Count Stacked Bar Chart (Timecourse)
# ==============================================================================
# Script: fig3d_deg_count.R
# Output: fig3d_deg_count.png/pdf, fig3d_deg_count_summary.csv
# Description: Reads volcano CSVs, counts Up/Down DEGs, and plots a dodged bar chart.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

message("Running: fig3d_deg_count")

# --- STEP 2: Read and Count DEGs ---
csv_files <- list(
    "WT Tm 6h vs WT Tm 0h"   = file.path(DESEQ_DIR, "fig3a_volcano_6h_vs_0h_data.csv"),
    "WT Tm 12h vs WT Tm 6h"  = file.path(DESEQ_DIR, "fig3b_volcano_12h_vs_6h_data.csv"),
    "WT Tm 24h vs WT Tm 12h" = file.path(DESEQ_DIR, "fig3c_volcano_24h_vs_12h_data.csv")
)

count_degs <- function(path, timepoint) {
    if (!file.exists(path)) {
        warning(paste("File not found:", path))
        return(NULL)
    }
    df <- read.csv(path)
    if (!"gene" %in% colnames(df)) df$gene <- df[, 1]

    sig <- df %>% filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_CUTOFF)

    n_up <- sum(sig$log2FoldChange > 0)
    n_down <- sum(sig$log2FoldChange < 0)

    data.frame(
        Timepoint = timepoint,
        Direction = c("UP", "DOWN"),
        Count     = c(n_up, n_down)
    )
}

deg_counts <- bind_rows(lapply(names(csv_files), function(tp) {
    count_degs(csv_files[[tp]], tp)
}))

if (nrow(deg_counts) == 0) {
    stop("No DEG data found. Ensure fig3a/b/c volcano scripts ran successfully.")
}

deg_counts$Timepoint <- factor(deg_counts$Timepoint, levels = c("WT Tm 6h vs WT Tm 0h", "WT Tm 12h vs WT Tm 6h", "WT Tm 24h vs WT Tm 12h"))
deg_counts$Direction <- factor(deg_counts$Direction, levels = c("UP", "DOWN"))

message("DEG Counts:")
for (tp in levels(deg_counts$Timepoint)) {
    sub <- deg_counts %>% filter(Timepoint == tp)
    message(sprintf("  %s: UP=%d, DOWN=%d", tp, sub$Count[sub$Direction == "UP"], sub$Count[sub$Direction == "DOWN"]))
}

write.csv(deg_counts, file.path(FIG3_DIR, "fig3d_deg_count_summary.csv"), row.names = FALSE)

# --- STEP 3: Create Dodged Bar Chart ---
bar_colors <- c("UP" = COLORS$UP_DOWN[["Up"]], "DOWN" = COLORS$UP_DOWN[["Down"]])

p <- ggplot(deg_counts, aes(x = Timepoint, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = Count),
        position = position_dodge(width = 0.8), vjust = -0.5,
        size = 3.5, fontface = "bold", color = "black"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = bar_colors, name = "Direction") +
    labs(
        title = "DEG Counts: Time course",
        x = "Comparison", y = "Number of DEGs"
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background    = element_rect(fill = "transparent", color = NA),
        panel.background   = element_rect(fill = "transparent", color = NA),
        plot.title         = element_text(size = PLOT_FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
        axis.text          = element_text(color = "black"),
        axis.title         = element_text(face = "plain"),
        panel.border       = element_rect(color = "black", fill = NA),
        aspect.ratio       = 1,
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "top",
        legend.title       = element_blank(),
        legend.text        = element_text(size = 10),
        plot.margin        = margin(20, 20, 20, 20)
    )

save_figure(p, "fig3d_deg_count", width = VOLCANO$WIDTH, height = VOLCANO$HEIGHT, dir = FIG3_DIR)

message("Saved: fig3d_deg_count")
message("Done: fig3d_deg_count")

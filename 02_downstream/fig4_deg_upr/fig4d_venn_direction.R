# ==============================================================================
# Figure 4d: Directional Venn Diagrams (WT vs KO at 6h)
# ==============================================================================
# Script: fig4d_venn_direction.R
# Output: fig4d_venn_up.png/pdf, fig4d_venn_down.png/pdf
# Description: Up/Down split Venn diagrams for WT vs KO DEG overlap.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

# --- Load Required Packages ---
suppressPackageStartupMessages({
    library(ggforce)
    library(ggplot2)
    library(dplyr)
})

message("Running: fig4d_venn_direction")

# --- LOCAL CONFIG ---
FILL_COLORS <- c("grey70", "black")
LINE_COLORS <- c("grey40", "black")

# --- Load DESeq2 Results ---
file_wt <- file.path(DESEQ_DIR, "deg_WT_6h_TmVsYPD.csv")
file_ko <- file.path(DESEQ_DIR, "deg_KO_6h_TmVsYPD.csv")

read_split_deg <- function(filepath) {
    if (!file.exists(filepath)) {
        warning("File not found: ", filepath)
        return(list(up = character(0), down = character(0)))
    }
    df <- read.csv(filepath)
    up <- df %>%
        filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange > LFC_CUTOFF) %>%
        pull(gene)
    down <- df %>%
        filter(!is.na(padj), padj < PADJ_CUTOFF, log2FoldChange < -LFC_CUTOFF) %>%
        pull(gene)
    return(list(up = up, down = down))
}

wt_lists <- read_split_deg(file_wt)
ko_lists <- read_split_deg(file_ko)

message("WT 6h: Up=", length(wt_lists$up), ", Down=", length(wt_lists$down))
message("KO 6h: Up=", length(ko_lists$up), ", Down=", length(ko_lists$down))

create_ggforce_venn <- function(wt_genes, ko_genes, name_wt, name_ko,
                                filename, title, fill_colors, line_colors) {
    n_wt   <- length(wt_genes)
    n_ko   <- length(ko_genes)
    n_both <- length(intersect(wt_genes, ko_genes))
    wt_only <- n_wt - n_both
    ko_only <- n_ko - n_both

    message("\nProcessing: ", filename)
    message("  ", name_wt, " only: ", wt_only)
    message("  Both: ", n_both)
    message("  ", name_ko, " only: ", ko_only)

    r <- VENN$CIRCLE_RADIUS
    d <- VENN$CIRCLE_DIST

    circles_df <- data.frame(
        x0    = c(-d / 2, d / 2),
        y0    = c(0, 0),
        r     = c(r, r),
        group = factor(c(name_wt, name_ko), levels = c(name_wt, name_ko))
    )

    left_center  <- -(d / 2 + r) / 2
    right_center <-  (d / 2 + r) / 2

    count_df <- data.frame(
        x     = c(left_center, 0, right_center),
        y     = c(0, 0, 0),
        label = c(wt_only, n_both, ko_only)
    )

    # --- Category labels (top or bottom, from config) ---
    cat_y <- if (VENN$CAT_POSITION == "bottom") -(r + 0.3) else (r + 0.3)
    cat_df <- data.frame(
        x     = c(-d / 2, d / 2),
        y     = c(cat_y, cat_y),
        label = c(name_wt, name_ko)
    )


    # --- Build plot ---
    p <- ggplot() +
        geom_circle(
            data = circles_df,
            aes(x0 = x0, y0 = y0, r = r, fill = group, color = group),
            alpha = VENN$CIRCLE_ALPHA,
            linewidth = VENN$EDGE_SIZE
        ) +
        geom_text(
            data = count_df,
            aes(x = x, y = y, label = label),
            size = VENN$COUNT_SIZE,
            family = PLOT_FONT_FAMILY,
            color = VENN$FONT_COLOR
        ) +
        geom_text(
            data = cat_df,
            aes(x = x, y = y, label = label),
            size = VENN$CAT_SIZE,
            fontface = "bold",
            family = PLOT_FONT_FAMILY,
            color = VENN$FONT_COLOR,
            parse = TRUE
        ) +
        scale_fill_manual(values = setNames(fill_colors, c(name_wt, name_ko))) +
        scale_color_manual(values = setNames(line_colors, c(name_wt, name_ko))) +
        coord_fixed() +
        labs(title = title) +
        theme_void(base_family = PLOT_FONT_FAMILY) +
        theme(
            plot.title = element_text(
                size = VENN$TITLE_SIZE, face = "bold", hjust = 0.5,
                color = PLOT_FONT_COLOR, margin = margin(b = 5)
            ),
            legend.position = "none",
            plot.margin = margin(20, 20, 20, 20)
        )

    save_figure(p, filename, width = VENN$PLOT_WIDTH, height = VENN$PLOT_HEIGHT, dir = FIG4_DIR)
    message("  Saved: ", filename)
}

# ==============================================================================
# Generate Venns
# ==============================================================================

create_ggforce_venn(
    wt_lists$up, ko_lists$up,
    "'WT'", "italic('hac1')*italic(Delta)",
    "fig4d_venn_up",
    bquote(bold("Upregulated Genes: WT vs")~bolditalic("hac1")*bolditalic(Delta)),
    FILL_COLORS, LINE_COLORS
)

create_ggforce_venn(
    wt_lists$down, ko_lists$down,
    "'WT'", "italic('hac1')*italic(Delta)",
    "fig4d_venn_down",
    bquote(bold("Downregulated Genes: WT vs")~bolditalic("hac1")*bolditalic(Delta)),
    FILL_COLORS, LINE_COLORS
)

# --- Save Gene Lists ---
max_len <- max(length(wt_lists$up), length(wt_lists$down), length(ko_lists$up), length(ko_lists$down))
pad <- function(x) c(x, rep(NA, max_len - length(x)))

df_export <- data.frame(
    WT_Up = pad(wt_lists$up), WT_Down = pad(wt_lists$down),
    KO_Up = pad(ko_lists$up), KO_Down = pad(ko_lists$down)
)
write.csv(df_export, file.path(FIG4_DIR, "fig4d_venn_direction_lists.csv"), row.names = FALSE)
message("Saved: fig4d_venn_direction_lists.csv")

message("Done: fig4d_venn_direction")

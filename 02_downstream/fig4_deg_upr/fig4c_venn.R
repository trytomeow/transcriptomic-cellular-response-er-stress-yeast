# ==============================================================================
# Figure 4c: Venn Diagram - WT vs hac1-KO Treatment Response
# ==============================================================================
# Script: fig4c_venn.R
# Output: fig4c_venn_all.png/pdf
# Description: 2-circle Venn (WT vs KO DEGs at 6h) to identify UPR-dependent
#              vs independent genes.
# ==============================================================================

# --- STEP 1: Load Shared Model ---
shared_path <- "fig4_volcano_shared.R"
if (!file.exists(shared_path)) shared_path <- "../fig4_deg_upr/fig4_volcano_shared.R"
source(shared_path)

suppressPackageStartupMessages({
    library(ggforce)
    library(ggplot2)
})

message("Running: fig4c_venn")

# --- LOCAL CONFIG ---
FILL_COLORS <- c("grey70", "black") # WT-grey, KO-black
LINE_COLORS <- c("grey40", "black") # WT-grey, KO-black
VENN_TITLE <- PLOT_TITLES$fig4c_venn$main
VENN_SUBTITLE_BASE <- PLOT_TITLES$fig4c_venn$sub

# --- STEP 2: Select Comparisons ---
wt_file <- file.path(DESEQ_DIR, "deg_WT_6h_TmVsYPD.csv")
ko_file <- file.path(DESEQ_DIR, "deg_KO_6h_TmVsYPD.csv")

if (!file.exists(wt_file) || !file.exists(ko_file)) {
    stop("Required DEG CSV files not found! Run fig4_volcano_shared.R first.")
}

wt_res_df <- read.csv(wt_file)
ko_res_df <- read.csv(ko_file)

# --- STEP 3: Filter DEGs based on VENN_MODE ---
message("Applying Venn Mode: ", VENN_MODE)

get_degs <- function(res_data, mode) {
    if (!("significance" %in% colnames(res_data))) {
        res_data <- res_data %>% mutate(significance = case_when(
            padj < PADJ_CUTOFF & log2FoldChange > LFC_CUTOFF ~ "Up",
            padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
            TRUE ~ "NS"
        ))
    }

    res_data %>%
        dplyr::filter(
            if (mode == "up") {
                significance == "Up"
            } else if (mode == "down") {
                significance == "Down"
            } else {
                significance %in% c("Up", "Down")
            }
        ) %>%
        dplyr::pull(gene) %>%
        unique()
}

wt_degs <- get_degs(wt_res_df, VENN_MODE)
ko_degs <- get_degs(ko_res_df, VENN_MODE)

message("  - WT DEGs (", VENN_MODE, "): ", length(wt_degs))
message("  - KO DEGs (", VENN_MODE, "): ", length(ko_degs))

# --- STEP 4: Calculate Overlaps ---
if (tolower(VENN_MODE) == "all") {
    wt_up   <- get_degs(wt_res_df, "up")
    wt_down <- get_degs(wt_res_df, "down")
    ko_up   <- get_degs(ko_res_df, "up")
    ko_down <- get_degs(ko_res_df, "down")

    upr_dependent   <- unique(c(setdiff(wt_up, ko_up), setdiff(wt_down, ko_down)))
    upr_independent <- unique(c(intersect(wt_up, ko_up), intersect(wt_down, ko_down)))
    ko_specific     <- unique(c(setdiff(ko_up, wt_up), setdiff(ko_down, wt_down)))
} else {
    upr_dependent   <- setdiff(wt_degs, ko_degs)
    upr_independent <- intersect(wt_degs, ko_degs)
    ko_specific     <- setdiff(ko_degs, wt_degs)
}

message("  - upr-dependent (WT only): ", length(upr_dependent))
message("  - upr-independent (Both): ", length(upr_independent))
message("  - KO-specific (KO only): ", length(ko_specific))

# --- STEP 5: Build ggforce Venn (WT=LEFT, KO=RIGHT always) ---
r <- VENN$CIRCLE_RADIUS
d <- VENN$CIRCLE_DIST

name_wt <- "'WT (Tm vs YPD)'"
name_ko <- "italic('hac1')*italic(Delta)~'(Tm vs YPD)'"

circles_df <- data.frame(
    x0    = c(-d / 2, d / 2),
    y0    = c(0, 0),
    r     = c(r, r),
    group = factor(c(name_wt, name_ko), levels = c(name_wt, name_ko))
)

# Count labels centered in each region
left_center <- -(d / 2 + r) / 2
right_center <- (d / 2 + r) / 2

count_df <- data.frame(
    x     = c(left_center, 0, right_center),
    y     = c(0, 0, 0),
    label = c(length(upr_dependent), length(upr_independent), length(ko_specific))
)

# Category labels (position from config)
cat_y <- if (VENN$CAT_POSITION == "bottom") -(r + 0.3) else (r + 0.3)
cat_df <- data.frame(
    x     = c(-d / 2, d / 2),
    y     = c(cat_y, cat_y),
    label = c(name_wt, name_ko)
)


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
    scale_fill_manual(values = setNames(FILL_COLORS, c(name_wt, name_ko))) +
    scale_color_manual(values = setNames(LINE_COLORS, c(name_wt, name_ko))) +
    coord_fixed() +
    labs(title = VENN_TITLE) +
    theme_void(base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.title = element_text(
            size = VENN$TITLE_SIZE, face = "bold", hjust = 0.5,
            color = PLOT_FONT_COLOR, margin = margin(b = 5)
        ),
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20)
    )

# --- STEP 6: Save ---
save_figure(p, "fig4c_venn", width = VENN$PLOT_WIDTH, height = VENN$PLOT_HEIGHT, dir = FIG4_DIR)

# --- STEP 7: Save Gene Lists ---
gene_lists_df <- data.frame(
    Category = c(
        rep("upr_Dependent", length(upr_dependent)),
        rep("upr_Independent", length(upr_independent)),
        rep("KO_Specific", length(ko_specific))
    ),
    ORF = c(upr_dependent, upr_independent, ko_specific)
)
gene_lists_df$Gene <- sapply(gene_lists_df$ORF, get_gene_name)

write.csv(gene_lists_df, file.path(FIG4_DIR, "fig4c_venn_gene_lists.csv"), row.names = FALSE)
message("  Saved Gene Lists CSV")

message("Done: fig4c_venn")

# ==============================================================================
# Config: Master Configuration File
# ==============================================================================
# Script: 00_config.R
# Output: (none — sourced by all downstream scripts)
# Description: Global settings for data paths, thresholds, plot styling,
#              sample exclusion, and execution flags.
# ==============================================================================

# ==============================================================================
# 1. GLOBAL EXECUTION CONTROL
# ==============================================================================

# --- Sample Exclusion ---
# Set to TRUE to remove outlier samples defined in SAMPLES_TO_EXCLUDE.
# Note on WT-0 baseline: All 0h samples have Treatment=YPD in metadata, so the group
# factor encodes them as WT_YPD_0 / KO_YPD_0. There are no WT_Tm_0 samples.
# For Fig 2 (treatment effect): WT_Tm_Xh vs WT_YPD_0 is the correct contrast.
# For Fig 3 (timecourse): WT_Tm_Xh vs WT_YPD_0 is also correct — WT-0 IS the
# Tm baseline (same physical sample, different experimental labelling context).
USE_EXCLUDED_SAMPLES <- TRUE  # TRUE = apply SAMPLES_TO_EXCLUDE filter; FALSE = use all samples

# --- Analysis Thresholds ---
PADJ_CUTOFF <- 0.05

# --- DESeq2 Settings ---
USE_LFC_SHRINK      <- FALSE  # TRUE = lfcShrink(ashr); FALSE = raw MLE

# cooksCutoff: TRUE (default outlier filtering) or Inf (disable filtering)
DESEQ_COOKS_CUTOFF <- TRUE

# Log2 Fold Change cutoff (Absolute value)
LFC_CUTOFF <- 1.0

GLOBAL_PLOT_WIDTH  <- 10
GLOBAL_PLOT_HEIGHT <- 8

# --- Venn Diagram Mode ---
# Options: "all" (All DEGs), "up" (Upregulated only), "down" (Downregulated only)
VENN_MODE <- "all"

# --- Fig 3D: Exclude WT-YPD at 12h ---
# TRUE: Drop the 12h YPD control point (YPD shows only 0h→6h→24h)
# FALSE: Keep all timepoints for both YPD and Tm
EXCLUDE_YPD_12H <- TRUE

# --- Fig 3H: Mfuzz Membership Threshold ---
# Genes below this value are faded in line plots. Range: 0–1. Default = 0.7.
MFUZZ_MIN_MEMBERSHIP <- 0.7

# --- GO Enrichment Plot Mode ---
# Options: "GeneRatio" (Proportion of genes) or "Count" (Number of genes)
GO_PLOT_X_AXIS <- "GeneRatio"

# --- GO Enrichment Cutoffs ---
GO_PVALUE_CUTOFF <- 0.05 # p-value threshold for GO enrichment
GO_QVALUE_CUTOFF <- 0.05 # q-value (adjusted p-value) threshold
GO_SIMPLIFY_CUTOFF <- 0.7 # cutoff for GO term simplification

# --- GO Top-N per Script ---
TOP_N_GO <- list(
    fig2e_treatment     = 10, # fig2e_go_treatment.R
    fig3e_temporal      = 10, # fig3e_go_temporal.R
    fig3i_mfuzz         = 10, # fig3i_mfuzz_go.R (per cluster)
    fig3i_mfuzz_summary = 40, # fig3i_mfuzz_go.R (summary heatmap)
    fig4e_direction     = 10 # fig4e_go_direction.R
)

# --- GO Global Settings ---
GO_ONTOLOGY <- "BP" # GO ontology: "BP", "MF", or "CC"

# --- KEGG Top-N per Script ---
TOP_N_KEGG <- list(
    fig2f_treatment      = 10, # fig2f_kegg_treatment.R
    fig3f_temporal       = 10, # fig3f_kegg_temporal.R
    fig3j_mfuzz          = 10, # fig3j_mfuzz_kegg.R (per cluster dot plot)
    fig3j_mfuzz_summary  = 40, # fig3j_mfuzz_kegg.R (summary heatmap)
    fig4f_upr            = 10 # fig4f_kegg_upr_direction.R
)

# --- KEGG Cutoffs ---
KEGG_PVALUE_CUTOFF <- 0.05
KEGG_QVALUE_CUTOFF <- 0.05

# --- Enrichment Map Network (emapplot) Settings ---
EMAPPLOT <- list(
    SHOW_CATEGORY = 30, # Max GO terms to display as nodes
    LAYOUT = "nicely", # Network layout: "nicely", "kk", "fr"
    SIMPLIFY_CUTOFF = 0.7, # Cutoff for simplified version
    WIDTH = GLOBAL_PLOT_WIDTH, HEIGHT = GLOBAL_PLOT_HEIGHT # Plot dimensions (inches)
)

GO_PLOT <- list(
    WIDTH = GLOBAL_PLOT_WIDTH,
    HEIGHT = GLOBAL_PLOT_HEIGHT,
    DOT_SIZE_RANGE = c(3, 8),
    FONT_FACE_TITLE = "bold",
    FONT_FACE_AXIS_X = "bold",
    FONT_FACE_AXIS_Y = "plain",
    FONT_FACE_STRIP = "bold",
    X_EXPAND = 0.6,
    Y_EXPAND = 0.05
)

# --- Execution Flags ---
EXECUTE_SECTIONS <- list(
    fig1_qc = TRUE,
    fig2_deg_treatment_0h = TRUE,
    fig3_deg_timecourse = TRUE,
    fig4_deg_upr = TRUE
)

# ==============================================================================
# 2. FILE PATHS & DIRECTORIES
# ==============================================================================

# --- Base Directory ---
# Auto-detected from this file's location when loaded via source().
# This file lives at: <project_root>/02_downstream/fig0_config/00_config.R
# Two dirname() calls walk up to <project_root>.
#
# If auto-detection fails, set BASE_DIR manually and uncomment the line below:
#   BASE_DIR <- "/path/to/repo"
BASE_DIR <- tryCatch({
    .ofile <- sys.frame(sys.nframe())$ofile
    if (is.null(.ofile)) stop("ofile is NULL")
    normalizePath(file.path(dirname(.ofile), "..", ".."), mustWork = FALSE)
}, error = function(e) {
    stop(
        "[CONFIG ERROR] Cannot auto-detect project root.\n",
        "  Set BASE_DIR manually in fig0_config/00_config.R:\n",
        "    BASE_DIR <- \"/path/to/repo\""
    )
})

# --- Input Directories ---
DATA_DIR      <- file.path(BASE_DIR, "data")
COUNTS_FILE   <- file.path(DATA_DIR, "gene_counts.csv")
METADATA_FILE <- file.path(DATA_DIR, "metadata.csv")
SCRIPTS_DIR   <- file.path(BASE_DIR, "02_downstream")

# --- Output Directories ---
OUTPUT_DIR <- file.path(BASE_DIR, "results")
DESEQ_DIR  <- file.path(OUTPUT_DIR, "deseq2_results")

# --- Output Sub-directories ---
DIRS <- list(
    FIG1 = file.path(OUTPUT_DIR, "fig1_qc"),
    FIG2 = file.path(OUTPUT_DIR, "fig2_deg_treatment_0h"),
    FIG3 = file.path(OUTPUT_DIR, "fig3_deg_timecourse"),
    FIG4 = file.path(OUTPUT_DIR, "fig4_deg_upr")
)

FIG1_DIR <- DIRS$FIG1
FIG2_DIR <- DIRS$FIG2
FIG3_DIR <- DIRS$FIG3
FIG4_DIR <- DIRS$FIG4

# Default output dir for save_figure() when dir= is omitted
FIGURE_DIR <- OUTPUT_DIR

# Create all directories if they don't exist
for (d in c(OUTPUT_DIR, DESEQ_DIR, unlist(DIRS))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}


# ==============================================================================
# 3. SAMPLE EXCLUSION LIST
# ==============================================================================
# Specific samples identified as outliers or technical failures
# Format: replicate-genotype-treatment-time
SAMPLES_TO_EXCLUDE <- c(
    "1-WT-Tm-6",    "1-hac1-YPD-6", "1-WT-YPD-24", "3-WT-0",
    "1-WT-Tm-12",   "1-WT-YPD-12",  "2-WT-YPD-12", "3-WT-YPD-12"
)

# ==============================================================================
# 4. PLOT STYLING & AESTHETICS
# ==============================================================================

# --- Global Style Settings ---
PLOT_FONT_FAMILY <- "sans" # "sans" is safe cross-platform (Arial-like)
PLOT_FONT_SIZE_TITLE <- 12
PLOT_FONT_SIZE_LABEL <- 10
PLOT_FONT_COLOR <- "black" # Global font color for plot text
PLOT_FILE_TYPE <- "both" # Options: "png", "pdf", "both"

# --- Color Palettes ---
# Consistent colors across all figures
COLORS <- list(
    TREATMENT = c(YPD = "#3498db", Tm = "#e74c3c"), # Blue / Red
    GENOTYPE = c(WT = "#2ecc71", KO = "#9b59b6"), # Green / Purple
    TIME = c("0h" = "#f1c40f", "6h" = "#e67e22", "12h" = "#e74c3c", "24h" = "#9b59b6"),
    UP_DOWN = c(Up = "#e74c3c", Down = "#3498db", NS = "#bdc3c7"),
    HEATMAP = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
)

VOLCANO <- list(
    WIDTH = GLOBAL_PLOT_WIDTH, # Standard figure width for all volcanos
    HEIGHT = GLOBAL_PLOT_HEIGHT, # Standard figure height for all volcanos
    XLIM = NULL, # NULL = auto-symmetric based on data
    YLIM = NULL, # NULL = auto-scale, or set c(0, 50)
    SHOW_LABELS = TRUE, # FALSE = don't show marker labels, TRUE = show them
    POINT_SIZE_NS = 1.0, # Non-significant points
    POINT_SIZE_SIG = 1.5, # Significant points
    POINT_ALPHA_NS = 0.4, # Transparency for NS
    POINT_ALPHA_SIG = 0.8, # Transparency for significant
    LABEL_SIZE = 3, # Gene label text size
    LABEL_COLOR = "black", # Color of the text and label segment
    LABEL_FONTFACE = "italic", # Font face of the text (e.g. "italic", "bold", "plain")
    LABEL_SEGMENT_COLOR = "gray30", # Color of the line pointing to the dot
    LABEL_SEGMENT_SIZE = 0.5, # Thickness of the pointing line
    LABEL_BOX_PADDING = 0.5, # Space around the text box
    THRESHOLD_LINE_COLOR = "gray40",
    THRESHOLD_LINE_TYPE = "dashed"
)

# --- Venn Diagram Settings ---
VENN <- list(
    PLOT_WIDTH       = GLOBAL_PLOT_WIDTH, # Plot width (inches)
    PLOT_HEIGHT      = GLOBAL_PLOT_HEIGHT, # Plot height (inches)
    LABEL_SIZE       = 6, # Count label size (legacy/ggVennDiagram)
    SET_SIZE         = 5, # Category name size (legacy/ggVennDiagram)
    EDGE_SIZE        = 1.2, # Circle outline thickness
    EDGE_LTY         = "solid", # Circle outline type
    TITLE_SIZE       = 16, # Plot title font size
    SUBTITLE_SIZE    = 12, # Subtitle font size
    FILL_LOW         = "white", # Gradient fill low end
    FILL_HIGH        = "gray70", # Gradient fill high end
    FONT_COLOR       = "black", # Label text color
    MARGIN           = c(30, 40, 30, 40), # Plot margins (t, r, b, l)
    # --- ggforce circle Venn config ---
    CIRCLE_RADIUS    = 1.5, # Circle radius (normalized units)
    CIRCLE_DIST      = 1.8, # Distance between circle centers (controls overlap)
    CIRCLE_ALPHA     = 0.3, # Circle fill transparency (0=transparent, 1=opaque)
    CAT_POSITION     = "top", # Category label position: "top" or "bottom"
    CAT_SIZE         = 5.5, # Category label font size (ggplot text size)
    COUNT_SIZE       = 7 # Count number font size inside regions
)

# --- ggplot2 Theme ---
if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)

    PUBLICATION_THEME <- theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
        theme(
            # Background & Grid
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
            panel.grid.minor = element_blank(),

            # Text & Titles
            plot.title = element_text(size = PLOT_FONT_SIZE_TITLE, face = "bold", hjust = 0.5, color = PLOT_FONT_COLOR, margin = margin(b = 10)),
            plot.subtitle = element_text(size = PLOT_FONT_SIZE_LABEL, hjust = 0.5, color = PLOT_FONT_COLOR, margin = margin(b = 10)),
            plot.margin = margin(20, 20, 20, 20), # Prevent clipping
            axis.title = element_text(face = "bold"),
            axis.text = element_text(color = "black"),

            # Legends
            legend.title = element_text(face = "bold"),
            legend.position = "right",

            # Facets
            strip.background = element_rect(fill = "#F0F0F0", color = NA),
            strip.text = element_text(face = "bold")
        )

    theme_set(PUBLICATION_THEME)
}


# ==============================================================================
# 5. CONFIGURABLE PLOT TITLES
# ==============================================================================
# Centralized control for all figure titles and subtitles.
# Access via: PLOT_TITLES$fig_id$main or PLOT_TITLES$fig_id$sub

PLOT_TITLES <- list(
    # --- Figure 1: QC ---
    fig1_pca = list(
        main = NULL,
        sub = NULL
    ),
    fig1_correlation = list(
        main = NULL,
        sub = NULL
    ),

    # --- Figure 2: Treatment effect against baseline 0h (Tm 6/12/24 vs YPD 0) ---
    fig2_volcano_Tm6h_vs_YPD0h = list(
        main = "WT Tm (6h) vs WT YPD (0h)"
    ),
    fig2_volcano_Tm12h_vs_YPD0h = list(
        main = "WT Tm (12h) vs WT YPD (0h)"
    ),
    fig2_volcano_Tm24h_vs_YPD0h = list(
        main = "WT Tm (24h) vs WT YPD (0h)"
    ),

    # --- Figure 3: Timecourse Volcanoes ---
    fig3a_volcano_6h_vs_0h = list(
        main = "WT Tm (6h) vs WT Tm (0h)"
    ),
    fig3b_volcano_12h_vs_6h = list(
        main = "WT Tm (12h) vs WT Tm (6h)"
    ),
    fig3c_volcano_24h_vs_12h = list(
        main = "WT Tm (24h) vs WT Tm (12h)"
    ),

    # --- Figure 4: UPR Comparisons ---
    fig4a_volcano_wt_6h = list(
        main = "WT Tm (6h) vs WT YPD (6h)"
    ),
    fig4b_volcano_ko_6h = list(
        main = bquote(bolditalic("hac1")*bolditalic(Delta)~bold("Tm vs")~bolditalic("hac1")*bolditalic(Delta)~bold("YPD (6h)"))
    ),
    fig4c_venn = list(
        main = "Venn Diagram: Overlap of DEGs",
        sub = paste0("Comparison of Gene Sets (Mode: ", VENN_MODE, ")")
    ),
    fig4d_venn_direction = list(
        main = "Directional Venn Diagram: Overlap of DEGs",
        sub = paste0("Comparison of Gene Sets (Mode: ", VENN_MODE, ")")
    )
)


# ==============================================================================
# 6. GENE LISTS & ANNOTATIONS
# ==============================================================================

# --- Marker Genes Definition ---
# Lists for specific figure panels
MARKERS <- list(
    VALIDATION = c(
        "YFL031W", "YJL034W", "YOL013C", "YCL043C", # CLUSTER1: HAC1, KAR2, HRD1, PDI1
        "YDL229W", "YPL131W", # CLUSTER2: SSB1, RPL5
        "YHR030C", "YKR097W", "YBR023C", # CLUSTER3: SLT2, PCK1, CHS3
        "YLR185W", "YOL127W" # CLUSTER4: RPL37A, RPL25
    )
)

# Combine for easy access
ALL_MARKERS <- unlist(MARKERS)

# --- Gene Name Mapping (ORF -> Common Name) ---
GENE_NAMES <- c(
    # CLUSTER 1
    YFL031W = "HAC1", YJL034W = "KAR2", YOL013C = "HRD1", YCL043C = "PDI1",
    # CLUSTER 2
    YDL229W = "SSB1", YPL131W = "RPL5",
    # CLUSTER 3
    YHR030C = "SLT2", YKR097W = "PCK1", YBR023C = "CHS3",
    # CLUSTER 4
    YLR185W = "RPL37A", YOL127W = "RPL25"
)

# Helper function to get gene name
get_gene_name <- function(orf) {
    if (orf %in% names(GENE_NAMES)) {
        return(GENE_NAMES[orf])
    }
    return(orf)
}

# Save figure to PNG and/or PDF
save_figure <- function(plot, filename, width = 10, height = 8, dpi = 300, dir = FIGURE_DIR) {
    if (PLOT_FILE_TYPE %in% c("png", "both")) {
        png_path <- file.path(dir, paste0(filename, ".png"))
        ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "transparent")
        message("  Saved PNG: ", png_path)
    }
    if (PLOT_FILE_TYPE %in% c("pdf", "both")) {
        pdf_path <- file.path(dir, paste0(filename, ".pdf"))
        ggsave(pdf_path, plot, width = width, height = height, bg = "transparent")
        message("  Saved PDF: ", pdf_path)
    }
}

# ==============================================================================
# CONFIGURATION LOADED SUCCESSFULLY
# ==============================================================================
message("Loaded 00_config.R:")
message("  - Sample Exclusion: ", USE_EXCLUDED_SAMPLES)
message("  - Cutoffs: padj < ", PADJ_CUTOFF, ", |log2FC| > ", LFC_CUTOFF)
message("  - Venn Mode: ", VENN_MODE)
message("  - Plot Titles: Loaded")

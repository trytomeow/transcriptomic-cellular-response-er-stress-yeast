# ==============================================================================
# Figure 4A: Volcano - WT Tm (6h) vs WT YPD (6h)
# ==============================================================================
# Script: fig4a_volcano_wt_6h.R
# Output: fig4a_volcano_wt_6h.png/pdf
# Description: UPR volcano for WT genotype at 6h (Tm vs YPD).
# ==============================================================================

# --- LOCAL CONFIG ---
SHOW_LABELS <- FALSE
SHOW_LEGEND <- TRUE

# --- STEP 1: Load Shared Components ---
shared_path <- "fig4_volcano_shared.R"
if (!file.exists(shared_path)) shared_path <- "../fig4_deg_upr/fig4_volcano_shared.R"
source(shared_path)

MARKER_GENES <- MARKERS$VALIDATION

message("Running: fig4a_volcano_wt_6h")

# --- STEP 2: Load Pre-Computed Results ---
res_file <- file.path(DESEQ_DIR, "deg_WT_6h_TmVsYPD.csv")
if (!file.exists(res_file)) {
    stop("CRITICAL ERROR: ", res_file, " not found. Run fig4_volcano_shared.R first.")
}
res <- read.csv(res_file)

# --- STEP 3: Generate Plot ---
p <- create_labeled_volcano(
    res         = res,
    title       = PLOT_TITLES$fig4a_volcano_wt_6h$main,
    label_genes = MARKER_GENES,
    show_labels = SHOW_LABELS,
    show_legend = SHOW_LEGEND
)

# --- STEP 4: Save Figure ---
save_figure(p, "fig4a_volcano_wt_6h", width = VOLCANO$WIDTH, height = VOLCANO$HEIGHT, dir = FIG4_DIR)

message("Done: fig4a_volcano_wt_6h")

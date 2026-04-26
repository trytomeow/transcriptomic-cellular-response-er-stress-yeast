# ==============================================================================
# Figure 1A: PCA - PC1 vs PC2
# ==============================================================================
# Script: fig1a_pca_pc1_pc2.R
# Output: fig1a_pca_pc1_pc2.png/pdf
# Description: Individual PCA plot for principal components 1 and 2.
# ==============================================================================

# --- STEP 1: Load Shared PCA Components ---
source("fig1_pca_shared.R")

message("Running: fig1a_pca_pc1_pc2")

# --- STEP 2: Generate Plot ---
p_pc1_pc2 <- create_pca_plot(
    data = pca_df, 
    pc_x = 1, 
    pc_y = 2, 
    var_exp = var_explained,
    title = PLOT_TITLES$fig1_pca$main,
    subtitle = PLOT_TITLES$fig1_pca$sub
)

# --- STEP 3: Save Figure ---
save_figure(p_pc1_pc2, "fig1a_pca_pc1_pc2", width = 8, height = 6, dir = FIG1_DIR)

message("Done: fig1a_pca_pc1_pc2")

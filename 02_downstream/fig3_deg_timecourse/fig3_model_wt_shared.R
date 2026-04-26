# ==============================================================================
# Figure 3: Model Verification (Core Model)
# ==============================================================================
# Script: fig3_model_wt_shared.R
# Output: (none — verification only)
# Description: Checks that core_dds.rds exists before running fig3 scripts.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

message("\n[INFO] --- VERIFYING CORE MODEL ---")

dds_path <- file.path(DESEQ_DIR, "core_dds.rds")

if (!file.exists(dds_path)) {
    stop("CRITICAL ERROR: core_dds.rds not found. Please run fig0_config/01_core_model.R first.")
}

message("[SUCCESS] Core Model verified. Ready for fig3 downstream scripts.\n")

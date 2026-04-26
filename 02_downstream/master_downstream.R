# ==============================================================================
# SCRIPT: master_downstream.R
# DESCRIPTION:
#   Master execution script for the Transcriptomics Senior Project Analysis Pipeline.
#   This script sequentially executes all analysis modules, ensuring dependencies
#   (configuration, shared functions) are loaded correctly.
#
# USAGE:
#   Run this script from RStudio or via command line: Rscript master_downstream.R
#
# DEPENDENCIES:
#   - 00_config.R (Must be located in 'fig0_config' or project root)
# ==============================================================================

rm(list = ls())
start_time <- Sys.time()

message("==============================================================================")
message("              TRANSCRIPTOMICS ANALYSIS PIPELINE: EXECUTION STARTED            ")
message("==============================================================================")
message("Start Time: ", start_time)
message("")

# ==============================================================================
# 1. PIPELINE CONFIGURATION
# ==============================================================================

# Search for 00_config.R relative to common working directories
POSSIBLE_CONFIG_PATHS <- c(
    "fig0_config/00_config.R",             # run from 02_downstream/
    "02_downstream/fig0_config/00_config.R" # run from repo root
)

config_loaded <- FALSE

for (path in POSSIBLE_CONFIG_PATHS) {
    if (file.exists(path)) {
        message("[INFO] Configuration file found at: ", path)
        source(path)
        config_loaded <- TRUE
        break
    }
}

if (!config_loaded) {
    stop(
        "\n[CRITICAL ERROR] Configuration file '00_config.R' not found.\n",
        "  Run this script from within the '02_downstream/' directory:\n",
        "    setwd('path/to/repo/02_downstream')\n",
        "    source('master_downstream.R')\n",
        "  Checked paths:\n",
        paste(POSSIBLE_CONFIG_PATHS, collapse = "\n  ")
    )
}

# ==============================================================================
# 2. ENVIRONMENT INITIALIZATION
# ==============================================================================

if (exists("SCRIPTS_DIR") && dir.exists(SCRIPTS_DIR)) {
    message("[INFO] Scripts Directory Verified: ", SCRIPTS_DIR)
} else {
    warning("[WARNING] SCRIPTS_DIR is not defined or invalid. Scripts may fail to load.")
}

message("[INFO] Current Working Directory: ", getwd())

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# Execute an R script with error handling and directory switching
# @param script_rel_path Relative path to the R script (from SCRIPTS_DIR).
# @return Status string: "Success", "Skipped", or "Failed".
execute_script <- function(script_rel_path) {
    script_full_path <- file.path(SCRIPTS_DIR, script_rel_path)

    if (!file.exists(script_full_path)) {
        message(sprintf("  [SKIP] File not found: %s", script_full_path))
        return("Skipped (Missing)")
    }

    script_name <- basename(script_full_path)
    script_dir <- dirname(script_full_path)

    message(sprintf("\n[EXEC] Running: %s...", script_name))

    # Store current working directory to restore later
    old_wd <- getwd()
    on.exit(setwd(old_wd))

    tryCatch(
        {
            setwd(script_dir)
            source(script_name, local = FALSE)
            message(sprintf("  [DONE] Completed: %s", script_name))
            return("Success")
        },
        error = function(e) {
            message(sprintf("  [ERROR] Execution failed for %s", script_name))
            message(sprintf("  [REASON] %s", e$message))
            return(paste("Failed:", e$message))
        }
    )
}

# ==============================================================================
# 4. SCRIPT EXECUTION QUEUE
# ==============================================================================

# Paths are RELATIVE to SCRIPTS_DIR.
PIPELINE_SCRIPTS <- c(
    # --- MODULE 0: CORE DATA MODEL ---
    "fig0_config/01_core_model.R",

    # --- MODULE 1: QC & PCA ---
    "fig1_qc/fig1_pca_shared.R",
    "fig1_qc/fig1a_pca_pc1_pc2.R",
    "fig1_qc/fig1b_correlation.R",

    # --- MODULE 2: Treatment 0h (Baseline Model) ---
    "fig2_deg_treatment_0h/fig2_model_wt_shared.R",
    "fig2_deg_treatment_0h/fig2a_volcano_Tm6h_vs_YPD0h.R",
    "fig2_deg_treatment_0h/fig2b_volcano_Tm12h_vs_YPD0h.R",
    "fig2_deg_treatment_0h/fig2c_volcano_Tm24h_vs_YPD0h.R",
    "fig2_deg_treatment_0h/fig2d_deg_count.R",
    "fig2_deg_treatment_0h/fig2e_go_treatment.R",
    "fig2_deg_treatment_0h/fig2f_kegg_treatment.R",

    # --- MODULE 3: Timecourse (Temporal Model) ---
    "fig3_deg_timecourse/fig3_model_wt_shared.R",
    "fig3_deg_timecourse/fig3a_volcano_6h_vs_0h.R",
    "fig3_deg_timecourse/fig3b_volcano_12h_vs_6h.R",
    "fig3_deg_timecourse/fig3c_volcano_24h_vs_12h.R",
    "fig3_deg_timecourse/fig3d_deg_count.R",
    "fig3_deg_timecourse/fig3e_go_temporal.R",
    "fig3_deg_timecourse/fig3f_kegg_temporal.R",
    "fig3_deg_timecourse/fig3g_mfuzz_elbow.R",
    "fig3_deg_timecourse/fig3h_mfuzz_clustering.R",
    "fig3_deg_timecourse/fig3i_mfuzz_go.R",
    "fig3_deg_timecourse/fig3j_mfuzz_kegg.R",

    # --- MODULE 4: UPR (Genotype Interaction) ---
    "fig4_deg_upr/fig4_volcano_shared.R",
    "fig4_deg_upr/fig4_upr_shared.R",
    "fig4_deg_upr/fig4a_volcano_wt_6h.R",
    "fig4_deg_upr/fig4b_volcano_ko_6h.R",
    "fig4_deg_upr/fig4c_venn.R",
    "fig4_deg_upr/fig4d_venn_direction.R",
    "fig4_deg_upr/fig4e_go_direction.R",
    "fig4_deg_upr/fig4f_kegg_upr_direction.R"
)


# ==============================================================================
# 5. EXECUTION & SUMMARY
# ==============================================================================

execution_log <- character(length(PIPELINE_SCRIPTS))
names(execution_log) <- PIPELINE_SCRIPTS

message("\n[INFO] Starting execution of ", length(PIPELINE_SCRIPTS), " scripts...")

for (i in seq_along(PIPELINE_SCRIPTS)) {
    script <- PIPELINE_SCRIPTS[i]
    execution_log[i] <- execute_script(script)
}

# --- Completion Summary ---
end_time <- Sys.time()
elapsed_min <- round(difftime(end_time, start_time, units = "mins"), 2)

message("\n==============================================================================")
message("              TRANSCRIPTOMICS ANALYSIS PIPELINE: COMPLETED                    ")
message("==============================================================================")
message("Total Runtime: ", elapsed_min, " minutes")
message("")
message("--- Execution Status Report ---")
print(data.frame(
    Script = names(execution_log),
    Status = execution_log,
    row.names = NULL
))
message("==============================================================================")

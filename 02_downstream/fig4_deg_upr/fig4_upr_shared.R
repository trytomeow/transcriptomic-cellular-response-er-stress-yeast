# ==============================================================================
# Figure 4: UPR Comparison - Shared Components
# ==============================================================================
# Script: fig4_upr_shared.R
# Output: (none — provides shared objects for fig4 downstream scripts)
# Description: Shared DESeq2 analysis for UPR-dependent vs Independent
#              classification. Identifies WT and KO DEGs at 6h.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
})

message("Running: fig4_upr_shared")

# --- STEP 2: Load Core Model (for enrichment universe) ---
dds_path <- file.path(DESEQ_DIR, "core_dds.rds")
if (!file.exists(dds_path)) {
    stop("CRITICAL ERROR: core_dds.rds not found. Run 01_core_model.R first.")
}
dds <- readRDS(dds_path)
message("  Loaded Core DESeq2 Model (for gene universe).")

# --- STEP 3: Load Pre-Computed DEG Results (from fig4_volcano_shared.R) ---
wt_csv <- file.path(DESEQ_DIR, "deg_WT_6h_TmVsYPD.csv")
ko_csv <- file.path(DESEQ_DIR, "deg_KO_6h_TmVsYPD.csv")
if (!file.exists(wt_csv) || !file.exists(ko_csv)) {
    stop("CRITICAL ERROR: DEG CSVs not found. Run fig4_volcano_shared.R (sourced by fig4a/fig4b) first.")
}
res_wt <- read.csv(wt_csv)
res_ko <- read.csv(ko_csv)
message("  Loaded pre-computed contrast results from CSV.")

deg_wt <- res_wt$gene[!is.na(res_wt$padj) & res_wt$padj < PADJ_CUTOFF & abs(res_wt$log2FoldChange) > LFC_CUTOFF]
deg_ko <- res_ko$gene[!is.na(res_ko$padj) & res_ko$padj < PADJ_CUTOFF & abs(res_ko$log2FoldChange) > LFC_CUTOFF]

message("DEGs Identified (|LFC| > ", LFC_CUTOFF, "):")
message("  - WT (6h): ", length(deg_wt))
message("  - KO (6h): ", length(deg_ko))

# --- STEP 4: Categorize Genes (UPR Dependency) ---
upr_dependent  <- setdiff(deg_wt, deg_ko)
upr_independent <- intersect(deg_wt, deg_ko)
ko_specific    <- setdiff(deg_ko, deg_wt)

message("Gene Categories:")
message("  - UPR-dependent (WT only): ", length(upr_dependent))
message("  - UPR-independent (Both): ", length(upr_independent))
message("  - KO-specific (KO only): ", length(ko_specific))

# --- Shared Helper: Split Gene List by Direction ---
get_direction_subsets <- function(gene_list, group_name, lfc_df) {
    df <- lfc_df %>% filter(ORF %in% gene_list)
    up_genes   <- df %>% filter(log2FoldChange >  LFC_CUTOFF) %>% pull(ORF)
    down_genes <- df %>% filter(log2FoldChange < -LFC_CUTOFF) %>% pull(ORF)
    list(
        list(name = "All",          group = group_name, genes = gene_list),
        list(name = "Upregulated",  group = group_name, genes = up_genes),
        list(name = "Downregulated",group = group_name, genes = down_genes)
    )
}

message("UPR comparison shared components loaded successfully.")

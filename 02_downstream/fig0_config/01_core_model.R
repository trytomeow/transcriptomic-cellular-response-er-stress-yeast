# ==============================================================================
# Config: Core DESeq2 Model Builder
# ==============================================================================
# Script: 01_core_model.R
# Output: core_dds.rds, core_vsd.rds, core_meta.rds
# Description: Builds the master DESeq2 model (~group) for all downstream
#              comparisons. Ensures identical dispersion and size factors.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(DESeq2)
    library(dplyr)
})

message("\n[INFO] --- BUILDING CORE MODEL (All Samples) ---")

# --- STEP 2: Load & Filter Data ---
counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
metadata <- read.csv(METADATA_FILE, row.names = 1)
colnames(metadata) <- tolower(colnames(metadata))

# Filter OUT bad samples if configured
if (exists("USE_EXCLUDED_SAMPLES") && USE_EXCLUDED_SAMPLES) {
    samples_to_keep <- !(rownames(metadata) %in% SAMPLES_TO_EXCLUDE)
    metadata <- metadata[samples_to_keep, ]
}

counts <- counts[, rownames(metadata)]
message("  Included Samples: ", ncol(counts))

# --- STEP 3: Setup Design Group ---
metadata$genotype  <- factor(metadata$genotype,  levels = c("WT", "KO"))
metadata$treatment <- factor(metadata$treatment, levels = c("YPD", "Tm"))
metadata$time      <- factor(metadata$time,      levels = c(0, 6, 12, 24))
metadata$group     <- factor(paste0(metadata$genotype, "_", metadata$treatment, "_", metadata$time))

message("  Groups defined: ", paste(levels(metadata$group), collapse = ", "))

# --- STEP 4: Build & Run DESeq2 ---
dds_core <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = metadata,
    design    = ~group
)

dds_core <- dds_core[rowSums(counts(dds_core)) >= 10, ]
message("  Genes passing pre-filter: ", nrow(dds_core))

dds_core <- DESeq(dds_core, quiet = TRUE)
vsd_core <- vst(dds_core, blind = FALSE)

# --- STEP 5: Save Centralized Data ---
saveRDS(dds_core, file.path(DESEQ_DIR, "core_dds.rds"))
saveRDS(vsd_core, file.path(DESEQ_DIR, "core_vsd.rds"))
saveRDS(metadata, file.path(DESEQ_DIR, "core_meta.rds"))

message("[SUCCESS] Core Model Built and Saved to ", DESEQ_DIR)
message("          Downstream scripts MUST use readRDS() to load this data.\n")

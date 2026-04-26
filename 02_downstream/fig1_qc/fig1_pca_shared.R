# ==============================================================================
# Figure 1: PCA - Shared Components
# ==============================================================================
# Script: fig1_pca_shared.R
# Output: (none — sourced by fig1a/b/c and fig1_pca_combined.R)
# Description: Data loading, VST normalization, PCA computation, and
#              shared plotting function / color palette for all PCA plots.
# ==============================================================================

config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

# --- STEP 1: Load Required Packages ---
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

message("Running: fig1_pca_shared")

# --- STEP 2: Load Raw Data ---
counts   <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
metadata <- read.csv(METADATA_FILE, row.names = 1)

colnames(metadata) <- tolower(colnames(metadata))
counts <- counts[, rownames(metadata)]

message("Loaded Data:")
message("  - Samples: ", ncol(counts))
message("  - Genes: ", nrow(counts))

# --- STEP 3: Sample Exclusion (Configurable) ---
if (exists("USE_EXCLUDED_SAMPLES") && USE_EXCLUDED_SAMPLES) {
    message("\n[FILTER] Applying Sample Exclusion List...")
    samples_to_keep <- !(rownames(metadata) %in% SAMPLES_TO_EXCLUDE)

    # Update objects
    metadata <- metadata[samples_to_keep, ]
    counts <- counts[, rownames(metadata)]

    message("  - Excluded: ", paste(SAMPLES_TO_EXCLUDE, collapse = ", "))
    message("  - Remaining Samples: ", ncol(counts))
} else {
    message("\n[FILTER] Using ALL samples (Exclusion Disabled).")
}

metadata$genotype  <- factor(metadata$genotype, levels = c("WT", "KO"))
metadata$treatment <- factor(metadata$treatment, levels = c("YPD", "Tm"))
metadata$time      <- factor(metadata$time)
metadata$group     <- factor(paste(metadata$genotype, metadata$treatment, metadata$time, sep = "_"))

# Initialize DESeq object (Design ~group is standard for versatile contrasts)
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~group
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

# blind=TRUE: unsupervised QC; prevents design-driven bias in variance stabilisation.
message("Applying VST transformation...")
vsd <- vst(dds, blind = TRUE)
message("  - Transformation complete. ", nrow(vsd), " genes retained.")

pca_data <- prcomp(t(assay(vsd)), scale. = TRUE)
pca_df   <- as.data.frame(pca_data$x)

pca_df$sample    <- rownames(pca_df)
pca_df$genotype  <- metadata[rownames(pca_df), "genotype"]
pca_df$treatment <- metadata[rownames(pca_df), "treatment"]
pca_df$time      <- metadata[rownames(pca_df), "time"]

pca_df$group <- ifelse(
    pca_df$time == 0,
    paste0(pca_df$genotype, "-0h"),
    paste(pca_df$genotype, pca_df$treatment, paste0(pca_df$time, "h"), sep = "-")
)

var_explained <- round(100 * summary(pca_data)$importance[2, ], 1)
message(
    "Variance Explained:",
    " PC1=", var_explained[1], "%",
    " PC2=", var_explained[2], "%",
    " PC3=", var_explained[3], "%"
)

PCA_COLORS <- c(
    "WT-0h" = "#aec7e8", "WT-YPD-6h" = "#98df8a", "WT-Tm-6h" = "#17becf",
    "WT-YPD-12h" = "#ffbb78", "WT-Tm-12h" = "#d62728",
    "WT-YPD-24h" = "#c5b0d5", "WT-Tm-24h" = "#ff9896",
    # KO samples
    "KO-0h" = "#bcbd22", "KO-YPD-6h" = "#9467bd", "KO-Tm-6h" = "#1f77b4"
)

# Standardized PCA Plot Function
create_pca_plot <- function(data, pc_x, pc_y, var_exp, title = NULL, subtitle = NULL, show_labels = TRUE) {
    x_var <- paste0("PC", pc_x)
    y_var <- paste0("PC", pc_y)

    p <- ggplot(data, aes(
        x = .data[[x_var]], y = .data[[y_var]],
        color = group
    )) +
        geom_point(size = 4, alpha = 0.85)

    if (show_labels) {
        p <- p + geom_text_repel(
            aes(label = sample),
            size = 2, color = "#000000",
            max.overlaps = 20, box.padding = 0.4,
            family = PLOT_FONT_FAMILY
        )
    }

    p <- p +
        scale_color_manual(values = PCA_COLORS, name = "Condition") +
        labs(
            x        = sprintf("%s (%s%%)", x_var, var_exp[pc_x]),
            y        = sprintf("%s (%s%%)", y_var, var_exp[pc_y]),
            title    = title,
            subtitle = subtitle
        ) +
        theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
        theme(
            plot.background  = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.title       = element_text(size = PLOT_FONT_SIZE_TITLE, face = "bold", hjust = 0.5),
            plot.subtitle    = element_text(hjust = 0.5, color = "#000000"),
            legend.text      = element_text(size = 9),
            panel.border     = element_rect(color = "#000000", fill = NA)
        )

    return(p)
}

message("PCA shared components loaded successfully.")

# ==============================================================================
# Figure 1B: Pearson Correlation Matrix
# ==============================================================================
# Script: fig1b_correlation.R
# Output: qc_correlation_heatmap.png, qc_correlation_barplot.png
# Description: Computes sample-to-sample Pearson correlation on VST data
#              and generates a heatmap and per-sample replicate correlation barplot.
# ==============================================================================
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

# --- LOCAL CONFIG ---
CORRELATION_CUTOFF <- 0.9  # Samples below this mean replicate correlation are flagged as outliers

# --- Load Required Packages ---
suppressPackageStartupMessages({
    library(DESeq2)
    library(pheatmap)
    library(ggplot2)
    library(RColorBrewer)
})

message("Running: fig1b_correlation")

# --- Load Data ---
counts <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
metadata <- read.csv(METADATA_FILE, row.names = 1)

# Standardize column names to lowercase
colnames(metadata) <- tolower(colnames(metadata))

# Ensure metadata and counts match
common_samples <- intersect(rownames(metadata), colnames(counts))
metadata <- metadata[common_samples, ]
counts <- counts[, common_samples]

metadata$group <- factor(paste(metadata$genotype, metadata$treatment, metadata$time, sep = "_"))

message("Loaded ", ncol(counts), " samples and ", nrow(counts), " genes")

# --- Sample Exclusion (Configurable) ---
if (exists("USE_EXCLUDED_SAMPLES") && USE_EXCLUDED_SAMPLES) {
    message("\n[FILTER] Applying Sample Exclusion List...")
    samples_to_keep <- !(rownames(metadata) %in% SAMPLES_TO_EXCLUDE)
    metadata <- metadata[samples_to_keep, ]
    counts   <- counts[, rownames(metadata)]
    message("  - Excluded: ", paste(SAMPLES_TO_EXCLUDE, collapse = ", "))
    message("  - Remaining Samples: ", ncol(counts))
} else {
    message("\n[FILTER] Using ALL samples (Exclusion Disabled).")
}

# --- Create DESeq2 Object ---
dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = metadata,
    design = ~1
)

# Filter low count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# blind=TRUE: unsupervised QC; prevents design-driven bias in variance stabilisation.
message("Applying VST transformation...")
vsd <- vst(dds, blind = TRUE)

# --- Compute Correlation ---
message("Computing Pearson correlation matrix...")
vst_mat <- assay(vsd)
cor_mat <- cor(vst_mat, method = "pearson")

# --- Plotting ---
output_file_png <- file.path(FIG1_DIR, paste0("fig1b_correlation_heatmap", ".png"))
output_file_pdf <- file.path(FIG1_DIR, paste0("fig1b_correlation_heatmap", ".pdf"))

# Define annotation for heatmap
annotation_col <- data.frame(
    Genotype = metadata$genotype,
    Treatment = metadata$treatment,
    Time = factor(metadata$time)
)
rownames(annotation_col) <- rownames(metadata)

# Define colors
ann_colors <- list(
    Genotype = c(WT = "grey70", KO = "black"),
    Treatment = c(YPD = "lightblue", Tm = "salmon"),
    Time = c("0" = "#EFF3FF", "6" = "#BDD7E7", "12" = "#6BAED6", "24" = "#2171B5")
)

# Define color palette (Blue -> Yellow -> Red) for correlation
cor_mat_plot <- cor_mat
cor_mat_plot[cor_mat_plot < 0.8] <- 0.8
cor_breaks <- seq(0.8, 1, length.out = 101)
cor_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# pheatmap is incompatible with ggsave(); saved via png()/pdf() directly.
draw_correlation_heatmap <- function() {
    pheatmap::pheatmap(cor_mat_plot,
        color = cor_colors,
        breaks = cor_breaks,
        legend_breaks = seq(0.8, 1, by = 0.05),
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        show_rownames = TRUE,
        show_colnames = TRUE,
        treeheight_row = 25,
        treeheight_col = 25,
        fontsize = 10,
        cellwidth = 18,
        cellheight = 18,
        border_color = "grey50"
    )
}

png(output_file_png, width = 11, height = 11, units = "in", res = 300, bg = "transparent")
draw_correlation_heatmap()
dev.off()
message("Correlation heatmap saved to: ", output_file_png)

pdf(output_file_pdf, width = 11, height = 11, bg = "transparent")
draw_correlation_heatmap()
dev.off()
message("Correlation heatmap saved to: ", output_file_pdf)

# --- Barplot of Mean Correlation with Replicates ---
message("Generating correlation barplot (vs replicates)...")

try(graphics.off(), silent = TRUE)

tryCatch(
    {
        # Calculate mean correlation for each sample with its own group members (replicates)
        samples <- rownames(metadata)
        mean_replicate_cors <- numeric(length(samples))
        names(mean_replicate_cors) <- samples

        for (s in samples) {
            # Find replicates (same group, but not self)
            s_group <- metadata[s, "group"]

            # Ensure group is valid
            if (is.na(s_group)) {
                mean_replicate_cors[s] <- NA
                next
            }

            # Identify replicates
            replicates <- rownames(metadata)[metadata$group == s_group & rownames(metadata) != s]

            # Check if any replicates found and they exist in cor_mat
            replicates <- intersect(replicates, colnames(cor_mat))

            if (length(replicates) > 0) {
                # Calculate mean correlation with replicates
                vals <- cor_mat[s, replicates]
                mean_replicate_cors[s] <- mean(vals, na.rm = TRUE)
            } else {
                # If no replicates, NA
                mean_replicate_cors[s] <- NA
            }
        }

        # Prepare Dataframe
        barplot_df <- data.frame(
            Sample = names(mean_replicate_cors),
            MeanCorrelation = mean_replicate_cors,
            Group = metadata[names(mean_replicate_cors), "group"]
        )

        # Remove NAs if any (samples without replicates)
        barplot_df <- barplot_df[!is.na(barplot_df$MeanCorrelation), ]

        # Sorting
        barplot_df <- barplot_df[order(barplot_df$Group, barplot_df$Sample), ]
        barplot_df$Sample <- factor(barplot_df$Sample, levels = barplot_df$Sample)

        p_bar <- ggplot(barplot_df, aes(x = Sample, y = MeanCorrelation, fill = Group)) +
            geom_bar(stat = "identity", color = "black", alpha = 0.9, width = 0.7) +
            geom_hline(aes(yintercept = CORRELATION_CUTOFF, linetype = "Cutoff"), color = "red", linewidth = 1) +
            scale_linetype_manual(
                name = "Guide", values = c("Cutoff" = "dashed"),
                labels = paste0("Cutoff < ", CORRELATION_CUTOFF)
            ) +
            scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
            labs(
                title = paste0("Mean Correlation with Group Replicates (",
                    if (exists("USE_EXCLUDED_SAMPLES") && USE_EXCLUDED_SAMPLES) "Filtered" else "All Samples",
                    ")"),
                subtitle = paste0("Samples below red line (< ", CORRELATION_CUTOFF, ") are potential outliers"),
                x = "Sample",
                y = "Mean Pearson Correlation (r)"
            ) +
            theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
            theme(
                plot.background     = element_rect(fill = "transparent", color = NA),
                panel.background    = element_rect(fill = "transparent", color = NA),
                plot.title          = element_text(face = "bold", hjust = 0.5),
                plot.subtitle       = element_text(color = "red", hjust = 0.5),
                axis.text.x         = element_text(angle = 45, hjust = 1, size = 10),
                panel.grid.major.x  = element_blank(),
                panel.grid.minor    = element_blank(),
                legend.position     = "right"
            )

        save_figure(p_bar, "fig1b_correlation_barplot", width = 12, height = 8, dir = FIG1_DIR)
    },
    error = function(e) {
        message("ERROR generating barplot: ", e$message)
    }
)

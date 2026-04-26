# ==============================================================================
# Figure 3F: KEGG Enrichment - Timecourse (Directional)
# ==============================================================================
# Script: fig3f_kegg_temporal.R
# Output: fig3f_kegg_{6h_vs_0h,12h_vs_6h,24h_vs_12h}.png/pdf,
#         fig3f_kegg_temporal_combined.csv
# Description: KEGG pathway enrichment (directional) for each timecourse transition.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(ggplot2)
    library(dplyr)
    library(stringr)
})

message("Running: fig3f_kegg_temporal")

csv_files <- list(
    "6h vs 0h"   = file.path(DESEQ_DIR, "fig3a_volcano_6h_vs_0h_data.csv"),
    "12h vs 6h"  = file.path(DESEQ_DIR, "fig3b_volcano_12h_vs_6h_data.csv"),
    "24h vs 12h" = file.path(DESEQ_DIR, "fig3c_volcano_24h_vs_12h_data.csv")
)

# --- STEP 2: Read and Filter DEGs ---
read_and_filter <- function(path, name) {
    if (!file.exists(path)) {
        warning(paste("File not found:", path))
        return(NULL)
    }
    df <- read.csv(path)
    if (!"gene" %in% colnames(df)) df$gene <- df[, 1]
    sig_df <- df %>%
        filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_CUTOFF) %>%
        dplyr::select(gene, log2FoldChange, padj)
    return(sig_df)
}

# --- STEP 3: Process Datasets ---
temp_data <- list()
for (n in names(csv_files)) {
    message(paste("Reading:", n))
    temp_data[[n]] <- read_and_filter(csv_files[[n]], n)
}

# --- STEP 4: Split into Directional Subsets ---
get_subsets <- function(df, group_name) {
    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }

    all_genes <- df$gene
    up_genes <- df %>%
        filter(log2FoldChange > LFC_CUTOFF) %>%
        pull(gene)
    down_genes <- df %>%
        filter(log2FoldChange < -LFC_CUTOFF) %>%
        pull(gene)

    list(
        list(name = "All", group = group_name, genes = all_genes),
        list(name = "Upregulated", group = group_name, genes = up_genes),
        list(name = "Downregulated", group = group_name, genes = down_genes)
    )
}

all_subsets <- c()
for (n in names(temp_data)) {
    all_subsets <- c(all_subsets, get_subsets(temp_data[[n]], n))
}

message(paste("Total subsets:", length(all_subsets)))


# --- STEP 5: Run KEGG Enrichment ---
# Note: enrichKEGG with organism='sce' uses KEGG's internal background; universe is ignored.
run_kegg <- function(subset_info) {
    genes <- subset_info$genes
    if (length(genes) < 5) {
        return(NULL)
    }

    tryCatch(
        {
            kegg_res <- enrichKEGG(
                gene          = genes,
                organism      = "sce",
                keyType       = "kegg",
                pAdjustMethod = "BH",
                pvalueCutoff  = KEGG_PVALUE_CUTOFF,
                qvalueCutoff  = KEGG_QVALUE_CUTOFF
            )

            if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
                return(NULL)
            }

            res_df <- kegg_res@result %>%
                filter(p.adjust < KEGG_PVALUE_CUTOFF) %>%
                arrange(p.adjust) %>%
                head(TOP_N_KEGG$fig3f_temporal) %>%
                mutate(
                    Group = subset_info$group,
                    Direction = subset_info$name
                )

            if (nrow(res_df) == 0) {
                return(NULL)
            }
            return(res_df)
        },
        error = function(e) {
            message(paste("Error in", subset_info$group, ":", e$message))
            return(NULL)
        }
    )
}

message("Running KEGG Enrichment...")
results_list <- lapply(all_subsets, run_kegg)
combined_res <- bind_rows(results_list)

if (nrow(combined_res) == 0) {
    stop("No significant KEGG pathways found in any timecourse comparison. This is expected if DEGs are too few.")
}

# --- STEP 6: Prepare for Plotting ---
combined_res <- combined_res %>%
    mutate(
        ratio_parts = strsplit(GeneRatio, "/"),
        GeneRatioNum = sapply(ratio_parts, function(x) as.numeric(x[1]) / as.numeric(x[2])),
        Description = str_wrap(Description, width = 40),
        Direction = factor(Direction, levels = c("All", "Upregulated", "Downregulated"))
    )

# --- STEP 7: Plotting Function ---
create_dotplot <- function(data, title, filename, plot_height = 8) {

    p <- ggplot(data, aes(x = Direction, y = reorder(Description, GeneRatioNum))) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        facet_wrap(~Group, scales = "free_y", ncol = 1) +
        scale_x_discrete(expand = expansion(add = GO_PLOT$X_EXPAND)) +
        scale_color_gradient(low = "#3498db", high = "#e74c3c", name = "-log10(p.adjust)") +
        scale_size_continuous(name = "Gene Ratio", range = GO_PLOT$DOT_SIZE_RANGE) +
        labs(
            title = paste0("KEGG: ", title),
            x = NULL, y = NULL
        ) +
        guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
        theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
        theme(
            plot.background  = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.title       = element_text(size = PLOT_FONT_SIZE_TITLE, face = GO_PLOT$FONT_FACE_TITLE, hjust = 0.5),
            strip.text       = element_text(size = 10, face = GO_PLOT$FONT_FACE_STRIP, color = "black"),
            strip.background = element_rect(fill = "gray95", color = NA),
            axis.text.x      = element_text(color = "black", face = GO_PLOT$FONT_FACE_AXIS_X),
            axis.text.y      = element_text(size = 7.5, color = "black", lineheight = 1, face = GO_PLOT$FONT_FACE_AXIS_Y),
            panel.border     = element_rect(color = "black", fill = NA),
            panel.grid.minor = element_blank(),
            legend.position  = "right",
            legend.title     = element_text(size = 8, face = "bold"),
            legend.text      = element_text(size = 7),
            legend.key.size  = unit(0.4, "cm"),
            plot.margin      = margin(20, 20, 20, 20)
        )

    save_figure(p, filename, width = 6.15, height = plot_height, dir = FIG3_DIR)
    message(paste("Saved:", filename))
}

# --- STEP 8: Generate Individual Plots ---
data_p1 <- combined_res %>% filter(Group == "6h vs 0h")
if (nrow(data_p1) > 0) create_dotplot(data_p1, "WT Tm (6h) vs WT Tm (0h)", "fig3f_kegg_6h_vs_0h", plot_height = 8)

data_p2 <- combined_res %>% filter(Group == "12h vs 6h")
if (nrow(data_p2) > 0) create_dotplot(data_p2, "WT Tm (12h) vs WT Tm (6h)", "fig3f_kegg_12h_vs_6h", plot_height = 8)

data_p3 <- combined_res %>% filter(Group == "24h vs 12h")
if (nrow(data_p3) > 0) create_dotplot(data_p3, "WT Tm (24h) vs WT Tm (12h)", "fig3f_kegg_24h_vs_12h", plot_height = 8)

# Save Summary CSV
write.csv(combined_res %>% dplyr::select(-ratio_parts),
    file.path(FIG3_DIR, "fig3f_kegg_temporal_combined.csv"),
    row.names = FALSE
)

message("Summary CSV saved: fig3f_kegg_temporal_combined.csv")
message("Done: fig3f_kegg_temporal")

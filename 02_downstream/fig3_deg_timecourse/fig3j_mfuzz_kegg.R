# ==============================================================================
# Figure 3J: Mfuzz Cluster KEGG Enrichment
# ==============================================================================
# Script: fig3j_mfuzz_kegg.R
# Output: fig3j_mfuzz_kegg_dotplot.png/pdf, fig3j_mfuzz_kegg_heatmap.png/pdf,
#         fig3j_mfuzz_kegg.csv
# Description: KEGG pathway enrichment for each Mfuzz cluster. Requires fig3h first.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

# --- LOCAL CONFIG ---
TOP_N_PER_CLUSTER <- TOP_N_KEGG$fig3j_mfuzz
TOP_N_SUMMARY_HM  <- TOP_N_KEGG$fig3j_mfuzz_summary

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(grid)
})

message("Running: fig3j_mfuzz_kegg")

# --- STEP 2: Load Cluster Data ---
cluster_file <- file.path(FIG3_DIR, "fig3h_mfuzz_genes.csv")
if (!file.exists(cluster_file)) {
    stop("Run fig3h_mfuzz_clustering.R first! File not found: ", cluster_file)
}

cluster_data <- read.csv(cluster_file)
n_clusters <- length(unique(cluster_data$cluster))
message("Loaded ", nrow(cluster_data), " genes in ", n_clusters, " clusters.")

# --- STEP 3: KEGG Enrichment per Cluster ---
# Note: enrichKEGG with organism='sce' uses KEGG's internal background; universe is ignored.
run_kegg_cluster <- function(genes, cluster_id) {
    if (length(genes) < 5) {
        message("  Cluster ", cluster_id, ": skipped (< 5 genes)")
        return(NULL)
    }
    tryCatch(
        {
            kegg_res <- enrichKEGG(
                gene          = genes,
                organism      = "sce",
                pAdjustMethod = "BH",
                pvalueCutoff  = KEGG_PVALUE_CUTOFF,
                qvalueCutoff  = KEGG_QVALUE_CUTOFF
            )
            if (is.null(kegg_res) || nrow(kegg_res@result) == 0) {
                message("  Cluster ", cluster_id, ": no significant terms")
                return(NULL)
            }
            res_df <- as.data.frame(kegg_res) %>%
                filter(p.adjust < KEGG_PVALUE_CUTOFF) %>%
                arrange(p.adjust) %>%
                head(TOP_N_PER_CLUSTER)
            if (nrow(res_df) == 0) {
                message("  Cluster ", cluster_id, ": 0 terms")
                return(NULL)
            }
            res_df <- res_df %>%
                mutate(
                    Cluster = paste0("Cluster ", cluster_id),
                    ClusterNum = cluster_id,
                    GeneRatioNum = as.numeric(sapply(
                        strsplit(GeneRatio, "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2])
                    ))
                )
            message("  Cluster ", cluster_id, ": ", nrow(res_df), " terms")
            return(res_df)
        },
        error = function(e) {
            message("  Cluster ", cluster_id, " error: ", e$message)
            return(NULL)
        }
    )
}

message("Running KEGG enrichment per cluster...")
results_list <- lapply(sort(unique(cluster_data$cluster)), function(clid) {
    genes <- cluster_data$gene[cluster_data$cluster == clid]
    run_kegg_cluster(genes, clid)
})
names(results_list) <- paste0("Cluster_", sort(unique(cluster_data$cluster)))
results_list <- results_list[!sapply(results_list, is.null)]
combined_res <- bind_rows(results_list)

if (nrow(combined_res) == 0) {
    stop("No significant KEGG pathways found in any Mfuzz cluster. Run fig3h first and check cluster sizes.")
}

write.csv(combined_res, file.path(FIG3_DIR, "fig3j_mfuzz_kegg.csv"), row.names = FALSE)
message("Saved: fig3j_mfuzz_kegg.csv")

# --- STEP 4: Dot Plot ---
dot_data <- combined_res %>%
    mutate(Description = str_wrap(Description, width = 45)) %>%
    group_by(Cluster) %>%
    arrange(GeneRatioNum) %>%
    ungroup() %>%
    mutate(Description = factor(Description, levels = unique(Description)))

p_dot <- ggplot(dot_data, aes(
    x    = reorder(Cluster, ClusterNum),
    y    = Description
)) +
    geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.85) +
    scale_color_gradient2(
        low      = "#3498DB",
        mid      = "#F39C12",
        high     = "#E74C3C",
        midpoint = 2,
        name     = "-log10(p.adj)"
    ) +
    scale_size_continuous(name = "Gene Ratio", range = GO_PLOT$DOT_SIZE_RANGE) +
    scale_x_discrete(expand = expansion(add = GO_PLOT$X_EXPAND)) +
    labs(
        title = "KEGG Pathway Enrichment per Mfuzz Cluster",
        subtitle = paste0("Top ", TOP_N_PER_CLUSTER, " pathways per cluster | padj < ", KEGG_PVALUE_CUTOFF),
        x = NULL, y = NULL
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, size = 9),
        axis.text.x      = element_text(face = "bold"),
        axis.text.y      = element_text(lineheight = 0.85),
        panel.border     = element_rect(color = "black", fill = NA),
        plot.margin      = margin(20, 20, 20, 20),
        legend.position  = "right"
    )

dot_h <- length(unique(dot_data$Description)) * 0.25 + 3
dot_h <- max(dot_h, 8)

save_figure(p_dot, "fig3j_mfuzz_kegg_dotplot",
    width = 12, height = dot_h, dir = FIG3_DIR
)
message("Saved: fig3j_mfuzz_kegg_dotplot.png/.pdf")

# --- STEP 5: Summary Heatmap ---
message("Building summary KEGG heatmap...")

summary_wide <- combined_res %>%
    dplyr::select(Description, ClusterNum, p.adjust) %>%
    mutate(neg_log_padj = -log10(p.adjust)) %>%
    dplyr::select(-p.adjust) %>%
    pivot_wider(
        names_from = ClusterNum,
        values_from = neg_log_padj,
        values_fill = 0,
        names_prefix = "Cluster_"
    ) %>%
    as.data.frame()

rownames(summary_wide) <- summary_wide$Description
summary_wide$Description <- NULL

row_max <- apply(summary_wide, 1, max)
top_paths <- names(sort(row_max, decreasing = TRUE))[1:min(TOP_N_SUMMARY_HM, nrow(summary_wide))]
summary_mat <- as.matrix(summary_wide[top_paths, ])

summary_mat_scaled <- t(scale(t(summary_mat)))
summary_mat_scaled[is.na(summary_mat_scaled)] <- 0
colnames(summary_mat_scaled) <- gsub("Cluster_", "Cluster", colnames(summary_mat_scaled))

color_fun <- colorRamp2(c(-2, 0, 2), c("#3182bd", "white", "#de2d26"))

ht_summary <- Heatmap(
    summary_mat_scaled,
    name = "Z-score",
    col = color_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = unit(14, "cm"),
    cluster_columns = FALSE,
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot = -90,
    column_title = "KEGG: Temporal Clusters",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        legend_height = unit(5, "cm")
    ),
    width = unit(n_clusters * 1.5, "cm"),
    height = unit(nrow(summary_mat) * 0.45, "cm")
)

summary_h <- nrow(summary_mat) * 0.25 + 3
summary_h <- max(summary_h, 8)

png(file.path(FIG3_DIR, "fig3j_mfuzz_kegg_heatmap.png"),
    width = 12, height = summary_h, units = "in", res = 300
)
draw(ht_summary, merge_legend = TRUE)
dev.off()

pdf(file.path(FIG3_DIR, "fig3j_mfuzz_kegg_heatmap.pdf"),
    width = 12, height = summary_h
)
draw(ht_summary, merge_legend = TRUE)
dev.off()
message("Saved: fig3j_mfuzz_kegg_heatmap.png/.pdf")

message("Done: fig3j_mfuzz_kegg")

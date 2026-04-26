# ==============================================================================
# Figure 3I: Mfuzz Cluster GO Enrichment
# ==============================================================================
# Script: fig3i_mfuzz_go.R
# Output: fig3i_mfuzz_go_dotplot.png/pdf, fig3i_mfuzz_go_heatmap.png/pdf,
#         fig3i_mfuzz_go.csv
# Description: GO BP enrichment for each Mfuzz cluster. Requires fig3h first.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

# --- LOCAL CONFIG ---
TOP_N         <- TOP_N_GO$fig3i_mfuzz
TOP_N_SUMMARY <- TOP_N_GO$fig3i_mfuzz_summary

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Sc.sgd.db)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(grid)
})

message("Running: fig3i_mfuzz_go")

# --- STEP 2: Load Universe ---
dds_path <- file.path(DESEQ_DIR, "core_dds.rds")
if (file.exists(dds_path)) {
    universe_genes <- rownames(readRDS(dds_path))
    message("Using Core Model expressed genes as universe (N=", length(universe_genes), ")")
} else {
    counts_raw     <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
    universe_genes <- rownames(counts_raw)
    warning("core_dds.rds not found — falling back to counts matrix universe")
}

# --- STEP 3: Load Cluster Data ---
cluster_file <- file.path(FIG3_DIR, "fig3h_mfuzz_genes.csv")
if (!file.exists(cluster_file)) {
    stop("Run fig3h_mfuzz_clustering.R first! File not found: ", cluster_file)
}

cluster_data <- read.csv(cluster_file)
n_clusters   <- length(unique(cluster_data$cluster))
message("Loaded ", nrow(cluster_data), " genes in ", n_clusters, " clusters.")



# --- STEP 4: GO Enrichment per Cluster ---
run_go_cluster <- function(genes, cluster_id, universe) {
    if (length(genes) < 5) {
        message("  Cluster ", cluster_id, ": skipped (< 5 genes)")
        return(NULL)
    }
    tryCatch(
        {
            go_res <- enrichGO(
                gene          = genes,
                universe      = universe,
                OrgDb         = org.Sc.sgd.db,
                keyType       = "ORF",
                ont           = GO_ONTOLOGY,
                pAdjustMethod = "BH",
                pvalueCutoff  = GO_PVALUE_CUTOFF,
                qvalueCutoff  = GO_QVALUE_CUTOFF
            )
            if (is.null(go_res) || nrow(go_res@result) == 0) {
                message("  Cluster ", cluster_id, ": no significant terms")
                return(NULL)
            }

            go_res <- simplify(go_res, cutoff = GO_SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)

            if (is.null(go_res) || nrow(go_res@result) == 0) {
                message("  Cluster ", cluster_id, ": no terms left after simplification")
                return(NULL)
            }

            res_df <- as.data.frame(go_res) %>%
                filter(p.adjust < GO_PVALUE_CUTOFF) %>%
                filter(Count >= 3) %>%
                arrange(p.adjust) %>%
                head(TOP_N)

            if (nrow(res_df) == 0) {
                message("  Cluster ", cluster_id, ": 0 terms passed internal filters")
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

message("Running GO enrichment per cluster...")
results_list <- lapply(sort(unique(cluster_data$cluster)), function(clid) {
    genes <- cluster_data$gene[cluster_data$cluster == clid]
    run_go_cluster(genes, clid, universe_genes)
})
names(results_list) <- paste0("Cluster_", sort(unique(cluster_data$cluster)))

# Remove NULL entries (clusters with no significant GO terms)
results_list <- Filter(Negate(is.null), results_list)
combined_res <- bind_rows(results_list)

if (nrow(combined_res) == 0) {
    stop("No significant GO terms found in any Mfuzz cluster. Run fig3h first and check cluster sizes.")
}

write.csv(combined_res, file.path(FIG3_DIR, "fig3i_mfuzz_go.csv"), row.names = FALSE)
message("Saved: fig3i_mfuzz_go.csv")


# --- STEP 5: Dot Plot ---
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
        title = paste0("GO ", GO_ONTOLOGY, " Enrichment per Temporal Cluster"),
        subtitle = paste0("Up to ", TOP_N, " simplified terms per cluster | padj < ", GO_PVALUE_CUTOFF, " | Count \u2265 3"),
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

save_figure(p_dot, "fig3i_mfuzz_go_dotplot",
    width = max(GO_PLOT$WIDTH, 12), height = dot_h, dir = FIG3_DIR
)
message("Saved: fig3i_mfuzz_go_dotplot.png/.pdf")

# --- STEP 6: Summary Heatmap ---
message("Building summary heatmap...")

summary_wide <- combined_res %>%
    dplyr::select(Description, ClusterNum, p.adjust) %>%
    mutate(neg_log_padj = -log10(p.adjust)) %>%
    dplyr::select(-p.adjust) %>%
    pivot_wider(
        names_from  = ClusterNum,
        values_from = neg_log_padj,
        values_fill = 0,
        names_prefix = "Cluster_"
    ) %>%
    as.data.frame()

rownames(summary_wide) <- summary_wide$Description
summary_wide$Description <- NULL

row_max   <- apply(summary_wide, 1, max)
top_paths <- names(sort(row_max, decreasing = TRUE))[1:min(TOP_N_SUMMARY, nrow(summary_wide))]
summary_mat <- as.matrix(summary_wide[top_paths, ])

summary_mat_scaled <- t(scale(t(summary_mat)))
summary_mat_scaled[is.na(summary_mat_scaled)] <- 0
colnames(summary_mat_scaled) <- gsub("Cluster_", "Cluster", colnames(summary_mat_scaled))

color_fun <- colorRamp2(c(-2, 0, 2), c("#3182bd", "white", "#de2d26"))

ht_summary <- Heatmap(
    summary_mat_scaled,
    name = "Z-score",
    col  = color_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_rows        = TRUE,
    show_row_dend       = FALSE,
    row_names_side      = "right",
    row_names_gp        = gpar(fontsize = 10),
    row_names_max_width = unit(14, "cm"),
    cluster_columns     = FALSE,
    column_names_gp     = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot    = -90,
    column_title        = "GO BP: Temporal Clusters",
    column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
        title_gp      = gpar(fontsize = 10, fontface = "bold"),
        labels_gp     = gpar(fontsize = 9),
        legend_height = unit(5, "cm")
    ),
    width  = unit(n_clusters * 1.5, "cm"),
    height = unit(nrow(summary_mat) * 0.45, "cm")
)

summary_h <- nrow(summary_mat) * 0.25 + 3
summary_h <- max(summary_h, 8)

png(file.path(FIG3_DIR, "fig3i_mfuzz_go_heatmap.png"),
    width = 12, height = summary_h, units = "in", res = 300
)
draw(ht_summary, merge_legend = TRUE)
dev.off()

pdf(file.path(FIG3_DIR, "fig3i_mfuzz_go_heatmap.pdf"),
    width = 12, height = summary_h
)
draw(ht_summary, merge_legend = TRUE)
dev.off()
message("Saved: fig3i_mfuzz_go_heatmap.png/.pdf")

message("Done: fig3i_mfuzz_go")

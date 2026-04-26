# ==============================================================================
# Figure 3H: Mfuzz Temporal Clustering
# ==============================================================================
# Script: fig3h_mfuzz_clustering.R
# Output: fig3h_mfuzz_clusters.png/pdf, fig3h_mfuzz_heatmap.png/pdf,
#         fig3h_mfuzz_genes.csv, fig3h_zscore_matrix.csv
# Description: Fuzzy C-means clustering of timecourse DEGs with line plots,
#              expression heatmap, and gene-cluster assignment export.
# ==============================================================================

# --- LOCAL CONFIG ---
N_CLUSTERS  <- 4
RANDOM_SEED <- 42

# --- Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(DESeq2)
    library(Mfuzz)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(grid)
    library(gridExtra)
})

message("Running: fig3h_mfuzz_clustering")

# --- STEP 1: Load Normalized Data (WT Timecourse) ---
vsd_path <- file.path(DESEQ_DIR, "core_vsd.rds")
if (!file.exists(vsd_path)) {
    stop("CRITICAL ERROR: core_vsd.rds not found. Run 01_core_model.R first.")
}
vsd <- readRDS(vsd_path)
vsd_mat <- assay(vsd)

# Load metadata for time grouping
meta_wt_path <- file.path(DESEQ_DIR, "core_meta.rds")
meta_wt <- readRDS(meta_wt_path)

message("  Loaded pre-calculated VST data. WT samples: ", ncol(vsd_mat), " | Genes: ", nrow(vsd_mat))

# --- STEP 2: Select DEGs for Clustering ---
deg_files <- list(
    "6h"  = file.path(DESEQ_DIR, "fig3a_volcano_6h_vs_0h_data.csv"),
    "12h" = file.path(DESEQ_DIR, "fig3b_volcano_12h_vs_6h_data.csv"),
    "24h" = file.path(DESEQ_DIR, "fig3c_volcano_24h_vs_12h_data.csv")
)

deg_genes <- character(0)
for (label in names(deg_files)) {
    fpath <- deg_files[[label]]
    if (!file.exists(fpath)) {
        warning("Missing: ", fpath)
        next
    }
    df <- read.csv(fpath)
    if (!"gene" %in% colnames(df)) df$gene <- df[, 1]
    if (!"padj" %in% colnames(df)) df$padj <- 1
    if (!"log2FoldChange" %in% colnames(df)) df$log2FoldChange <- 0
    sig <- df$gene[!is.na(df$padj) & df$padj < PADJ_CUTOFF &
        abs(df$log2FoldChange) > LFC_CUTOFF]
    message(
        "  ", label, ": ", length(sig), " DEGs (padj < ", PADJ_CUTOFF,
        ", |LFC| > ", LFC_CUTOFF, ")"
    )
    deg_genes <- union(deg_genes, sig)
}

available_genes <- intersect(deg_genes, rownames(vsd_mat))
message("Union DEGs: ", length(deg_genes), " | In VST: ", length(available_genes))
if (length(available_genes) < 10) stop("Too few DEGs. Run fig3a/b/c first.")

# --- STEP 3: Calculate Group Mean Expression ---
timepoints <- c("0", "6", "12", "24")
time_labels <- c("0h", "6h", "12h", "24h")

mean_expr <- matrix(
    NA,
    nrow = length(available_genes),
    ncol = length(timepoints),
    dimnames = list(available_genes, time_labels)
)

for (i in seq_along(timepoints)) {
    t <- timepoints[i]
    if (t == "0") {
        samps <- rownames(meta_wt)[meta_wt$time == "0" & meta_wt$treatment == "YPD" & meta_wt$genotype == "WT"]
    } else {
        samps <- rownames(meta_wt)[meta_wt$time == t & meta_wt$treatment == "Tm" & meta_wt$genotype == "WT"]
    }
    if (length(samps) > 0) {
        mean_expr[available_genes, i] <- rowMeans(vsd_mat[available_genes, samps, drop = FALSE])
    }
}

mean_expr <- mean_expr[complete.cases(mean_expr), ]
message("Genes with complete data: ", nrow(mean_expr))

eset <- new("ExpressionSet", exprs = mean_expr)
eset <- standardise(eset)

# --- STEP 4: Estimate Fuzzifier ---
m_opt <- mestimate(eset)
message("  Fuzzifier: m = ", round(m_opt, 3))

# --- STEP 5: Run Mfuzz Clustering ---
message("\nRunning Mfuzz with k=", N_CLUSTERS, ", m=", round(m_opt, 3), " ...")
set.seed(RANDOM_SEED)
cl <- mfuzz(eset, c = N_CLUSTERS, m = m_opt)

# Summary
for (i in 1:N_CLUSTERS) {
    n <- sum(cl$cluster == i)
    message("  Cluster ", i, ": ", n, " genes")
}

# Set DO_REORDER=FALSE to inspect raw cluster IDs, then remap.
DO_REORDER <- TRUE
NEW_ORDER_MAP <- c(
    "2" = 1,
    "3" = 2,
    "1" = 3,
    "4" = 4
)

if (DO_REORDER) {
    message("\nApplying manual cluster reordering...")

    actual_ids   <- as.character(sort(unique(cl$cluster)))
    map_keys     <- names(NEW_ORDER_MAP)
    missing_keys <- setdiff(map_keys, actual_ids)
    if (length(missing_keys) > 0) {
        stop(
            "DO_REORDER map contains invalid cluster IDs: ",
            paste(missing_keys, collapse = ", "),
            "\nActual cluster IDs from this run: ", paste(actual_ids, collapse = ", "),
            "\nRun once with DO_REORDER=FALSE to inspect actual cluster numbers, then update NEW_ORDER_MAP."
        )
    }

    old_clusters <- as.character(cl$cluster)
    new_clusters <- integer(length(old_clusters))
    for (old_idx in names(NEW_ORDER_MAP)) {
        new_idx <- NEW_ORDER_MAP[[old_idx]]
        new_clusters[old_clusters == old_idx] <- new_idx
    }
    cl$cluster <- new_clusters

    old_membership <- cl$membership
    new_membership <- matrix(0, nrow = nrow(old_membership), ncol = ncol(old_membership))
    rownames(new_membership) <- rownames(old_membership)
    colnames(new_membership) <- 1:N_CLUSTERS

    for (old_idx_str in names(NEW_ORDER_MAP)) {
        old_idx <- as.integer(old_idx_str)
        new_idx <- NEW_ORDER_MAP[[old_idx_str]]
        new_membership[, new_idx] <- old_membership[, old_idx]
    }
    cl$membership <- new_membership

    message("  Successfully reordered clusters.")
    for (i in 1:N_CLUSTERS) {
        n <- sum(cl$cluster == i)
        message("  New Cluster ", i, ": ", n, " genes")
    }
}

# --- STEP 6: Prepare Plot Data ---
zscore_df <- as.data.frame(exprs(eset)) %>%
    mutate(gene = rownames(exprs(eset))) %>%
    pivot_longer(cols = -gene, names_to = "time", values_to = "zscore") %>%
    mutate(time_num = as.numeric(gsub("h", "", time)))

cluster_df <- data.frame(
    gene       = rownames(exprs(eset)),
    cluster    = cl$cluster,
    membership = apply(cl$membership, 1, max)
)

plot_data <- zscore_df %>% left_join(cluster_df, by = "gene")

# Cluster means
cluster_means <- plot_data %>%
    group_by(cluster, time, time_num) %>%
    summarise(mean_zscore = mean(zscore), .groups = "drop")

# Membership-based alpha: strong members (>= MFUZZ_MIN_MEMBERSHIP) are opaque,
# weak/boundary genes are faded — consistent with standard Mfuzz visualization.
plot_data$line_alpha <- ifelse(plot_data$membership >= MFUZZ_MIN_MEMBERSHIP, 0.65, 0.15)

# Cluster labels with n
cluster_sizes <- cluster_df %>%
    count(cluster) %>%
    mutate(label = paste0("Cluster ", cluster, "\n(n=", n, ")"))
label_map <- setNames(cluster_sizes$label, cluster_sizes$cluster)
plot_data$cluster_label <- label_map[as.character(plot_data$cluster)]
cluster_means$cluster_label <- label_map[as.character(cluster_means$cluster)]

# --- STEP 7: Mfuzz Line Plot ---
mem_palette <- colorRampPalette(c("#4575B4", "#ABD9E9", "#FEE090", "#D73027"))(100)

p_lines <- ggplot(plot_data, aes(x = time_num, y = zscore)) +
    geom_line(aes(group = gene, color = membership, alpha = line_alpha), linewidth = 0.35) +
    scale_alpha_identity() +
    geom_line(
        data = cluster_means,
        aes(x = time_num, y = mean_zscore),
        color = "#D73027", linewidth = 1.8
    ) +
    geom_point(
        data = cluster_means,
        aes(x = time_num, y = mean_zscore),
        color = "#D73027", size = 3, stroke = 0.8
    ) +
    scale_color_gradientn(
        colors = mem_palette,
        limits = c(0, 1),
        name   = "Membership"
    ) +
    facet_wrap(~cluster_label, ncol = 2) +
    scale_x_continuous(
        breaks = c(0, 6, 12, 24),
        labels = c("0h", "6h", "12h", "24h")
    ) +
    labs(
        title = "Temporal Gene Clustering",
        subtitle = paste0(
            "Fuzzy c-means | k=", N_CLUSTERS,
            " | m=", round(m_opt, 2), " | n=", nrow(mean_expr), " DEGs"
        ),
        x = "Time (hours)", y = "Z-score"
    ) +
    theme_minimal(base_size = 11, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, size = 9),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        panel.border     = element_rect(color = "black", fill = NA),
        legend.position  = "right"
    )

save_figure(p_lines, "fig3h_mfuzz_clusters",
    width = 10, height = ceiling(N_CLUSTERS / 2) * 3.5 + 2, dir = FIG3_DIR
)
message("Saved: fig3h_mfuzz_clusters.png/.pdf")

# ==============================================================================
# --- STEP 8: Expression Heatmap (ordered by cluster) — ggplot2 facet_wrap layout
#         WITH marker gene annotations
# ==============================================================================
# Layout matches the cluster line plot: ncol=2, same dimensions

# Sort genes by cluster then by membership (descending)
gene_order <- cluster_df %>%
    arrange(cluster, desc(membership)) %>%
    pull(gene)

heatmap_mat <- exprs(eset)[gene_order, ]

hm_df <- as.data.frame(heatmap_mat) %>%
    mutate(gene = rownames(heatmap_mat)) %>%
    left_join(cluster_df, by = "gene") %>%
    pivot_longer(
        cols = all_of(time_labels),
        names_to  = "time",
        values_to = "zscore"
    )

hm_df$gene         <- factor(hm_df$gene, levels = rev(gene_order))
hm_df$time         <- factor(hm_df$time, levels = time_labels)
hm_df$cluster_label <- label_map[as.character(hm_df$cluster)]
hm_df$zscore_clamped <- pmin(pmax(hm_df$zscore, -2.5), 2.5)

# --- Marker Gene Annotation Data ---
# Find which VALIDATION markers are in the heatmap and which cluster they belong to
marker_orfs <- MARKERS$VALIDATION
marker_in_hm <- data.frame(
    gene = marker_orfs[marker_orfs %in% gene_order],
    stringsAsFactors = FALSE
)

if (nrow(marker_in_hm) > 0) {
    marker_in_hm$gene_name <- sapply(marker_in_hm$gene, get_gene_name)
    marker_in_hm <- marker_in_hm %>%
        left_join(cluster_df[, c("gene", "cluster")], by = "gene")
    marker_in_hm$cluster_label <- label_map[as.character(marker_in_hm$cluster)]

    # odd clusters (1,3,...) = left column; even (2,4,...) = right column
    marker_in_hm$label_side <- ifelse(marker_in_hm$cluster %% 2 == 1, "left", "right")

    n_times <- length(time_labels)
    marker_in_hm$x_pos      <- ifelse(marker_in_hm$label_side == "left", 0.1, n_times + 0.9)
    marker_in_hm$hjust_val  <- ifelse(marker_in_hm$label_side == "left", 1, 0)
    marker_in_hm$x_seg_start <- ifelse(marker_in_hm$label_side == "left", 0.5, n_times + 0.5)
    marker_in_hm$x_seg_end  <- ifelse(marker_in_hm$label_side == "left", 0.15, n_times + 0.85)

    # Convert gene to factor with same levels as hm_df
    marker_in_hm$gene <- factor(marker_in_hm$gene, levels = levels(hm_df$gene))

    message("Marker genes found in heatmap: ", paste(marker_in_hm$gene_name, collapse = ", "))
} else {
    message("No marker genes found in heatmap data.")
}

# --- Build Heatmap ---
p_heatmap <- ggplot(hm_df, aes(x = time, y = gene, fill = zscore_clamped)) +
    geom_tile() +
    scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = c(-2.5, 2.5),
        name     = "Z-score"
    ) +
    facet_wrap(~cluster_label, ncol = 2, scales = "free_y") +
    scale_x_discrete(expand = c(0, 0)) +
    labs(
        title = "Temporal Gene Expression Heatmap",
        subtitle = paste0(
            "Genes ordered by Mfuzz membership | k=", N_CLUSTERS,
            " | n=", nrow(mean_expr), " DEGs"
        ),
        x = "Time", y = NULL
    ) +
    theme_minimal(base_size = PLOT_FONT_SIZE_LABEL, base_family = PLOT_FONT_FAMILY) +
    theme(
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle    = element_text(hjust = 0.5, size = 9),
        strip.text       = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "#EEEEEE", color = NA),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        panel.grid       = element_blank(),
        panel.border     = element_rect(color = "black", fill = NA),
        legend.position  = "right"
    )

if (nrow(marker_in_hm) > 0) {
    p_heatmap <- p_heatmap +
        coord_cartesian(clip = "off") +
        theme(plot.margin = margin(10, 50, 10, 50))

    p_heatmap <- p_heatmap +
        geom_segment(
            data = marker_in_hm,
            aes(
                x = x_seg_start, xend = x_seg_end,
                y = gene, yend = gene
            ),
            inherit.aes = FALSE,
            color = "gray30", linewidth = 0.3
        ) +
        geom_text(
            data = marker_in_hm,
            aes(x = x_pos, y = gene, label = gene_name, hjust = hjust_val),
            inherit.aes = FALSE,
            size = 2.5, fontface = "italic", color = "black",
            family = PLOT_FONT_FAMILY
        )
}

# Use same dimensions as line plot
hm_height <- ceiling(N_CLUSTERS / 2) * 3.5 + 2

save_figure(p_heatmap, "fig3h_mfuzz_heatmap",
    width = 10, height = hm_height, dir = FIG3_DIR
)
message("Saved: fig3h_mfuzz_heatmap.png/.pdf")

# --- STEP 9: Export Data ---
output_df <- cluster_df
for (i in 1:N_CLUSTERS) {
    output_df[[paste0("mem_cluster_", i)]] <- cl$membership[output_df$gene, i]
}
write.csv(output_df, file.path(FIG3_DIR, "fig3h_mfuzz_genes.csv"), row.names = FALSE)

zscore_out <- as.data.frame(exprs(eset))
zscore_out$gene <- rownames(zscore_out)
zscore_out$cluster <- cl$cluster[zscore_out$gene]
write.csv(zscore_out, file.path(FIG3_DIR, "fig3h_zscore_matrix.csv"), row.names = FALSE)

message("Saved: fig3h_mfuzz_genes.csv, fig3h_zscore_matrix.csv")
message("Done: fig3h_mfuzz_clustering")

# ==============================================================================
# Figure 2E: GO Enrichment - Treatment Effect (Tm vs YPD 0h)
# ==============================================================================
# Script: fig2e_go_treatment.R
# Output: fig2e_go_Tm{6,12,24}h_vs_YPD0h.png/pdf,
#         fig2e_go_treatment_combined.csv
# Description: GO BP enrichment (directional) for each treatment timepoint.
# ==============================================================================

# --- STEP 1: Load Configuration ---
config_path <- "../fig0_config/00_config.R"
if (!file.exists(config_path)) config_path <- "fig0_config/00_config.R"
source(config_path)

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Sc.sgd.db)
    library(ggplot2)
    library(dplyr)
    library(stringr)
})

message("Running: fig2e_go_treatment")

csv_files <- list(
    "Tm 6h vs YPD 0h"  = file.path(DESEQ_DIR, "fig2a_volcano_Tm6h_vs_YPD0h_data.csv"),
    "Tm 12h vs YPD 0h" = file.path(DESEQ_DIR, "fig2b_volcano_Tm12h_vs_YPD0h_data.csv"),
    "Tm 24h vs YPD 0h" = file.path(DESEQ_DIR, "fig2c_volcano_Tm24h_vs_YPD0h_data.csv")
)

# --- STEP 2: Read and Filter DEGs ---
read_and_filter <- function(path, name) {
    if (!file.exists(path)) {
        warning(paste("File not found:", path))
        return(NULL)
    }
    df <- read.csv(path)

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

message(paste("Total subsets to analyze:", length(all_subsets)))

# --- STEP 5: Run GO Enrichment ---
dds_path <- file.path(DESEQ_DIR, "core_dds.rds")
if (file.exists(dds_path)) {
    dds            <- readRDS(dds_path)
    universe_genes <- rownames(dds)
    message("Using Core Model expressed genes as universe (N=", length(universe_genes), ")")
} else {
    counts_raw     <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
    universe_genes <- rownames(counts_raw)
    message("Using all genes in counts_matrix as universe.")
}

run_go <- function(subset_info) {
    genes <- subset_info$genes
    if (length(genes) < 5) return(NULL)

    tryCatch(
        {
            go_res <- enrichGO(
                gene = genes,
                universe = universe_genes,
                OrgDb = org.Sc.sgd.db,
                keyType = "ORF",
                ont = GO_ONTOLOGY,
                pAdjustMethod = "BH",
                pvalueCutoff = GO_PVALUE_CUTOFF,
                qvalueCutoff = GO_QVALUE_CUTOFF
            )

            if (is.null(go_res) || nrow(go_res@result) == 0) {
                return(NULL)
            }

            go_res <- simplify(go_res, cutoff = GO_SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
            if (is.null(go_res) || nrow(go_res@result) == 0) {
                return(NULL)
            }

            res_df <- go_res@result %>%
                filter(p.adjust < GO_PVALUE_CUTOFF) %>%
                filter(qvalue   < GO_QVALUE_CUTOFF) %>%
                arrange(p.adjust) %>%
                head(TOP_N_GO$fig2e_treatment) %>%
                mutate(
                    Group = subset_info$group,
                    Direction = subset_info$name
                )
            return(res_df)
        },
        error = function(e) {
            message(paste("Error in", subset_info$group, ":", e$message))
            return(NULL)
        }
    )
}

message("Running GO Enrichment (Top 10 BP)...")
results_list <- lapply(all_subsets, run_go)
combined_res <- bind_rows(results_list)

if (nrow(combined_res) == 0) {
    stop("No significant GO terms found in any comparison. This is expected if DEGs are too few.")
}

# --- STEP 6: Prepare for Plotting ---
plot_size_label <- switch(GO_PLOT_X_AXIS,
    "FoldEnrichment" = "Fold Enrichment",
    "Count"          = "Gene Count",
    "Gene Ratio"
)

combined_res <- combined_res %>%
    mutate(
        ratio_parts  = strsplit(as.character(GeneRatio), "/"),
        GeneRatioNum = sapply(ratio_parts, function(x) as.numeric(x[1]) / as.numeric(x[2])),
        bg_parts     = strsplit(as.character(BgRatio), "/"),
        BgRatioNum   = sapply(bg_parts, function(x) as.numeric(x[1]) / as.numeric(x[2])),
        FoldEnrichment = GeneRatioNum / BgRatioNum,
        PlotSizeVar  = switch(GO_PLOT_X_AXIS,
            "FoldEnrichment" = FoldEnrichment,
            "Count"          = as.numeric(Count),
            GeneRatioNum
        ),
        Description = str_wrap(Description, width = 40),
        Direction = factor(Direction, levels = c("All", "Upregulated", "Downregulated"))
    )

# --- STEP 7: Plotting Function ---
create_dotplot <- function(data, title, filename, plot_height = 8) {

    p <- ggplot(data, aes(x = Direction, y = reorder(Description, GeneRatioNum))) +
        geom_point(aes(size = PlotSizeVar, color = -log10(p.adjust)), alpha = 0.8) +
        facet_wrap(~Group, scales = "free_y", ncol = 1) +
        scale_x_discrete(expand = expansion(add = GO_PLOT$X_EXPAND)) +
        scale_color_gradient(low = "#3498db", high = "#e74c3c", name = "-log10(p.adjust)") +
        scale_size_continuous(name = plot_size_label, range = GO_PLOT$DOT_SIZE_RANGE) +
        labs(
            title = paste0("GO BP: ", title),
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

    save_figure(p, filename, width = 6.15, height = plot_height, dir = FIG2_DIR)
    message(paste("Saved:", filename))
}

# --- STEP 8: Generate Individual Plots ---
data_p1 <- combined_res %>% filter(Group == "Tm 6h vs YPD 0h")
if (nrow(data_p1) > 0) create_dotplot(data_p1, "WT Tm (6h) vs WT YPD (0h)", "fig2e_go_Tm6h_vs_YPD0h", plot_height = 8)

data_p2 <- combined_res %>% filter(Group == "Tm 12h vs YPD 0h")
if (nrow(data_p2) > 0) create_dotplot(data_p2, "WT Tm (12h) vs WT YPD (0h)", "fig2e_go_Tm12h_vs_YPD0h", plot_height = 8)

data_p3 <- combined_res %>% filter(Group == "Tm 24h vs YPD 0h")
if (nrow(data_p3) > 0) create_dotplot(data_p3, "WT Tm (24h) vs WT YPD (0h)", "fig2e_go_Tm24h_vs_YPD0h", plot_height = 8)

# Save Summary CSV
write.csv(combined_res %>% dplyr::select(-ratio_parts, -bg_parts),
    file.path(FIG2_DIR, "fig2e_go_treatment_combined.csv"),
    row.names = FALSE
)

message("Summary CSV saved: fig2e_go_treatment_combined.csv")
message("Done: fig2e_go_treatment")

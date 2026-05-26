# ==============================================================================
# Figure 4f: KEGG Enrichment - UPR-Dependent vs Independent
# ==============================================================================
# Script: fig4f_kegg_upr_direction.R
# Output: fig4f_kegg_upr_*.png/pdf, fig4f_kegg_upr_results.csv
# Description: KEGG pathway enrichment for UPR-Dependent vs Independent genes,
#              split by direction (All/Up/Down).
# ==============================================================================

# --- Load Shared UPR Analysis ---
shared_path <- "fig4_upr_shared.R"
if (!file.exists(shared_path)) shared_path <- "../fig4_deg_upr/fig4_upr_shared.R"
source(shared_path)

# --- LOCAL CONFIG ---
TOP_N_LIMIT <- TOP_N_KEGG$fig4f_upr

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(patchwork)
})

message("Running: fig4f_kegg_upr_direction")

# --- STEP 1: Prepare Direction Subsets ---
lfc_wt <- res_wt %>% dplyr::select(gene, log2FoldChange) %>% dplyr::rename(ORF = gene)
lfc_ko <- res_ko %>% dplyr::select(gene, log2FoldChange) %>% dplyr::rename(ORF = gene)

lists_dependent   <- get_direction_subsets(upr_dependent,  "UPR-Dependent",   lfc_wt)
lists_independent <- get_direction_subsets(upr_independent, "UPR-Independent", lfc_wt)
lists_ko          <- get_direction_subsets(ko_specific,     "KO-Specific",     lfc_ko)
all_subsets <- c(lists_dependent, lists_independent, lists_ko)

message("Subset Sizes:")
for (s in all_subsets) {
    message(sprintf("  %s - %s: %d genes", s$group, s$name, length(s$genes)))
}

# --- STEP 2: Run KEGG Enrichment ---
# Note: enrichKEGG with organism='sce' uses KEGG's internal background; universe is ignored.
run_kegg <- function(subset_info) {
    genes <- subset_info$genes
    if (length(genes) < 5) {
        message(sprintf("  [%s] %s: skipped (< 5 genes)", subset_info$group, subset_info$name))
        return(NULL)
    }

    tryCatch(
        {
            res <- enrichKEGG(
                gene          = genes,
                organism      = "sce",
                keyType       = "kegg",
                pAdjustMethod = "BH",
                pvalueCutoff  = KEGG_PVALUE_CUTOFF,
                qvalueCutoff  = KEGG_QVALUE_CUTOFF
            )

            if (is.null(res) || nrow(res@result) == 0) {
                return(NULL)
            }

            df <- as.data.frame(res) %>%
                filter(p.adjust < KEGG_PVALUE_CUTOFF) %>%
                arrange(p.adjust) %>%
                head(TOP_N_LIMIT)

            if (nrow(df) == 0) {
                return(NULL)
            }

            df <- df %>%
                mutate(
                    Group = subset_info$group,
                    Direction = subset_info$name,
                    GeneRatioNum = as.numeric(sapply(
                        strsplit(as.character(GeneRatio), "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2])
                    )),
                    BgRatioNum = as.numeric(sapply(
                        strsplit(as.character(BgRatio), "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2])
                    )),
                    FoldEnrichment = GeneRatioNum / BgRatioNum,
                    PlotSizeVar = switch(GO_PLOT_X_AXIS,
                        "FoldEnrichment" = FoldEnrichment,
                        "Count"          = as.numeric(Count),
                        GeneRatioNum
                    )
                )
            df$GeneRatioNum <- as.numeric(df$GeneRatioNum)

            return(df)
        },
        error = function(e) {
            message(sprintf("  Error in [%s] %s: %s", subset_info$group, subset_info$name, e$message))
            return(NULL)
        }
    )
}

message("\nRunning KEGG enrichment...")
results_list <- lapply(all_subsets, run_kegg)
combined <- bind_rows(results_list)

if (nrow(combined) == 0) {
    stop("No significant KEGG pathways found for any UPR gene subset. This is expected if subsets are too small.")
}

combined$Direction   <- factor(combined$Direction,   levels = c("All", "Upregulated", "Downregulated"))
combined$Group       <- factor(combined$Group,       levels = c("UPR-Dependent", "UPR-Independent", "KO-Specific"))
combined$Description <- str_wrap(combined$Description, width = 40)

write.csv(combined, file.path(FIG4_DIR, "fig4f_kegg_upr_results_combined.csv"), row.names = FALSE)
message("Saved: fig4f_kegg_upr_results_combined.csv")

# --- STEP 3: Generate Plots ---
plot_size_label <- switch(GO_PLOT_X_AXIS,
    "FoldEnrichment" = "Fold Enrichment",
    "Count"          = "Gene Count",
    "Gene Ratio"
)

create_dotplot <- function(data, title, filename, plot_height = 8) {
    if (is.null(data) || nrow(data) == 0) {
        message(paste("  No sig terms for", title, "- Skipping plot."))
        return()
    }

    p <- ggplot(data, aes(x = Direction, y = reorder(Description, GeneRatioNum))) +
        geom_point(aes(size = PlotSizeVar, color = -log10(p.adjust)), alpha = 0.8) +
        facet_wrap(~Group, ncol = 1, scales = "free_y", labeller = as_labeller(c(
            "UPR-Dependent" = "bold('UPR-Dependent')",
            "UPR-Independent" = "bold('UPR-Independent')",
            "KO-Specific" = "bolditalic('hac1')*bolditalic(Delta)*bold('-Specific')"
        ), default = label_parsed)) +
        scale_color_gradient(low = "#3498db", high = "#e74c3c", name = "-log10(p.adjust)") +
        scale_x_discrete(expand = expansion(add = GO_PLOT$X_EXPAND)) +
        scale_size_continuous(name = plot_size_label, range = GO_PLOT$DOT_SIZE_RANGE) +
        labs(
            title = if (is.character(title)) paste0("KEGG: ", title) else bquote(bold("KEGG:") ~ .(title)),
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

    save_figure(p, filename, width = 6.15, height = plot_height, dir = FIG4_DIR)
    message("  Saved plot: ", filename)
}

create_dotplot(combined, "All Directional Pathways", "fig4f_kegg_upr_direction_combined", plot_height = 30)

dep_data <- combined %>% filter(Group == "UPR-Dependent")
create_dotplot(dep_data, "UPR-Dependent", "fig4f_kegg_upr_direction_dependent", plot_height = 8)

indep_data <- combined %>% filter(Group == "UPR-Independent")
create_dotplot(indep_data, "UPR-Independent", "fig4f_kegg_upr_direction_independent", plot_height = 8)

ko_data <- combined %>% filter(Group == "KO-Specific")
create_dotplot(ko_data, bquote(bolditalic("hac1")*bolditalic(Delta)*bold("-Specific")), "fig4f_kegg_upr_direction_ko_specific", plot_height = 8)

message("Done: fig4f_kegg_upr_direction")

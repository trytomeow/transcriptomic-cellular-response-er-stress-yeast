# ==============================================================================
# Figure 4e: GO Enrichment Directional (All/Up/Down)
# ==============================================================================
# Script: fig4e_go_direction.R
# Output: fig4e_go_direction_combined.png/pdf, fig4e_go_direction_combined.csv
# Description: GO BP enrichment comparing UPR-Dependent vs Independent genes
#              broken down by direction (All/Up/Down).
# ==============================================================================

# --- STEP 1: Load Shared Components ---
shared_path <- "fig4_upr_shared.R"
if (!file.exists(shared_path)) shared_path <- "../fig4_deg_upr/fig4_upr_shared.R"
source(shared_path)
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Sc.sgd.db)
    library(ggplot2)
    library(dplyr)
    library(stringr)
})

message("Running: fig4e_go_direction")

# --- STEP 2: Validate Inputs ---
if (!exists("upr_dependent") || !exists("upr_independent")) {
    stop("Error: Gene lists not found.")
}
if (!exists("res_wt")) {
    stop("Error: 'res_wt' (DESeq2 results) not found in fig4_upr_shared.R")
}

# --- STEP 3: Prepare Direction Subsets ---
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

# --- STEP 4: Run Enrichment ---
universe_genes <- rownames(dds)

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
                arrange(p.adjust) %>%
                filter(p.adjust < GO_PVALUE_CUTOFF) %>%
                filter(qvalue   < GO_QVALUE_CUTOFF) %>%
                head(TOP_N_GO$fig4e_direction) %>%
                mutate(
                    Group = subset_info$group,
                    Direction = subset_info$name
                )
            return(res_df)
        },
        error = function(e) {
            message(sprintf("    Error in %s-%s: %s", subset_info$group, subset_info$name, e$message))
            return(NULL)
        }
    )
}

message("\nRunning GO Enrichment (BP)...")
results_list <- lapply(all_subsets, run_go)
combined_res <- bind_rows(results_list)

if (nrow(combined_res) == 0) {
    stop("No significant GO terms found for any subset.")
}

# --- STEP 5: Process Data for Plotting ---
combined_res <- combined_res %>%
    mutate(
        ratio_parts  = strsplit(GeneRatio, "/"),
        GeneRatioNum = sapply(ratio_parts, function(x) as.numeric(x[1]) / as.numeric(x[2])),
        Description  = str_wrap(Description, width = 40)
    )
combined_res$Direction <- factor(combined_res$Direction, levels = c("All", "Upregulated", "Downregulated"))
combined_res$Group     <- factor(combined_res$Group,     levels = c("UPR-Dependent", "UPR-Independent", "KO-Specific"))

# --- STEP 6: Plotting Function ---
create_dotplot <- function(data, title_suffix = "", filename, plot_height = 16) {
    p <- ggplot(data, aes(
        x = Direction,
        y = reorder(Description, GeneRatioNum)
    )) +
        geom_point(aes(size = GeneRatioNum, color = -log10(p.adjust)), alpha = 0.8) +
        facet_wrap(~Group, scales = "free_y", ncol = 1, labeller = as_labeller(c(
            "UPR-Dependent" = "bold('UPR-Dependent')",
            "UPR-Independent" = "bold('UPR-Independent')",
            "KO-Specific" = "bolditalic('hac1')*bolditalic(Delta)*bold('-Specific')"
        ), default = label_parsed)) +
        scale_color_gradient(
            low = "#3498db", high = "#e74c3c",
            name = "-log10(p.adjust)"
        ) +
        scale_x_discrete(expand = expansion(add = GO_PLOT$X_EXPAND)) +
        scale_size_continuous(name = "Gene Ratio", range = GO_PLOT$DOT_SIZE_RANGE) +
        labs(
            title = if (is.character(title_suffix) && title_suffix == "") {
                "GO BP"
            } else if (is.character(title_suffix)) {
                paste0("GO BP: ", title_suffix)
            } else {
                bquote(bold("GO BP:") ~ .(title_suffix))
            },
            x = NULL,
            y = NULL
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

# --- STEP 7: Generate Plots ---
create_dotplot(combined_res, "", "fig4e_go_direction_combined", plot_height = 35)

dep_data <- combined_res %>% filter(Group == "UPR-Dependent")
create_dotplot(dep_data, "UPR-Dependent", "fig4e_go_direction_dependent", plot_height = 8)

indep_data <- combined_res %>% filter(Group == "UPR-Independent")
create_dotplot(indep_data, "UPR-Independent", "fig4e_go_direction_independent", plot_height = 8)

ko_data <- combined_res %>% filter(Group == "KO-Specific")
create_dotplot(ko_data, bquote(bolditalic("hac1")*bolditalic(Delta)*bold("-Specific")), "fig4e_go_direction_ko_specific", plot_height = 8)

# --- STEP 8: Save Results ---
export_df <- combined_res %>% dplyr::select(-ratio_parts)
write.csv(export_df, file.path(FIG4_DIR, "fig4e_go_direction_combined.csv"), row.names = FALSE)

message("Done: fig4e_go_direction")
message("  Saved 4 Plots: Combined, Dependent, Independent, KO-Specific")
message("  Saved CSV: fig4e_go_direction_combined.csv")

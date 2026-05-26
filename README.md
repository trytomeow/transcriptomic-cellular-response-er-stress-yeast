# ER Stress Transcriptomics Analysis Pipeline

## Overview

R pipeline for generating figures for transcriptomic analysis of ER stress response in *S. cerevisiae* (WT vs *hac1*Δ).

---

## Quick Start

```r
# 1. Clone or download this repository
# 2. Set working directory to the downstream folder
setwd("path/to/02_downstream")

# 3. Run the complete downstream pipeline
source("master_downstream.R")
```

> **Path auto-detection:** `00_config.R` automatically detects the project root from its
> own file location. No manual path configuration is needed as long as you `source()` the
> pipeline from within `02_downstream/`.
> If auto-detection fails, open `fig0_config/00_config.R` and set `BASE_DIR` manually:
>
> ```r
> BASE_DIR <- "/path/to/repo"
> ```

---

## Repository Structure

```
repo/
├── data/
│   ├── raw/                           ← Raw FASTQ files (input to upstream)
│   ├── gene_counts.csv                ← Processed count matrix (input to downstream)
│   └── metadata.csv                   ← Sample metadata (genotype, treatment, time)
│
├── refs/                              ← Reference genome files (not committed)
│   ├── hisat2_index/                    HISAT2 splice-aware index
│   ├── Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
│   └── Saccharomyces_cerevisiae.R64.gtf
│
├── 01_upstream/                       ← Bash pipeline (runs on Linux/HPC)
│   ├── 00_config.sh                     Shared config & helper functions
│   ├── master_upstream.sh               Orchestrates all upstream steps
│   ├── 01_fastqc_pretrim.sh
│   ├── 02_trimmomatic.sh
│   ├── 03_fastqc_posttrim.sh
│   ├── 04_hisat2_alignment.sh
│   ├── 05_variant_calling.sh
│   ├── 06_stringtie.sh
│   ├── 07_featurecounts.sh
│   └── 08_multiqc.sh
│
├── 02_downstream/                     ← R pipeline (runs locally in RStudio)
│   ├── master_downstream.R              Run this to execute the full pipeline
│   ├── fig0_config/
│   │   ├── 00_config.R                  Global settings (paths, colors, thresholds)
│   │   └── 01_core_model.R              Builds the shared DESeq2 model
│   ├── fig1_qc/                         PCA & sample correlation
│   ├── fig2_deg_treatment_0h/           Treatment effect vs baseline (0h)
│   ├── fig3_deg_timecourse/             Temporal DEG dynamics & Mfuzz clustering
│   └── fig4_deg_upr/                    UPR pathway analysis (WT vs hac1Δ)
│
└── results/                           ← Auto-created on first run
    ├── upstream/                        Upstream QC, BAM, VCF, StringTie outputs
    │   ├── 01_fastqc_pretrim/
    │   ├── 02_trimmed/
    │   ├── 03_fastqc_posttrim/
    │   ├── 04_sorted_bam/
    │   ├── 05_variants/
    │   ├── 06_stringtie/
    │   ├── 07_counts/                   Raw featureCounts output
    │   └── 08_multiqc/
    ├── fig1_qc/
    ├── fig2_deg_treatment_0h/
    ├── fig3_deg_timecourse/
    ├── fig4_deg_upr/
    └── deseq2_results/                  RDS objects & CSV tables
```

> **Upstream → Downstream handoff:**
> After running the upstream pipeline, process the featureCounts output and place it in `data/`:
> ```bash
> # Convert tab-delimited featureCounts output to CSV for R
> # results/upstream/07_counts/gene_counts.txt → data/gene_counts.csv
> ```
> The `metadata.csv` file must be created manually (sample names, genotype, treatment, timepoint).


---

## Script Folder Structure

| Phase                     | Folder                                   | Description                                           |
| ------------------------- | ---------------------------------------- | ----------------------------------------------------- |
| **Upstream** (Bash) | `01_upstream/`                         | FastQC, Trimmomatic, HISAT2, StringTie, featureCounts |
| **Downstream** (R)  | `02_downstream/fig0_config/`           | Configuration & Setup                                 |
|                           | `02_downstream/fig1_qc/`               | PCA & Quality Control                                 |
|                           | `02_downstream/fig2_deg_treatment_0h/` | Treatment vs Base 0h                                  |
|                           | `02_downstream/fig3_deg_timecourse/`   | Temporal DEG Dynamics                                 |
|                           | `02_downstream/fig4_deg_upr/`          | UPR DEG Analysis                                      |

---

## Dependencies

### Upstream (Bash Pipeline)

Tools must be available in your `$PATH` (recommended: install via conda).

| Tool                    | Purpose                                |
| ----------------------- | -------------------------------------- |
| FastQC                  | Raw & trimmed read quality control     |
| Trimmomatic             | Adapter trimming and quality filtering |
| HISAT2                  | Splice-aware RNA-seq alignment         |
| SAMtools                | BAM sorting, indexing, and processing  |
| BCFtools                | Variant calling and VCF filtering      |
| StringTie               | Transcript assembly and quantification |
| featureCounts (Subread) | Gene-level read count matrix           |
| MultiQC                 | Aggregate QC report                    |

### Downstream (R Pipeline)

**R** ≥ 4.3 and **Bioconductor** ≥ 3.18 are recommended.

#### Bioconductor Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",           # Differential expression analysis
  "Biobase",          # ExpressionSet (required by Mfuzz)
  "clusterProfiler",  # GO & KEGG enrichment analysis
  "org.Sc.sgd.db",    # S. cerevisiae gene annotation database
  "ComplexHeatmap",   # Expression heatmaps
  "Mfuzz"             # Fuzzy c-means temporal clustering
))
```

#### CRAN Packages

```r
install.packages(c(
  # Data manipulation
  "dplyr",
  "tidyr",
  "stringr",

  # Plotting
  "ggplot2",
  "ggrepel",      # Gene label repulsion on volcano plots
  "patchwork",    # Multi-panel figure composition
  "pheatmap",     # Correlation heatmaps
  "RColorBrewer", # Color palettes
  "ggforce",      # Circle-based Venn diagrams
  "circlize",     # Color mapping (required by ComplexHeatmap)
  "grid",
  "gridExtra"
))
```

---

## Figure 0: Configuration

| Script              | Description                            |
| ------------------- | -------------------------------------- |
| `00_config.R`     | Global settings, paths, colors, themes |
| `01_core_model.R` | Builds unified DESeq2 model            |

---

## Figure 1: Data Quality & PCA

| Script                  | Content                    |
| ----------------------- | -------------------------- |
| `fig1_pca_shared.R`   | Shared functions for PCA   |
| `fig1a_pca_pc1_pc2.R` | PCA: PC1 vs PC2            |
| `fig1b_correlation.R` | Sample correlation heatmap |

---

## Figure 2: DEG Treatment Baseline (0h)

**Statistical Model:** Centralized Core Model (`core_dds.rds`)

| Script                             | Comparison                 |
| ---------------------------------- | -------------------------- |
| `fig2_model_wt_shared.R`         | Verifies Core Model exists |
| `fig2a_volcano_Tm6h_vs_YPD0h.R`  | WT-Tm-6h vs WT-YPD-0h      |
| `fig2b_volcano_Tm12h_vs_YPD0h.R` | WT-Tm-12h vs WT-YPD-0h     |
| `fig2c_volcano_Tm24h_vs_YPD0h.R` | WT-Tm-24h vs WT-YPD-0h     |
| `fig2d_deg_count.R`              | DEG Count Barplot          |
| `fig2e_go_treatment.R`           | GO Enrichment              |
| `fig2f_kegg_treatment.R`         | KEGG Enrichment            |

---

## Figure 3: DEG Timecourse

**Statistical Model:** Centralized Core Model (`core_dds.rds`)

| Script                         | Comparison                 |
| ------------------------------ | -------------------------- |
| `fig3_model_wt_shared.R`     | Verifies Core Model exists |
| `fig3a_volcano_6h_vs_0h.R`   | WT-Tm-6h vs WT-0h          |
| `fig3b_volcano_12h_vs_6h.R`  | WT-Tm-12h vs 6h            |
| `fig3c_volcano_24h_vs_12h.R` | WT-Tm-24h vs 12h           |
| `fig3d_deg_count.R`          | Temporal DEG Count         |
| `fig3e_go_temporal.R`        | GO Enrichment (temporal)   |
| `fig3f_kegg_temporal.R`      | KEGG Enrichment (temporal) |
| `fig3g_mfuzz_elbow.R`        | Mfuzz Elbow plot           |
| `fig3h_mfuzz_clustering.R`   | Soft clustering            |
| `fig3i_mfuzz_go.R`           | Mfuzz GO enrichment        |
| `fig3j_mfuzz_kegg.R`         | Mfuzz KEGG enrichment      |

---

## Figure 4: UPR DEG Analysis

**Statistical Model:** Centralized Core Model (`core_dds.rds`)

| Script                         | Content                               |
| ------------------------------ | ------------------------------------- |
| `fig4_upr_shared.R`          | Shared UPR DEG functions              |
| `fig4_volcano_shared.R`      | Extracts 6h contrasts from Core Model |
| `fig4a_volcano_wt_6h.R`      | WT-Tm-6h vs WT-YPD-6h                 |
| `fig4b_volcano_ko_6h.R`      | hac1Δ-Tm-6h vs hac1Δ-YPD-6h         |
| `fig4c_venn.R`               | Venn Diagram                          |
| `fig4d_venn_direction.R`     | Venn Diagram Direction                |
| `fig4e_go_direction.R`       | GO Enrichment Direction               |
| `fig4f_kegg_upr_direction.R` | KEGG Enrichment Direction             |

---

## Output Directories

| Directory                          | Contents                       |
| ---------------------------------- | ------------------------------ |
| `results/fig1_qc/`               | PCA & QC plots                 |
| `results/fig2_deg_treatment_0h/` | Baseline 0h plots              |
| `results/fig3_deg_timecourse/`   | Temporal & Clustering plots    |
| `results/fig4_deg_upr/`          | UPR Pathway plots              |
| `results/deseq2_results/`        | DESeq2 RDS objects & CSV files |

---

## Parameters

| Parameter       | Value |
| --------------- | ----- |
| `PADJ_CUTOFF` | 0.05  |
| `LFC_CUTOFF`  | 1.0   |
| Base font size  | 12pt  |
| Figure DPI      | 300   |

---

## Excluded Samples

| Sample           | Reason  |
| ---------------- | ------- |
| `1-WT-YPD-12`  | Outlier |
| `2-WT-YPD-12`  | Outlier |
| `3-WT-YPD-12`  | Outlier |
| `1-WT-YPD-24`  | Outlier |
| `1-hac1-YPD-6` | Outlier |
| `3-WT-0`       | Outlier |
| `1-WT-Tm-12`   | Outlier |
| `1-WT-Tm-6`    | Outlier |


---

## License

This project is licensed under the **MIT License** — see [`LICENSE`](LICENSE) for details.

---

**Author:** Thanatan W.
**Date:** May 2026

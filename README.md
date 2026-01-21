# ARFD RNA-seq Analysis

This repository contains an R Markdown analysis pipeline for differential expression analysis of ARF (Auxin Response Factor) deletion mutants.

## Overview

The analysis pipeline performs:
1. **Data import** from featureCounts output
2. **DESeq2 normalization** and differential expression analysis
3. **DEG classification** (Up/Down regulated genes)
4. **Heatmaps** with hierarchical clustering
5. **Volcano plots** for each comparison
6. **UpSet plots** showing gene overlaps
7. **GO term enrichment analysis**
8. **Correlation analysis** for finding co-regulated genes

## Requirements

### R Packages

```r
# Bioconductor packages
BiocManager::install(c("DESeq2", "org.At.tair.db", "clusterProfiler", "GOSemSim"))

# CRAN packages
install.packages(c("dplyr", "ggplot2", "pheatmap", "svglite", "ggrepel", "tidyr", "UpSetR", "eulerr"))
```

## Setup Instructions

### 1. Add your featureCounts file

Copy your featureCounts output file to this directory and name it `gene_counts.txt`:

```bash
cp /path/to/your/gene_counts.txt ./gene_counts.txt
```

### 2. Verify sample metadata

The `sample_metadata.csv` file contains the mapping between sample numbers and construct names:

| Sample | Construct | Description |
|--------|-----------|-------------|
| 01, 08, 15, 22 | GUS | Control |
| 02, 09, 16, 23 | GUS_aux | GUS + auxin |
| 03, 10, 17, 24 | 5D | ARF5 delta |
| 04, 11, 18, 25 | 6D | ARF6 delta |
| 05, 12, 19, 26 | 7D | ARF7 delta |
| 06, 13, 20, 27 | 8D | ARF8 delta |
| 07, 14, 21, 28 | 19D | ARF19 delta |

Modify `sample_metadata.csv` if your sample layout is different.

### 3. Configure analysis parameters

In `ARFD_RNAseq_Analysis.Rmd`, modify the Configuration section:

```r
# DEG thresholds
padj_threshold <- 0.01      # Adjusted p-value cutoff
lfc_threshold <- 2          # Log2 fold change cutoff
min_avg_reads <- 20         # Minimum average normalized reads filter
```

### 4. Run the analysis

Open the R Markdown file in RStudio and click "Knit", or run from command line:

```r
rmarkdown::render("ARFD_RNAseq_Analysis.Rmd")
```

## Output Files

The analysis generates:

| File | Description |
|------|-------------|
| `ARFD_normalized_counts_with_averages.csv` | Normalized counts with per-construct averages |
| `ARFD_DE_results_LFC2_p01.csv` | Full DE results with Up/Down classification |
| `ARFD_DE_results_LFC2_p01_clustered.csv` | DE results with gene cluster assignments |
| `gene_clusters.csv` | Gene-to-cluster mapping |
| `All_DEGs_heatmap.svg` | Heatmap of all DEGs |
| `Volcano_*.svg` | Volcano plots for each comparison |
| `UpSet_UP_regulated.svg` | UpSet plot of up-regulated genes |
| `UpSet_DOWN_regulated.svg` | UpSet plot of down-regulated genes |
| `GO_enrichment_UP.svg` | GO term enrichment bubble plot |

## Experimental Design

- **Control**: GUS (samples 01, 08, 15, 22)
- **Treatments**: GUS_aux, 5D, 6D, 7D, 8D, 19D (ARF deletion mutants)
- **Replicates**: 4 biological replicates per condition
- **Total samples**: 28

## Analysis Details

### Differential Expression Criteria

Genes are classified as differentially expressed if:
- Adjusted p-value < 0.01
- |log2 fold change| > 2
- Average normalized reads > 20 in either treatment or control

### Clustering

Hierarchical clustering is performed on DEGs using row-scaled (z-score) average normalized reads. The optimal number of clusters is determined using the elbow method.

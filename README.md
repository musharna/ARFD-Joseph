# ARFD RNA-seq Analysis

This repository contains a comprehensive R Markdown analysis pipeline for differential expression analysis of ARF (Auxin Response Factor) deletion mutants, including downstream motif analysis, external dataset integration, and transcription factor enrichment.

## Overview

The analysis pipeline performs:
1. **Data import** from featureCounts output
2. **DESeq2 normalization** and differential expression analysis
3. **DEG classification** (Up/Down regulated genes)
4. **Heatmaps** with hierarchical clustering
5. **Volcano plots** for each comparison
6. **UpSet plots** and **Venn diagrams** showing gene overlaps
7. **GO term enrichment analysis** with semantic similarity clustering
8. **Cluster-specific GO enrichment**
9. **Correlation analysis** for finding co-regulated genes
10. **Promoter motif analysis** (AuxRE clusters, repeat orientation)
11. **Coupling element identification**
12. **SEA motif enrichment comparison**
13. **TF Deacon enrichment analysis**
14. **External dataset integration** (Okushima, cell-type specific)

## Requirements

### R Packages

```r
# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "org.At.tair.db",
  "clusterProfiler",
  "GOSemSim",
  "Biostrings"
))

# CRAN packages
install.packages(c(
  "dplyr",
  "ggplot2",
  "pheatmap",
  "svglite",
  "ggrepel",
  "tidyr",
  "UpSetR",
  "eulerr",
  "ggVennDiagram",
  "scales"
))
```

### External Tools (Optional for Advanced Analyses)

- **MEME Suite** (for MCAST and SEA motif analysis): http://meme-suite.org/
- **TF Deacon**: For transcription factor family enrichment analysis

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

### 4. Optional: Add external datasets (for advanced analyses)

Create the following directory structure for optional analyses:

```
ARFD-Joseph/
├── ARFD_RNAseq_Analysis.Rmd
├── gene_counts.txt
├── Clustering/
│   ├── Motif_Analysis/
│   │   ├── mcast_output.tsv (MCAST results)
│   │   └── Coupling_Elements/
│   │       └── *_coupling_SEA*.tsv (SEA results)
│   └── promoter_sequences.fasta (promoter sequences for motif analysis)
├── Other_Datasets/
│   ├── Okushima2005.csv (published ARF mutant data)
│   ├── Auxin_Induced_Okushima.csv (known auxin-responsive genes)
│   └── Cell_Type_Specific_Data_auxin_responsive_BB.csv
└── TF_Deacon/
    └── ARFs/
        └── *.csv (TF Deacon output files)
```

**Notes:**
- These external files are **optional** and only required for specific downstream analyses
- The main RNA-seq analysis (sections 1-8) will run without these files
- Advanced analyses (motif, external datasets, TF enrichment) will be skipped if files are not present
- Set `eval=FALSE` chunks to `eval=TRUE` when you have the required files

#### Getting External Datasets

1. **Okushima 2005 data**: Download from GEO (GSE627)
2. **Cell-type specific data**: Contact original authors or use published supplementary data
3. **Promoter sequences**: Extract using `bedtools` or BioconductorR packages
4. **MCAST output**: Run MEME Suite MCAST on promoter sequences with AuxRE motif
5. **SEA output**: Run MEME Suite SEA on coupling element sequences
6. **TF Deacon**: Run DEACoN analysis on your gene lists

### 5. Run the analysis

Open the R Markdown file in RStudio and click "Knit", or run from command line:

```r
rmarkdown::render("ARFD_RNAseq_Analysis.Rmd")
```

**Note:** The Rmd file has several code chunks marked with `eval=FALSE` for optional analyses. These will be skipped automatically if required files are not present.

## Output Files

### Core RNA-seq Outputs

| File | Description |
|------|-------------|
| `ARFD_normalized_counts_with_averages.csv` | Normalized counts with per-construct averages |
| `ARFD_DE_results_LFC2_p01.csv` | Full DE results with Up/Down classification |
| `ARFD_DE_results_LFC2_p01_clustered.csv` | DE results with gene cluster assignments |
| `gene_clusters.csv` | Gene-to-cluster mapping |
| `All_DEGs_heatmap.svg` | Heatmap of all DEGs with clustering |
| `All_DEGs_heatmap_clustered.svg` | Heatmap with cluster annotations |
| `Auxin_Responsive_heatmap.svg` | Heatmap of auxin-responsive genes |
| `K-means_all_DEGs.svg` | Elbow plot for cluster number selection |

### Visualization Outputs

| File | Description |
|------|-------------|
| `Volcano_*.svg` | Volcano plots for each construct comparison |
| `UpSet_UP_regulated.svg` | UpSet plot of up-regulated gene overlaps |
| `UpSet_DOWN_regulated.svg` | UpSet plot of down-regulated gene overlaps |
| `Venn_UP_regulated.svg` | Euler diagram of up-regulated genes |
| `Venn_DOWN_regulated.svg` | Euler diagram of down-regulated genes |
| `VennDiagram_UP_ggVenn.svg` | Alternative Venn diagram (4-way comparison) |
| `VennDiagram_DOWN_ggVenn.svg` | Alternative Venn diagram (4-way comparison) |

### GO Enrichment Outputs

| File | Description |
|------|-------------|
| `GO_enrichment_UP_semantically_clustered.svg` | GO terms clustered by semantic similarity |
| `GO_enrichment_by_construct.csv` | GO enrichment results per construct |
| `GO_enrichment_by_cluster.svg` | GO enrichment for gene clusters |
| `GO_enrichment_by_cluster.csv` | GO cluster enrichment data |

### Motif Analysis Outputs (Optional)

| File | Description |
|------|-------------|
| `Motif_Analysis_All_Hits.csv` | All AuxRE motif hits with positions |
| `Motif_Analysis_Summary.csv` | Summary of motifs per cluster |
| `Repeat_Orientation_Barplot.svg` | DR/IR/ER orientation patterns |
| `auxre_coupling_elements.fasta` | Extracted coupling element sequences |
| `SEA_Motif_Enrichment_Volcano.svg` | Differential motif enrichment |
| `SEA_Motif_Enrichment_Barplot.svg` | Top enriched motifs |
| `SEA_Motif_Comparison.csv` | Motif enrichment comparison data |

### TF Enrichment Outputs (Optional)

| File | Description |
|------|-------------|
| `TF_Deacon_Enrichment.svg` | TF family enrichment heatmap |
| `TF_Deacon_Summary.csv` | TF enrichment statistics |

### External Dataset Outputs (Optional)

| File | Description |
|------|-------------|
| `Cell_Type_Specific_Expression.svg` | DEGs in cell-type specific data |
| `Okushima_Overlap_Genes.csv` | Genes overlapping with Okushima data |
| `Dataset_Overlap_Venn.svg` | Venn diagram of dataset overlaps |

## Experimental Design

- **Control**: GUS (samples 01, 08, 15, 22)
- **Treatments**: GUS_aux, 5D, 6D, 7D, 8D, 19D (ARF deletion mutants)
- **Replicates**: 4 biological replicates per condition
- **Total samples**: 28

## Analysis Details

### Differential Expression Criteria

Genes are classified as differentially expressed if:
- Adjusted p-value < 0.01 (default, configurable)
- |log2 fold change| > 2 (default, configurable)
- Average normalized reads > 20 in either treatment or control

### Clustering

Hierarchical clustering is performed on DEGs using row-scaled (z-score) average normalized reads. The optimal number of clusters is determined using the elbow method (k-means WSS).

### GO Term Enrichment

- Performed using `clusterProfiler` with Arabidopsis annotation (`org.At.tair.db`)
- GO terms are clustered by **semantic similarity** using `GOSemSim` (Wang method)
- Enrichment is calculated both **per construct** and **per gene cluster**
- Multiple testing correction: Benjamini-Hochberg (FDR)

### Motif Analysis Workflow (Optional)

1. **Extract promoter sequences** (-1000 to +200 bp from TSS)
2. **Run MCAST** to find AuxRE motif clusters (TGTCNN)
3. **Parse MCAST output** to identify individual motifs
4. **Calculate spacing** between adjacent motifs
5. **Classify orientation**: DR (direct repeat), IR (inverted repeat), ER (everted repeat)
6. **Identify coupling elements**: clusters of ≥2 AuxREs within 30bp
7. **Run SEA** on coupling elements to find enriched TF binding sites
8. **Compare enrichment** between gene clusters

### TF Deacon Analysis (Optional)

- Identifies enriched transcription factor families in DEG sets
- Uses Fisher's combined p-value to assess family-level significance
- Normalizes by family size (fraction of significant TFs)
- Hierarchical clustering of TF families across gene lists

### External Dataset Integration (Optional)

- **Okushima 2005**: Compare with published ARF7/ARF19 mutant data
- **Cell-type specific**: Overlay DEGs with tissue/cell-type expression
- **Statistical overlap testing**: Fisher's exact test for dataset comparisons
- **Venn diagrams**: Visualize overlaps between datasets

## Troubleshooting

### Common Issues

1. **"undefined columns selected" error**
   - Ensure `gene_counts.txt` exists in the project directory
   - Check that column names match expected pattern (default: `24250R-01-XX_*`)

2. **Missing packages**
   ```r
   # Install all required packages
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("DESeq2", "org.At.tair.db", "clusterProfiler", "GOSemSim", "Biostrings"))
   install.packages(c("dplyr", "ggplot2", "pheatmap", "svglite", "ggrepel", "tidyr", "UpSetR", "eulerr", "ggVennDiagram", "scales"))
   ```

3. **Out of memory errors**
   - Reduce the number of genes analyzed (increase `min_avg_reads` threshold)
   - Process fewer constructs at once
   - Increase R memory limit: `memory.limit(size=16000)` (Windows)

4. **GO enrichment fails**
   - Ensure at least 5 genes per category
   - Check internet connection (for GO database access)
   - Verify gene IDs are in TAIR format (AT1G01010, etc.)

5. **Optional analyses not running**
   - Check that `eval=FALSE` chunks have required input files
   - Create necessary directory structure
   - Set `eval=TRUE` when files are ready

## Citation

If you use this pipeline, please cite the relevant packages:
- **DESeq2**: Love, Huber, and Anders (2014) Genome Biology
- **clusterProfiler**: Yu et al. (2012) OMICS
- **GOSemSim**: Yu (2010) Bioinformatics
- **MEME Suite**: Bailey et al. (2015) Nucleic Acids Research

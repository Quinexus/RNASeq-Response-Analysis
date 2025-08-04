# RNASeq Pipeline

This repository contains an R-based pipeline for RNA-Seq response analysis, including differential expression, pathway enrichment, deconvolution, and machine learning.

## Features

- **Differential Expression Analysis:** Uses DESeq2 for responder/nonresponder comparisons.
- **Gene Annotation:** Supports Ensembl, Entrez, and gene symbol annotation via biomaRt and AnnotationDbi.
- **Pathway Enrichment:** GSEA for Hallmark, NABA, and custom gene sets.
- **Bulk Deconvolution:** xCell, EPIC, MCPcounter, Decosus methods.
- **Matrisome Index Calculation:** For ECM-related gene signatures.
- **Automated Reporting:** Generates annotated PDF reports with plots and tables.
- **Combination Analysis:** Integrates results across multiple datasets.
- **Machine Learning:** Model training and feature selection for response prediction.

## Usage

1. **Install Required Packages:**  
   The pipeline uses Bioconductor and CRAN packages. See the top of `Pipeline.R` for the full list.

2. **Prepare Input Data:**  
   - Counts matrix (genes x samples)
   - Metadata (samples x attributes)

3. **Run Analysis:**  
   Source `Pipeline.R` and use the main functions:
   - `response_analysis()` for single dataset analysis
   - `analyse_datasets()` for batch analysis
   - `printer_analysis()` for report generation

4. **Customisation:**  
   Edit variables at the top of `Pipeline.R` to adjust gene lists, pathways, and plotting options.

## Example

```r
source("Pipeline.R")

# Example usage
results <- response_analysis(
  name = "ExampleDataset",
  counts = "path/to/counts.csv",
  metadata = "path/to/metadata.csv",
  species = "human",
  annotation = "biomart"
)

printer_analysis(results$response)
```

## Output

- Annotated tables and plots (volcano, PCA, heatmaps, pathway dotplots)
- PDF reports in `printer_analyses/`
- Combined analyses in `combined_analyses/`

## File Required

- `Pipeline.R` – Main pipeline script
- `Matrisome_Hs_MasterList_SHORT.csv` – ECM gene list
- `InnateDB_genes.csv` – Immuno gene list
- `bagaev.csv` – Bagaev gene signatures

# TCGA PAAD Analysis

The repository includes `tcga_paad.R`, an R script for downloading, processing, and analyzing pancreatic adenocarcinoma (PAAD) data from TCGA using the TCGAbiolinks package. This script supports:

- **Data Download & Filtering:** Automated download and filtering of TCGA PAAD gene expression and clinical data.
- **Survival Analysis:** Functions for Kaplan-Meier survival analysis by gene, gene set, or signature score.
- **NMF Clustering:** Non-negative matrix factorization (NMF) for sample clustering, cluster-defining gene identification, and downstream analyses.
- **Cluster Prediction:** Assigns clusters to new datasets using gene signatures or NNLS.
- **Pathway Analysis:** GSEA for cluster-specific gene lists.

See `tcga_paad.R` for usage examples and function

# Pseudobulk Preparation

For single-cell RNA-seq data, you can generate pseudobulk counts and metadata using the provided Python script in `pseudobulk_prep/pseudobulk-prep.py`. This script uses [Scanpy](https://scanpy.readthedocs.io/) to aggregate single-cell counts by sample, producing a counts matrix and metadata compatible with the R pipeline.

**Usage:**

1. Prepare your AnnData object (`adata`) with sample information in `adata.obs['sample']`.
2. Edit `required_metadata_obs` to include relevant metadata columns.
3. Run the script to generate:
   - `pseudobulk_counts.csv` (genes x samples)
   - `pseudobulk_metadata.csv` (samples x attributes)

These files can be used as input for the main RNASeq pipeline.

See `pseudobulk_prep/pseudobulk-prep.py

# GEO Database Download

To automate downloading datasets from the NCBI GEO database, use the provided Python script in `GEO database Download/GSE_download.py`. This script connects to the GEO FTP server and downloads all files for a specified GSE accession.

**Usage:**

```bash
python GEO_download/GSE_download.py <GSE_accession>
```

Replace `<GSE_accession>` with the desired GEO Series accession (e.g., `GSE123456`). The script will create a folder for the accession and download all associated files, decompressing `.gz` files automatically.

See `GEO_download/GSE_download.py` for details.
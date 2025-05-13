# RNA Seq Pipeline Useful Functions

# 1. Read in data - counts matrix + metadata
  # Create function to prepare pseudobulk data (separate in python)
# 2. Differential Expression Analysis
  # Using DESeq2 (most important for now)
  # Using edgeR, limma and others (for extra functionality)
# 3. Perform gene annotation (add gene counts for genes with same gene symbol)
# 4. Visualise PCA
# 5. Visualise volcano plot
# 6. Pathway Analysis
# 7. Deconvolution (not important for now)

# Load in required packages
check_install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (source == "CRAN") {
      install.packages(pkg)
    } else if (source == "Bioconductor") {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else if (source == "source" && !is.null(url)) {
      install.packages(url, repos = NULL, type = "source")
    } else if (source == "github" && !is.null(url)) {
      if (!require("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github(url)
    }
  }
  library(pkg, character.only = TRUE)
}

packages <- c("tidyverse", "DESeq2", "EnhancedVolcano")

for (package in packages) {
  check_install_and_load(package)
}

# 1. Read in data
load_rna_dataset <- function(counts, metadata) {
  if (is.character(counts)) read.csv(counts)
  if (is.character(metadata)) read.csv(metadata)
  return(list(counts = counts, metadata = metadata))
}

# 2. Differential Expression Analysis
# Using DESeq2
run_deseq <- function(counts, metadata, design) {
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = filtered_bcr_drug_data,
                                design = design)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  dds <- DESeq(dds)

  res <- results(dds)
  
  return(list(res=res, dds=dds))
}

# 3. Gene Annotation
# Using biomart
biomart_annotation <- function(res, species) {
  
}
# Using AnnotationDbi
annodbi_annotation <- function(res, species) {
  
}

# 4/5. Visualise Data
visualise_data <- function(dds, res, design) {
  vsd <- vst(dds, blind=FALSE)
  # plotPCA(vsd, intgroup=c(as.character(design)))
  
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj')
}


# THE PIPELINE
rna_seq_pipeline <- function(counts, metadata, design, species) {
  raw_data <- load_rna_dataset(counts, metadata)
  de_res <- run_deseq(raw_data$counts, raw_data$metadata, design)
  anno_res <- biomart_annotation(de_res$res, species)
  visualise_data(de_res$dds, de_res$res, design)
}





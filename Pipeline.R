# RNA Seq Pipeline Useful Functions

# Load in required packages
library(DESeq2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(tidyverse)

# To mouse case
to_mouse_case <- function(symbols) {
  sapply(symbols, function(g) {
    if (nchar(g) == 0) return("")
    first <- substr(g, 1, 1)
    rest <- tolower(substr(g, 2, nchar(g)))
    paste0(first, rest)
  })
}


# 1. Read in data
load_rna_dataset <- function(counts, metadata, is_txi=FALSE) {
  if (is.character(counts)) read.csv(counts, row.names = 1)
  if (is.character(metadata)) read.csv(metadata, row.names = 1)
  
  if (is_txi) counts_data <- counts[["counts"]] %>% as.data.frame()
  else counts_data <- counts
  
  metadata <- metadata[colnames(counts_data), ]
  stopifnot(all(colnames(counts_data) == rownames(metadata)))
  
  return(list(counts = counts, metadata = metadata))
}

# 2. Gene Annotation
# Using biomart
biomart_annotation <- function(ensembl_list, species) {
  # Determine dataset needed
  if (species == "human") ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  else if (species == "mouse") ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  else stop("Species must be mouse or human!")
  
  # Remove version number
  ids <- sub("\\..*", "", ensembl_list)
  
  # Check if transcript or gene
  id_type <- ifelse(all(grepl("^ENST", ids)), "ensembl_transcript_id", "ensembl_gene_id")
  
  # Attributes needed
  attributes <- c(id_type, "external_gene_name")
  
  # Run biomaRt
  gene_names <- getBM(
    attributes = attributes,
    filters = id_type,
    values = ids,
    mart = ensembl
  )
  
  # Ensure all input IDs are returned, even if not annotated
  df_input <- data.frame(input_id = ids, original = ensembl_list)
  df_merged <- merge(df_input, gene_names, by.x = "input_id", by.y = id_type, all.x = TRUE)
  df_merged$gene_name <- df_merged$external_gene_name

  return(df_merged)
}

# Using AnnotationDbi
annodbi_annotation <- function(ensembl_list, species) {
  # Remove version numbers
  ids <- sub("\\..*", "", ensembl_list)
  
  if(species == "human") orgdb <- org.Hs.eg.db
  else if (species == "mouse") orgdb <- org.Mm.eg.db
  else stop("Species must be mouse or human!")
  
  # Annotate only valid Ensembl IDs
  annotations_orgDb <- AnnotationDbi::select(
    orgdb,
    keys = ids,
    columns = c("ENSEMBL", "SYMBOL"),
    keytype = "ENSEMBL"
  )
  
  # If duplicates exist (one ENSEMBL maps to many symbols), keep first
  annotations_orgDb <- annotations_orgDb[!duplicated(annotations_orgDb$ENSEMBL), ]
  
  df_input <- data.frame(input_id = ids, original = ensembl_list, stringsAsFactors = FALSE)
  df_merged <- merge(df_input, annotations_orgDb, by.x = "input_id", by.y = "ENSEMBL", all.x = TRUE)
  df_merged$gene_name <- df_merged$SYMBOL
  
  return(df_merged)
}

# Merge res with annotation of choice
merge_annotate <- function(df, species, method = "annodbi") {
  # Get ensemble ids
  df$original <- rownames(df)
  
  # Run annotation
  if (method == "biomart") annotation <- biomart_annotation(rownames(df), species)
  else if (method == "annodbi") annotation <- annodbi_annotation(rownames(df), species)
  else if (method == "none") {
    df$gene_name <- rownames(df)
    return(df)
  }
  else stop("method not available")
  
  # merge based on original unedited ensembl ids
  merged_df <- df %>%
    dplyr::left_join(
      annotation %>% dplyr::select(original, gene_name),
      by = "original"
    )
  
  # restore rownames
  rownames(merged_df) <- merged_df$original
  merged_df$original <- NULL
  
  return(merged_df)
}

# 3. Differential Expression Analysis

# Using DESeq2

## Prepping dds
prepare_dds <- function(counts, metadata, design, is.txi=FALSE, use_custom_function=TRUE, parameter="untreated") {
  # Use given formula or create based on column
  if (is_formula(design)) formula <- design
  else if (!use_custom_function) formula <- as.formula(paste0("~ ", design))
  else {
    design_data <- design_formula(metadata, design, parameter)
    formula <- design_data$formula
    metadata <- design_data$metadata
  }
  
  # Prepare DDS object
  if (!is.txi) dds <- DESeqDataSetFromMatrix(counts, metadata, formula)
  else dds <- DESeq2::DESeqDataSetFromTximport(txi, metadata, formula)
  
  if(use_custom_function & !is.null(design_data$reference)) dds$group <- relevel(dds$group, ref=design_data$reference)
  
  # Remove low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  return(dds)
}

design_formula <- function(metadata, design, reference) {
  if (design == "response") {
    formula <- as.formula("~ response")
    reference <- NULL
  }
  else if (design == "response + treatment") {
    metadata$group <- paste0(tolower(metadata$response), tolower(metadata$treatment))
    formula <- as.formula("~ group")
    reference <- paste0("nonresponder", reference)
  }
  return(list(formula=formula, metadata=metadata, reference=reference))
}

## Running DESeq2
run_deseq <- function(dds, species, method) {
  # Run DESeq2
  dds <- DESeq(dds)
  resnames <- print(resultsNames(dds))
  
  if (length(resnames) == 2) title <- resnames[2] else title <- "group_responderuntreated_vs_nonresponderuntreated"
  
  res <- results(dds, name = title)
  res <- res[order(res$padj),] %>% as.data.frame() %>% merge_annotate(species, method)
  return(list(res=res, dds=dds, title=title))
}

## Analysis Pipeline
rna_seq_analysis <- function(counts, metadata, design, 
                             species="human",
                             is.txi=FALSE, 
                             annotation="biomart", 
                             treatment="untreated",
                             plot.volcano=FALSE,
                             plot.pca=FALSE) {
  
  raw_data <- load_rna_dataset(counts, metadata, is.txi)
  dds <- prepare_dds(raw_data$counts, raw_data$metadata, design, is.txi, parameter=treatment)
  de_res <- run_deseq(dds, species, annotation)
  
  if(plot.volcano) visualise_volcano(de_res)
  if(plot.pca) visualise_pca(de_res)
  
  return(de_res)
}

# 4. Visualise Data

## Volcano plot
visualise_volcano <- function(de_res, proteins=NULL) {
  res <- de_res$res
  
  if (!is.null(proteins)) res <- res[rownames(res) %in% proteins, ]
  
  plot(EnhancedVolcano(de_res$res, lab = de_res$res$gene_name, x = 'log2FoldChange', y = 'padj', title = de_res$title, pCutoff = 0.05))
}

## PCA plot
visualise_pca <- function(de_res) {
    vsd <- vst(de_res$dds, blind=FALSE)
    plot(plotPCA(vsd, intgroup=c("response")))
}

visualise_heatmap <- function(de_res) {
  
}

view_res_table <- function(de_res, direction=NULL) {
  res <- de_res$res %>% as.data.frame()
  res <- res[!is.na(res$padj),]

  if(!is.null(direction)) {
    if (tolower(direction) == "up") res <- res[res$log2FoldChange > 1,]
    else res <- res[res$log2FoldChange < -1,]
  }
  print(res[res$padj <= 0.05,])
}



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

# 1. Read in data
load_rna_dataset <- function(counts, metadata) {
  if (is.character(counts)) read.csv(counts, row.names = 1)
  if (is.character(metadata)) read.csv(metadata, row.names = 1)
  
  metadata <- metadata[colnames(counts), ]
  
  stopifnot(all(colnames(counts) == rownames(metadata)))
  return(list(counts = counts, metadata = metadata))
}

# 2. Gene Annotation
# Using biomart
biomart_annotation <- function(ensembl_list, species, use_species_attr=FALSE) {
  # Determine dataset needed
  if (species == "human") {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    spec_attr <- "hgnc_symbol"
  }
  else if (species == "mouse") {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    spec_attr <- "mgi_symbol"
  }
  else stop("Species must be mouse or human!")
  
  # Remove version number
  ids <- sub("\\..*", "", ensembl_list)
  
  # Check if transcript or gene
  id_type <- ifelse(all(grepl("^ENST", ids)), "ensembl_transcript_id", "ensembl_gene_id")
  
  # Attributes needed
  attributes <- ifelse(use_species_attr, c(id_type, spec_attr), c(id_type, "external_gene_name"))
  
  # Run biomaRt
  gene_names <- getBM(
    attributes = attributes,
    filters = id_type,
    values = ids,
    mart = ensembl
  )
  
  print(gene_names)
  
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
  else orgdb <- org.Mm.eg.db
  
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
  else if (method == "none") return(df)
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
prepare_dds <- function(counts, metadata, design) {
  # Use given formula or create based on column
  if (!is_formula(design)) formula <- as.formula(paste0("~ ", design))
  else formula <- design
  
  # Prepare DDS object
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = formula)
  
  # Remove low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  return(dds)
}

## Prepping tximport
prepare_txi <- function(txi, metadata, design) {
  dds <- DESeq2::DESeqDataSetFromTximport(txi, metadata, design)
  return(dds)
}

## Running DESeq2
run_deseq <- function(dds, print_res_names = FALSE, contrast=NULL) {
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Print results names to see possible contrasts
  if(print_res_names) print(resultsNames(dds))
  
  if(is.null(contrast)) res <- results(dds)
  else res <- results(dds, contrast = contrast)
  
  res <- res[order(res$padj),]

  return(list(res=res, dds=dds))
}

# 4. Visualise Data

## Volcano plot
visualise_volcano <- function(res, title="Volcano Plot") {
  if("gene_name" %in% colnames(res)) labels <- res$gene_name else labels <- rownames(res)
  plot(EnhancedVolcano(res, lab = labels, x = 'log2FoldChange', y = 'padj', title = title))
}

## PCA plot
visualise_pca <- function(dds, design) {
  if (!is_formula(design)) {
    vsd <- vst(dds, blind=FALSE)
    plot(plotPCA(vsd, intgroup=c(as.character(design))))
  }
  else print("VSD not plotted, please plot mannually")
}

## Processing to loop through resultsNames
process_visualise <- function(dds, species, annotation) {
  results_names <- resultsNames(dds)
  for (i in 2:(length(results_names))) {
    res <- results(dds, name = results_names[i]) %>% as.data.frame() %>% merge_annotate(species, annotation)
    res <- res[order(res$padj),]
    visualise_volcano(res, results_names[i])
  }
}

# THE PIPELINE
rna_seq_pipeline <- function(counts, metadata, design, species, annotation="annodbi", txi=FALSE, contrast=NULL) {
  if(!txi){
    raw_data <- load_rna_dataset(counts, metadata)
    dds <- prepare_dds(raw_data$counts, raw_data$metadata, design)
  } else {
    dds <- prepare_txi(counts, metadata, design)
  }

  de_res <- run_deseq(dds, FALSE, contrast) 
  process_visualise(de_res$dds, species, annotation)
  
  # visualise_pca(de_res$dds, design)
  return(de_res$dds)
}

# TREATMENT RESPONSE PIPELINE without interaction
treatment_response_pipeline <- function(counts, metadata, species, annotation, treatment, txi=FALSE) {
  # Prepare dds
  if (txi) dds <- prepare_txi(counts, metadata, ~response)
  else dds <- prepare_dds(counts, metadata, ~response)
  
  # Filter by treatment
  dds_treatment <- dds[, dds$treatment == treatment]
  dds_treatment[[treatment]] <- droplevels(dds_treatment$response)
  design(dds_treatment) <- as.formula(paste0("~ ", treatment))
  
  # Run DESeq2
  de_res <- run_deseq(dds_treatment, FALSE, NULL) 
  
  # Plot volcano
  process_visualise(de_res$dds, species, annotation)
}

# TREATMENT RESPONSE PIPELINE with interaction
treatment_response_interaction_pipeline <- function(counts, metadata, species, annotation_method = "none", txi=FALSE) {
  # Set factor levels explicitly to control reference levels
  metadata$treatment <- factor(metadata$treatment, levels = c("Treated", "Untreated"))
  metadata$response <- factor(metadata$response, levels = c("Nonresponder", "Responder"))
  
  # Prepare DESeqDataSet with interaction model
  if(!txi) dds <- prepare_dds(counts, metadata, ~ treatment + response + treatment:response)
  else dds <- prepare_txi(counts, metadata, ~ treatment + response + treatment:response)
  
  # Run DESeq
  dds <- DESeq(dds)
  print(resultsNames(dds))  # Show available contrast names for transparency
  
  # Extract contrasts
  res_treated <- results(dds, name = "response_Responder_vs_Nonresponder")
  res_untreated <- results(dds, contrast = list(c("response_Responder_vs_Nonresponder", "treatmentUntreated.responseResponder")))
  res_interaction <- results(dds, name = "treatmentUntreated.responseResponder")
  
  # Annotate results
  res_treated <- merge_annotate(as.data.frame(res_treated), species, method = annotation_method)
  res_untreated <- merge_annotate(as.data.frame(res_untreated), species, method = annotation_method)
  res_interaction <- merge_annotate(as.data.frame(res_interaction), species, method = annotation_method)
  
  # Order by adjusted p-value
  res_treated <- res_treated[order(res_treated$padj), ]
  res_untreated <- res_untreated[order(res_untreated$padj), ]
  res_interaction <- res_interaction[order(res_interaction$padj), ]
  
  # Visualize
  visualise_volcano(res_treated, "Responder vs Nonresponder (Treated)")
  visualise_volcano(res_untreated, "Responder vs Nonresponder (Untreated)")
  visualise_volcano(res_interaction, "Interaction Effect")
  
  # Return a named list of results
  return(list(
    dds = dds,
    treated_vs_nonresponder = res_treated,
    untreated_vs_nonresponder = res_untreated,
    interaction = res_interaction
  ))
}

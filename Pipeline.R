# RNA Seq Pipeline Useful Functions

# Load in required packages
library(DESeq2)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(tidyverse)
library(grid)
library(gridExtra)
library(png)
library(enrichplot)
library(Orthology.eg.db)

# To mouse case
# Required if using homologene mapping
to_mouse_case <- function(symbols, method = "tomouse", source = "biomart") {
  if (method == "case") {
    return(sapply(symbols, function(g) {
      if (nchar(g) == 0) return("")
      first <- substr(g, 1, 1)
      rest <- tolower(substr(g, 2, nchar(g)))
      paste0(first, rest)
    }))
  }
  
  if (!method %in% c("tomouse", "tohuman")) {
    stop("Invalid method. Use 'case', 'tomouse', or 'tohuman'.")
  }
  
  if (source=="biomart") {
    human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
    
    if(method=="tomouse") {
      genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = symbols, 
                       mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
      mousex <- unique(genesV2[, 2])
    } else {
      genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = symbols, 
                       mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      humanx <- unique(genesV2[, 2])
    }
    return(unname(mapped[symbols])) 
    
  } else if (source=="annodbi") {
    horg <- org.Hs.eg.db
    morg <- org.Mm.eg.db
    orth <- Orthology.eg.db
    if(method=="tomouse") {
      humang <- mapIds(horg, symbols, "ENTREZID", "SYMBOL") %>% na.omit()
      
      mapped <- AnnotationDbi::select(orth, humang, "Mus_musculus", "Homo_sapiens")
      names(mapped) <- c("Homo_egid", "Mus_egid")
      musyhb <- AnnotationDbi::select(morg, as.character(mapped[,2]), "SYMBOL", "ENTREZID")
      return(setNames(musyhb[,2], symbols))
    } else {
      mouseg_list <- mapIds(morg, keys = symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "list")
      
      # Flatten by keeping the first valid mapping per symbol
      mouseg <- sapply(mouseg_list, function(x) if (length(x) > 0) x[1] else NA)
      
      # Now mouseg is a named vector
      valid_idx <- which(!is.na(mouseg))
      mouseg <- mouseg[valid_idx]
      symbols <- symbols[valid_idx]
      
      mapped <- AnnotationDbi::select(orth, mouseg, "Homo_sapiens","Mus_musculus")
      names(mapped) <- c("Mus_egid","Homo_egid")
      husymb <- AnnotationDbi::select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
      return(setNames(husymb[,2], symbols))
    }
  } else {
    stop("Invalid source. Use biomart or annodbi")
  }
}


# 1. Read in data
load_rna_dataset <- function(counts, metadata, is_txi=FALSE) {
  # read in csv files
  if (is.character(counts)) counts <- read.csv(counts, row.names = 1)
  if (is.character(metadata)) metadata <- read.csv(metadata, row.names = 1)
  
  # check if txi set
  if (is_txi) counts_data <- counts[["counts"]] %>% as.data.frame()
  else counts_data <- counts
  
  # order metadata
  metadata <- metadata[colnames(counts_data), ]
  # check if columns deseq2 appropriate
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
prepare_dds <- function(counts, 
                        metadata, 
                        design, 
                        is.txi=FALSE, 
                        use_custom_function=TRUE, 
                        parameter="untreated") { 
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
  else dds <- DESeq2::DESeqDataSetFromTximport(counts, metadata, formula)
  
  # relevel as required
  if(use_custom_function) {
    if(!is.null(design_data$reference)) dds$group <- relevel(dds$group, ref=design_data$reference)
  }
  
  # Remove low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  return(dds)
}

design_formula <- function(metadata, design, reference) {
  # if response return response
  if (design == "response") {
    formula <- as.formula("~ response")
    reference <- NULL
  }
  # if multifactorial create new column and deal with it
  else if (design == "response + treatment") {
    metadata$group <- paste0(tolower(metadata$response), tolower(metadata$treatment))
    formula <- as.formula("~ group")
    reference <- paste0("nonresponder", reference)
  }
  return(list(formula=formula, metadata=metadata, reference=reference))
}

## Running DESeq2
run_deseq <- function(dds, species, method, title) {
  # Run DESeq2
  dds <- DESeq(dds)
  
  # get resultsNames, results and title
  resnames <- resultsNames(dds)
  if (length(resnames) == 2) title <- resnames[2]
  res <- results(dds, name = title)
  
  # order res and annotate
  res <- res[order(res$padj),] %>% as.data.frame() %>% merge_annotate(species, method)
  
  return(list(res=res, dds=dds, title=title))
}

## Analysis Pipeline
rna_seq_analysis <- function(counts, metadata, design, 
                             species="human",
                             is.txi=FALSE, 
                             annotation= "biomart", 
                             treatment="untreated",
                             title=NULL,
                             print=FALSE) {
  
  # 1. Load in data
  raw_data <- load_rna_dataset(counts, metadata, is.txi)
  # 2. Prepare dds object
  if(is_formula(design)) design_is_forumla <- TRUE else design_is_forumla <- FALSE
  dds <- prepare_dds(raw_data$counts, raw_data$metadata, design, is.txi, parameter=treatment, use_custom_function = !design_is_forumla)
  # 3. Run deseq2 and annotate
  if(is.null(title)) title <- paste0("group_responder", treatment, "_vs_nonresponder", treatment)
  de_res <- run_deseq(dds, species, annotation, title)
  # 4. Visualise as necessary
  doPathway <- ifelse(species == "mouse", FALSE, TRUE)
  if(print) printer_analysis(de_res, do_pathway = doPathway, species = species)
    
  de_res$species <- species
  # return list of dds, res and title
  return(de_res)
}

# 4. Visualise Data

## Volcano plot
visualise_volcano <- function(de_res, proteins=NULL, plot=TRUE, save=FALSE, save_name="volcano.png") {
  # get res
  res <- de_res$res
  # use only proteins passed in
  if (!is.null(proteins)) res <- res[res$gene_name %in% proteins, ]
  
  volcanoPlot <- EnhancedVolcano(res, 
                       lab = res$gene_name, 
                       x = 'log2FoldChange', 
                       y = 'padj', 
                       title = de_res$title, 
                       pCutoff = 0.05)
  if(plot) plot(volcanoPlot)
  if(save) ggsave(save_name, volcanoPlot)
}

## PCA plot
visualise_pca <- function(de_res, plot=TRUE, save=FALSE, save_name="pca.png") {
  # get vsd
  vsd <- vst(de_res$dds, blind=FALSE)
  
  # check if treatment column exists
  if ("treatment" %in% (colData(de_res$dds) %>% as.data.frame %>% colnames())) pcaPlot <- plotPCA(vsd, intgroup=c("response", "treatment")) 
  else if ("batch" %in% (colData(de_res$dds) %>% as.data.frame %>% colnames())) pcaPlot <- plotPCA(vsd, intgroup=c("batch", "response")) 
  else pcaPlot <- plotPCA(vsd, intgroup=c("response"))
  
  if(plot) plot(pcaPlot)
  if(save) ggsave(save_name, pcaPlot)
}

## Heatmap plot
visualise_heatmap <- function(de_res, proteins=NULL, plot=TRUE, save=FALSE, save_name="heatmap.png") {
  # get relevant data
  res <- de_res$res
  vsd <- vst(de_res$dds, blind = FALSE)
  
  # filter out needed data
  if (is.null(proteins)) res <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  else res <- res[res$gene_name %in% proteins & !is.na(res$gene_name), ]
  
  # remove na and duplicates
  res <- res[!is.na(res$gene_name), ]
  res <- res[rownames(res) %in% rownames(vsd), ]
  res <- res[!duplicated(res$gene_name), ]
  
  # get expression data
  expr_data <- assay(vsd)[rownames(res), , drop = FALSE]
  rownames(expr_data) <- res$gene_name
  
  scale.dat <- t(scale(t(expr_data)))
  meta_cols <- c("response")
  if ("treatment" %in% colnames(colData(de_res$dds))) {
    meta_cols <- c(meta_cols, "treatment")
  }
  df <- as.data.frame(colData(de_res$dds)[, meta_cols, drop = FALSE])
  
  # Reorder columns by response (e.g., responders grouped together)
  ordered_idx <- order(df$response)
  scale.dat <- scale.dat[, ordered_idx]
  df <- df[ordered_idx, , drop = FALSE]
  
  heatmap <- pheatmap(scale.dat,
           cluster_rows = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_cols = FALSE,
           annotation_col = df,
           silent=!plot)
  
  if(save) save_pheatmap(heatmap, save_name)
}

# save heatmap function taken from: https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e
save_pheatmap <- function(x, filename, width=12, height=12){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(grepl(".png",filename)){
    png(filename, width=width, height=height, units = "in", res=300)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else if(grepl(".pdf",filename)){
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else{
    print("Filename did not contain '.png' or '.pdf'")
  }
}

view_res_table <- function(de_res, direction=NULL) {
  # get data
  res <- de_res$res %>% as.data.frame()
  res <- res[!is.na(res$padj),]

  # get direction
  if(!is.null(direction)) {
    if (tolower(direction) == "up") res <- res[res$log2FoldChange > 1,]
    else res <- res[res$log2FoldChange < -1,]
  } else res <- res[abs(res$log2FoldChange) > 1,]
  return(res[res$padj <= 0.05,])
}

# save df as png
df_png <- function(df, save_name="x.png", row_width=50, col_width=200) {
  png(save_name, height = row_width*nrow(df), width = col_width*ncol(df))
  grid.table(df)
  dev.off()
}

# 5. Pathway Enrichment (GSEA with MSigDB Hallmark sets)
pathway_enrichment <- function(de_res, species = "human", category = "H") {
  # get DESeq2 results
  res <- de_res$res
  # remove NAs
  res <- res[!is.na(res$gene_name) & !is.na(res$stat), ]
  # remove duplicates (keep highest stat per gene)
  res <- res[!duplicated(res$gene_name), ]

  # Create ranked gene list
  geneList <- res$stat
  names(geneList) <- res$gene_name
  geneList <- sort(geneList, decreasing = TRUE)
  
  # For enricher
  sig_genes <- res$gene_name[res$padj < 0.05 & abs(res$log2FoldChange) > 1]
  
  # Get MSigDB hallmark gene sets
  msigdb_sets <- msigdbr(species = species, collection = category)
  msigdbr_t2g <- dplyr::distinct(msigdb_sets, gs_name, gene_symbol)
  
  # Run enricher
  gseaRes <- GSEA(geneList, TERM2GENE = msigdbr_t2g)
  enrichRes <- enricher(sig_genes, TERM2GENE = msigdbr_t2g)
  
  return(gseaRes)
}

visualise_pathways <- function(de_res, species="human", plot=TRUE, save=FALSE, save_name="pathways.png") { 
 pathways <- pathway_enrichment(de_res, species)
 pathways_dotplot <- dotplot(pathways)
 
 if(plot) plot(pathways_dotplot)
 if(save) ggsave(save_name, pathways_dotplot)
}

# PNG printer pipeline
printer_analysis <- function(de_res, 
                             title="analysis", 
                             do_pathway=TRUE, 
                             extra_info=NULL, 
                             species="human") {
  dir.create(title, showWarnings = FALSE)
  metadata <- de_res$dds@colData %>% as.data.frame()
  matrisome_genes <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")
  if(species == "mouse") matrisome_genes <- to_mouse_case(matrisome_genes, method = "tomouse", source = "annodbi")
   
  ### Slide 1: Title
  title_slide <- function(title) {
    grid.newpage()
    grid.text(title, gp = gpar(fontsize = 24, fontface = "bold"))
  }
  
  ### Slide 2: Extra Info (Adding text to slide)
  text <- function(text) {
    grid.newpage()
    grid.text(text, x = 0.05, y = 0.95, just = c("left", "top"),
              gp = gpar(fontsize = 12))
  }
  
  ### Slide 3+: Heading Slide
  add_heading_slide <- function(main_title, subtitle) {
    grid.newpage()
    grid.text(main_title, y = 0.65, gp = gpar(fontsize = 20, fontface = "bold"))
    grid.text(subtitle, y = 0.45, gp = gpar(fontsize = 14))
  }
  
  ### Slide 3: metadata
  df_png(metadata, paste0(title, "/metadata.png"), col_width = 200)
  
  ### Slide 3: PCA
  visualise_pca(de_res, plot=FALSE, save=TRUE, save_name=paste0(title, "/pca.png")) # saved as pca.png
  
  ### Slide 4: Volcano
  visualise_volcano(de_res, plot=FALSE, save=TRUE, save_name=paste0(title, "/volcano.png")) # saved as volcano.png
  
  ### Slide 5: Significant Genes
  significant_genes <- view_res_table(de_res) %>% dplyr::slice(1:15)
  df_png(significant_genes, save_name=paste0(title, "/significant_genes.png"), col_width = 110)
  
  ### Slide 6: Heatmap matrisome
  visualise_heatmap(de_res, save=TRUE, proteins=matrisome_genes, save_name=paste0(title, "/heatmap_matrisome.png"), plot=FALSE)
  
  ### Slide 7: Heatmap significant
  visualise_heatmap(de_res, save=TRUE, save_name=paste0(title, "/heatmap_significant.png"), plot=FALSE)
  
  ### Slide 8: Pathway (as required)
  if (do_pathway) {
    visualise_pathways(de_res, species=species, plot=FALSE, save=TRUE, save_name=paste0(title, "/pathways.png"))
  }
  
  # Load in data
  slide_files <- c(
    "metadata.png",
    "pca.png",
    "volcano.png",
    "significant_genes.png",
    "heatmap_matrisome.png",
    "heatmap_significant.png"
  )
  if(do_pathway) slide_files <- append(slide_files, "pathways.png")
  slide_files <- file.path(title, slide_files)
  
  
  # Create the PDF!
  pdf(file = paste0(title, ".pdf"), width = 8, height = 6)
  
  ## Slide 1: Title
  title_slide(title)
  
  ## Slide 2: Extra Info
  if (!is.null(extra_info) && nzchar(extra_info)) {
    text(extra_info)
  }
  
  
  # Slide 3+: Annotated image slides
  slide_titles <- c(
    "Metadata Table",
    "PCA Plot",
    "Volcano Plot",
    "Top 15 Significant Genes",
    "Heatmap: Matrisome Genes",
    "Heatmap: All Significant Genes"
  )
  if(do_pathway) slide_titles <- append(slide_titles, "Dotplot Pathways")
  
  for (i in seq_along(slide_files)) {
    # Add heading slide
    add_heading_slide(title, slide_titles[i])
    
    # Then plot
    img <- readPNG(slide_files[i])
    grid.newpage()
    grid.draw(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")))
  }
  dev.off()
}

# Find commonly upregulated and downregulated genes
common_genes <- function(res_vector, direction="up") {
  genes <- c("s")
  for (i in seq_along(res_vector)) {
    res_genes <- view_res_table(res_vector[[i]], direction)$gene_name
    if(res_vector[[i]]$species == "mouse") res_genes <- to_mouse_case(res_genes, "tohuman", "annodbi")
    res_genes <- res_genes[!is.na(res_genes)]
    genes <- append(genes, res_genes)
  }
  gene_counts <- table(genes)
  common <- gene_counts[gene_counts > 1]
  return(sort(common, decreasing = TRUE))
}

# datasets vector format
# c(batch_name, counts_matrix=NULL, metadata=NULL, de_res = NULL)
merge_human_datasets <- function(datasets_vector, annotation="annodbi") {
  merged_counts <- NULL
  merged_metadata <- NULL
  
  for (i in seq_along(datasets_vector)) {
    # 1. Load dataset
    title <- datasets_vector[[i]]$title
    if(!is.null(datasets_vector[[i]]$deres.dds)) {
      deres <- datasets_vector[[i]]$deres.dds
      counts_matrix <- assay(deres)  # stays as matrix
      metadata <- as.data.frame(colData(deres))
    } else {
      counts_matrix <- datasets_vector[[i]]$counts
      metadata <- datasets_vector[[i]]$metadata
    }

    # 2. Reverse annotation
    counts_matrix <- reverse_annotate(counts_matrix, annotation)
    
    # 3. Filter lowly expressed genes
    keep <- rowSums(counts_matrix >= 10) >= 3
    counts_matrix <- counts_matrix[keep, ]
    
    # 4. Merge counts (by rownames = Ensembl IDs)
    if (is.null(merged_counts)) {
      merged_counts <- counts_matrix
    } else {
      # Full join and fill NAs with 0
      merged_counts <- merge(
        merged_counts,
        counts_matrix,
        by = "row.names",
        all = TRUE
      )
      rownames(merged_counts) <- merged_counts$Row.names
      merged_counts$Row.names <- NULL
      merged_counts[is.na(merged_counts)] <- 0
    }
    
    # 5. Prepare and merge metadata
    metadata$sample <- rownames(metadata)
    metadata$response <- tolower(metadata$response)
    metadata$batch <- title  # use provided title
    metadata <- metadata[, c("sample", "response", "batch")]
    
    merged_metadata <- if (is.null(merged_metadata)) metadata else rbind(merged_metadata, metadata)
  }

  # 6. Run DESeq2 on merged datasets
  deres_merged <- rna_seq_analysis(merged_counts, merged_metadata, ~batch+response, annotation="annodbi", title="response_responder_vs_nonresponder")
  
  return(deres_merged)
}

reverse_annotate <- function(counts_matrix, method) {
  # Ensure input is a matrix
  if (!is.matrix(counts_matrix)) {
    counts_matrix <- as.matrix(counts_matrix)
  }
  
  # If Ensembl IDs with version are present
  if (grepl("^ENSG", rownames(counts_matrix)[1])) {
    ids <- sub("\\..*", "", rownames(counts_matrix))  # strip versions
    counts_matrix <- rowsum(counts_matrix, group = ids)  # collapse duplicates
  } else {
    ids <- rownames(counts_matrix)
    
    if (method == "biomart") {
      ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      attributes <- c("external_gene_name", "ensembl_gene_id")
      
      ensembl_names <- biomaRt::getBM(
        attributes = attributes,
        filters = "external_gene_name",
        values = ids,
        mart = ensembl
      )
      
      ensembl_names <- ensembl_names[!duplicated(ensembl_names$external_gene_name), ]
      matched_ids <- ensembl_names$ensembl_gene_id[match(ids, ensembl_names$external_gene_name)]
    } else {
      orgdb <- org.Hs.eg.db::org.Hs.eg.db
      
      annotations_orgDb <- AnnotationDbi::select(
        orgdb,
        keys = ids,
        columns = c("ENSEMBL", "SYMBOL"),
        keytype = "SYMBOL"
      )
      
      annotations_orgDb <- annotations_orgDb[!duplicated(annotations_orgDb$SYMBOL), ]
      matched_ids <- annotations_orgDb$ENSEMBL[match(ids, annotations_orgDb$SYMBOL)]
    }
    
    # Safe annotation & duplicate handling
    counts_df <- as.data.frame(counts_matrix)
    counts_df$mapped_id <- matched_ids
    
    counts_df <- counts_df[!is.na(counts_df$mapped_id), ]
    counts_df <- aggregate(. ~ mapped_id, data = counts_df, FUN = sum)
    
    rownames(counts_df) <- counts_df$mapped_id
    counts_df$mapped_id <- NULL
    
    counts_matrix <- as.matrix(counts_df)
  }
  
  return(counts_matrix)
}

  
  
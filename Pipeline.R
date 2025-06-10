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
library(sva)
library(qpdf)

# ---- 1. Load all packages  ----
load_rna_dataset <- function(counts, metadata) {
  # check if csv then import
  if (is.character(counts)) counts <- read.csv(counts, row.names = 1)
  if (is.character(metadata)) metadata <- read.csv(metadata, row.names = 1)
  
  # check if txi set
  is_txi <- all(c("abundance", "counts", "length") %in% names(counts))
  if (is_txi) counts_data <- counts[["counts"]] %>% as.data.frame()
  else counts_data <- counts
  
  # order metadata
  metadata <- metadata[colnames(counts_data), ]
  # check if columns deseq2 appropriate
  stopifnot(all(colnames(counts_data) == rownames(metadata)))
  message("Dataset is valid!")
  
  return(list(counts = counts, metadata = metadata))
}

# ---- 2. Prepare DDS  ----

# prepare dds object
prepare_dds <- function(counts, metadata, formula, relevel_ref=NULL) {
  # internal txi check
  is_txi <- all(c("abundance", "counts", "length") %in% names(counts))
  
  # dds prep
  if (!is_txi) dds <- DESeqDataSetFromMatrix(counts, metadata, formula)
  else dds <- DESeq2::DESeqDataSetFromTximport(counts, metadata, formula)
  
  # check if needing to relevel
  if(!is.null(relevel_ref)) { dds$group <- relevel(dds$group, ref=relevel_ref) }
  
  # remove low counts
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  
  
  return(dds)
}

# prepare dds list (including handling treatment responses)
response_treatment_dds <- function(counts, metadata) {
  dds_list <- list()
  
  # response analysis
  response_dds <- prepare_dds(counts, metadata, ~response)
  dds_list$response <- response_dds
  
  # check treatment
  if (any(grepl("^treatment$", colnames(metadata) %>% tolower()))) {
    # created merged column
    metadata$group <- paste0(tolower(metadata$response), tolower(metadata$treatment))
    
    # dds for UNTREATED response vs nonresponse
    untreated_dds <- prepare_dds(counts, metadata, ~group, "nonresponderuntreated")
    dds_list$untreated <- untreated_dds
    
    # dds for TREATED response vs nonresponse
    treated_dds <- prepare_dds(counts, metadata, ~group, "nonrespondertreated")
    dds_list$treated <- treated_dds
  }
  message("dds objects prepared!")
  return(dds_list)
}

# ---- 3. Run DESeq2  ----
run_deseq <- function(dds, title=NULL) {
  # run DESeq2
  dds <- DESeq(dds)
  
  # assign title for pure response analysis
  if(is.null(title)) title <- resultsNames(dds)[2] 
  
  # get results according to requirements
  res <- results(dds, name=title)
  
  # order res
  res <- res[order(res$padj),]
  
  return(list(res=res, dds=dds, title=title))
}

# ---- 4. Annotate Datasets  ----
# check annotation status
annotation_status <- function(gene_list) {
  percent_ensembl <- sum(grepl("^ENS|ERCC", gene_list))/length(gene_list)
  if (percent_ensembl > 0.5) return("ensembl")
  else if(all(grepl("^\\d+$", gene_list))) return("entrez")
  else return("none")
}

# annotate genes list and produce list
annotate_genes <- function(gene_list, method, species) {
  annotation_status <- annotation_status(gene_list)
  
  if (annotation_status == "ensembl") {
    # remove ensembl version number
    ids <- sub("\\..*", "", gene_list)
    # check if transcript or gene
    biomart_id_type <- ifelse(all(grepl("^ENST", ids)), "ensembl_transcript_id", "ensembl_gene_id")
    annodbi_keytype <- ifelse(all(grepl("^ENST", ids)), "ENSEMBLTRANS", "ENSEMBL")
  } else if (annotation_status == "entrez") {
    ids <- gene_list
    biomart_id_type <- "entrezgene_id"
    annodbi_keytype <- "ENTREZID"
  } else return(data.frame(original = gene_list, gene_name = gene_list))
  
  
  if(method == "biomart") {
    # get appropriate mart for species
    if (species == "human") ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    else ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    # attributes needed
    attributes <- c(biomart_id_type, "external_gene_name")
    
    # run biomart
    gene_names <- getBM(
      attributes = attributes,
      filters = biomart_id_type,
      values = ids,
      mart = ensembl
    )
    
    # ensure all input IDs are returned, even if not annotated
    df_input <- data.frame(input_id = ids, original = gene_list)
    df_merged <- merge(df_input, gene_names, by.x = "input_id", by.y = biomart_id_type, all.x = TRUE)
    df_merged$gene_name <- df_merged$external_gene_name
    
    # if duplicate mapping exists, keep first
    df_merged <- df_merged[!duplicated(df_merged$input_id),]
  } else if (method == "annodbi") {
    # get org.db for species
    if(species == "human") orgdb <- org.Hs.eg.db
    else orgdb <- org.Mm.eg.db
    
    # run annodbi
    annotations_orgDb <- AnnotationDbi::select(
      orgdb,
      keys = ids,
      columns = c(annodbi_keytype, "SYMBOL"),
      keytype = annodbi_keytype
    )
    
    df_input <- data.frame(input_id = ids, original = gene_list, stringsAsFactors = FALSE)
    df_merged <- merge(df_input, annotations_orgDb, by.x = "input_id", by.y = annodbi_keytype, all.x = TRUE)
    df_merged$gene_name <- df_merged$SYMBOL
  }
  
  # if duplicate mapping exists, keep first
  df_merged <- df_merged[!duplicated(df_merged$input_id),]
  
  return(df_merged)
}

# annotate df and add annotation column (gene_name) to df
annotate_df <- function(df, method, species) {
  ids <- rownames(df)
  df$original <- ids
  
  annotation <- annotate_genes(ids, method, species)
  
  # merge based on original unedited ids
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

# ---- The Analysis Pipeline  ----
response_analysis <- function(name, counts, metadata, species="human", annotation="biomart") {
  # 1. load in data
  raw_data <- load_rna_dataset(counts, metadata)
  
  # 2. prepare dds object
  dds_list <- response_treatment_dds(counts, metadata)
  
  # 3. run deseq2
  deres_list <- list()
  deres_list$response <- run_deseq(dds_list$response)
  if (any(grepl("^treatment$", colnames(metadata) %>% tolower()))) {
    deres_list$untreated <- run_deseq(dds_list$untreated, "group_responderuntreated_vs_nonresponderuntreated")
    deres_list$treated <- run_deseq(dds_list$treated, "group_respondertreated_vs_nonrespondertreated")
  }
  message("Finished running DESeq2")
  
  # 4. annotate res and add other metadata
  deres_list <- lapply(names(deres_list), function(deres_name) {
    deres <- deres_list[[deres_name]]
    deres$species <- species
    deres$name <- paste0(name, "_", deres_name)
    
    # annotation
    deres$res <- deres$res %>% as.data.frame() %>% annotate_df(annotation, species)
    # for combination analysis
    deres$res$dataset <- name
    return(deres)
  })
  message("Annotated results")
  
  message("Completed analyses returning deres_list")
  return(deres_list)
}

# ---- 5. Visualise Data  ---- 

plot_or_save <- function(plot, deres, plot_bool, save_bool, save_name, plot_type) {
  if(plot_bool) plot(plot)
  if(is.null(save_name)) save_name <- paste0(deres$name, "_", plot_type, ".png")
  if(save_bool)  ggsave(save_name, plot)
}

# volcano plot
visualise_volcano <- function(deres, proteins=NULL, plot=TRUE, save=FALSE, save_name=NULL) {
  # get res
  res <- deres$res
  
  # use only proteins passed in
  if (!is.null(proteins)) res <- res[res$gene_name %in% proteins, ]
  
  volcano_plot <- EnhancedVolcano(res, lab=res$gene_name,
                                  x='log2FoldChange', y='padj',
                                  title=deres$title, pCutoff=0.05)
  
  plot_or_save(volcano_plot, deres, plot, save, save_name, "volcano")
}

# PCA plot
visualise_pca <- function(deres, plot=TRUE, save=FALSE, save_name=NULL) {
  # get vsd
  vsd <- vst(deres$dds, blind=FALSE)
  
  # check if treatment column exists
  if (
    all(grepl("^treatment$", 
          colData(deres$dds) %>% as.data.frame %>% colnames() %>% tolower()))
    ) pcaPlot <- plotPCA(vsd, intgroup=c("response", "treatment")) 
  else pcaPlot <- plotPCA(vsd, intgroup=c("response"))
  
  plot_or_save(pcaPlot, deres, plot, save, save_name, "pca")
}

# heatmap plot
visualise_heatmap <- function(deres, proteins=NULL, plot=TRUE, save=FALSE, save_name=NULL) {
  # get relevant data
  res <- deres$res
  vsd <- vst(deres$dds, blind=FALSE)
  
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
  if ("treatment" %in% colnames(colData(deres$dds))) {
    meta_cols <- c(meta_cols, "treatment")
  }
  
  df <- as.data.frame(colData(deres$dds)[, meta_cols, drop = FALSE])
  
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
  
  if(is.null(save_name)) save_name <- paste0(deres$name, "_heatmap", ".png")
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

# print res table
view_res_table <- function(deres, direction=NULL) {
  # get data
  res <- deres$res %>% as.data.frame()
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


# ---- 6. Pathway Enrichment ----
pathway_enrichment <- function(deres, species = "human", category = "H") {
  # get DESeq2 results
  res <- deres$res
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
  gseaRes <- GSEA(geneList, TERM2GENE = msigdbr_t2g, pvalueCutoff = 1)

  return(gseaRes)
}

dotplot_pathways <- function(deres, plot=TRUE, save=FALSE, save_name=NULL) { 
  pathways <- pathway_enrichment(deres, deres$species)
  pathways_dotplot <- dotplot(pathways, x = "GeneRatio",
                              color = "NES", size = "p.adjust")
  
  plot_or_save(pathways_dotplot, deres, plot, save, save_name, "pathways")
}

# ---- 7. Bulk Deconvolution  ----

# ---- 8. Printer Analysis  ----
printer_analysis <- function(deres, extra_info=NULL) {
  # get title
  title <- deres$name
  
  # create dir to store all printer_analyses files
  base_dir <- "printer_analyses/"
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  file_path <- paste0(base_dir, title)
  if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
  
  # get metadata from deres object
  metadata <- deres$dds@colData %>% as.data.frame()
  
  # matrisome_genes
  matrisome_file <- read_csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
  matrisome_index <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")
  if(deres$species == "mouse") {
    matrisome_index <- to_mouse_case(matrisome_index, method = "tomouse", source = "annodbi")
    matrisome_genes <- matrisome_file$mouse
  } else {
    matrisome_genes <- matrisome_file$GeneSymbol
  }
  
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
  df_png(metadata, paste0(file_path, "/metadata.png"), col_width = 200)
  
  ### Slide 4: PCA
  visualise_pca(deres, plot=FALSE, save=TRUE, save_name=paste0(file_path, "/pca.png")) # saved as pca.png
  
  ### Slide 5: Volcano
  visualise_volcano(deres, plot=FALSE, save=TRUE, save_name=paste0(file_path, "/volcano.png")) # saved as volcano.png
  
  ### Slide 6: Heatmap matrisome index
  visualise_heatmap(deres, save=TRUE, proteins=matrisome_index, save_name=paste0(file_path, "/heatmap_matrisome_index.png"), plot=FALSE)
  
  ### Slide 7: Heatmap matrisome genes
  visualise_heatmap(deres, save=TRUE, proteins=matrisome_genes, save_name=paste0(file_path, "/heatmap_matrisome.png"), plot=FALSE)
  
  
  significant_genes <- view_res_table(deres) %>% dplyr::slice(1:15)
  if(nrow(significant_genes) > 1) {
    ### Slide 8: Significant Genes
    df_png(significant_genes, save_name=paste0(file_path, "/significant_genes.png"), col_width = 110)
    
    ### Slide 9: Heatmap significant
    visualise_heatmap(deres, save=TRUE, save_name=paste0(file_path, "/heatmap_significant.png"), plot=FALSE)
    
    plot_significant <- TRUE
  } else plot_significant <- FALSE
  
  ### Slide 9: Pathways
  #visualise_pathways(deres, species=species, plot=FALSE, save=TRUE, save_name=paste0(file_path, "/pathways.png"))
  
  # Load in data
  slide_files <- c(
    "metadata.png",
    "pca.png",
    "volcano.png",
    "heatmap_matrisome_index.png",
    "heatmap_matrisome.png"
  )
  if (plot_significant) slide_files <- c(slide_files, "significant_genes.png",
                                         "heatmap_significant.png")
  slide_files <- file.path(file_path, slide_files)
  
  
  # Create the PDF!
  pdf(file = paste0(file_path, ".pdf"), width = 8, height = 6)
  
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
    "Heatmap: Matrisome Index Genes",
    "Heatmap: Matrisome Genes"
  )
  if (plot_significant) slide_titles <- c(slide_titles, "Top 15 Significant Genes",
                                          "Heatmap: All Significant Genes")
  
  for (i in seq_along(slide_files)) {
    # Add heading slide
    add_heading_slide(title, slide_titles[i])
    
    # Then plot
    img <- readPNG(slide_files[i])
    grid.newpage()
    grid.draw(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")))
  }
  dev.off()
  
  message(paste0("Completed ", deres$name))
}

recursive_printer <- function(deres_list) {
  for (i in seq_along(deres_list)) printer_analysis(deres_list[[i]])
}

# ---- 9. Combination Analysis  ----
list_analyse_print <- function(datasets_list) {
  deres_list <- lapply(datasets_list, function(dataset) {
    deres <- response_analysis(dataset$name, dataset$counts, dataset$metadata, dataset$species, "annodbi")
    recursive_printer(deres)
    return(deres)
  })
  combine_pdfs(output_file = "publicdata.pdf")
  return(deres_list)
}

matrisome_combined <- function(deres_list) {
  merged_res_list <- list()
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      if(grepl("_response$", entry$name)) {
        res <- entry$res
        if (entry$species == "human") res$gene_name <- to_mouse_case(res$gene_name, method = "case")
        merged_res_list[[entry$name]] <- res
      }
    }
  }
  
  matrisome_index <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB") %>% to_mouse_case(., source = "annodbi")
  
  combined_res_df <- do.call(rbind, merged_res_list) %>% filter(!is.na(gene_name)) %>%
    filter(gene_name %in% matrisome_index) %>%
    mutate(Significance = ifelse(is.na(padj), "NA", ifelse(padj < 0.05, "< 0.05", "≥ 0.05")))
  
  # Plot
  ggplot(combined_res_df, aes(x = gene_name, y = dataset)) +
    geom_point(aes(size = Significance, color = log2FoldChange)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5, "NA" = 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
}

gsea_combined <- function(deres_list) {
  merged_gsea_list <- list()
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      if(grepl("_response$", entry$name)) {
        gsea <- pathway_enrichment(entry, entry$species)@result
        gsea$dataset <- gsub("_response$", "", entry$name)
        merged_gsea_list[[entry$name]] <- gsea
      }
    }
  }
  
  combined_gsea_df <- do.call(rbind, merged_gsea_list) %>% 
    mutate(Significance = ifelse(p.adjust < 0.05, "< 0.05", "≥ 0.05"))
  
  ggplot(combined_gsea_df, aes(x = ID, y = dataset)) +
    geom_point(aes(size = Significance, color = NES)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- 10. Other Useful Functions  ----
# To mouse case
# Required if using homologene mapping
to_mouse_case <- function(symbols, method = "tomouse", source = "biomart") {
  if (method == "case") {
    return(sapply(symbols, function(g) {
      if (is.na(g) || nchar(g) == 0) return("")
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

# Combine pdfs
combine_pdfs <- function(pdfs_dir="printer_analyses", output_file="combined.pdf") {
  pdf_files <- list.files(pdfs_dir, pattern = "\\.pdf$", full.names = TRUE)
  if (length(pdf_files) > 0) {
    # Combine PDFs using qpdf
    pdf_combine(input = pdf_files, output = output_file)
    
    cat("PDFs combined successfully into:", output_file, "\n")
  } else {
    cat("No PDF files found in the directory.\n")
  }
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
merge_datasets <- function(datasets_vector, species="human", annotation="annodbi") {
  merged_counts <- NULL
  merged_metadata <- NULL
  
  for (i in seq_along(datasets_vector)) {
    # 1. Load dataset
    title <- datasets_vector[[i]]$title
    if(!is.null(datasets_vector[[i]]$deres)) {
      deres <- datasets_vector[[i]]$deres$dds
      counts_matrix <- assay(deres)  # stays as matrix
      metadata <- as.data.frame(colData(deres))
    } else {
      counts_matrix <- datasets_vector[[i]]$counts
      metadata <- datasets_vector[[i]]$metadata
    }

    # 2. Reverse annotation
    counts_matrix <- reverse_annotate(counts_matrix, species, annotation)
    
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
  deres_merged <- rna_seq_analysis(merged_counts, merged_metadata, ~batch+response, annotation="annodbi", title="response_responder_vs_nonresponder", species = species)
  
  return(deres_merged)
}

reverse_annotate <- function(counts_matrix, species, method) {
  # Ensure input is a matrix
  if (!is.matrix(counts_matrix)) {
    counts_matrix <- as.matrix(counts_matrix)
  }
  
  # If Ensembl IDs with version are present
  if (grepl("^ENS", rownames(counts_matrix)[1])) {
    ids <- sub("\\..*", "", rownames(counts_matrix))  # strip versions
    counts_matrix <- rowsum(counts_matrix, group = ids)  # collapse duplicates
  } else {
    ids <- rownames(counts_matrix)
    
    if (method == "biomart") {
      if(species == "human") ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      else ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      
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
      if(species == "human") orgdb <- org.Hs.eg.db
      else orgdb <- org.Mm.eg.db
      
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

read_GEO_matrix_file <- function(location) {
  lines <- readLines(location)
  empty_line_index <- which(lines=="")[1]
  df <- read.delim(location, skip=empty_line_index, row.names = NULL) %>% t() %>% as.data.frame()
  colnames(df) <- df[1, ] %>% as.character() %>% make.unique()
  colnames(df) <- sub("^!", "", colnames(df))
  required_columns <- colnames(df) %>% grepl("Sample_geo_accession|Sample_characteristics", .)
  
  df <- df[-1,] %>% select(colnames(df)[required_columns])

  return(df)
}



  
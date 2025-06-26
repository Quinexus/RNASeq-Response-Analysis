# RNA Seq Pipeline Useful Functions

# ----- Load all packages  ----

# annotation
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(msigdbr)

# Differential Expression
library(DESeq2)

# pathway analysis
library(clusterProfiler)

# deconvolution Packages
library(xCell)
library(EPIC)
#library(MCPcounter)
library(Decosus)

# batch correction
library(sva)

# data handling
library(tidyverse)
library(SummarizedExperiment)

# ggplot
library(ggh4x)
library(ggtext)
library(ggpubr)

# visualisation
library(pheatmap)
library(EnhancedVolcano)
library(enrichplot)

# printer analysis 
library(grid)
library(gridExtra)
library(png)
library(qpdf)

# ML packages
library(caret)
library(randomForest)
library(pROC)
library(FactoMineR)
library(factoextra)
library(glmnet)

# ---- 1. Load all datasets  ----
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

reverse_human_annotate <- function(counts, method="annodbi") {
  counts_df <- as.data.frame(counts)
  ids <- rownames(counts_df)
  annotation_type <- annotation_status(ids)
  
  # handle Ensembl case — already good
  if (annotation_type == "ensembl") {
    counts_df$ensembl_id <- sub("\\..*", "", ids)
    rownames(counts_df) <- counts_df$ensembl_id
    counts_df$ensembl_id <- NULL
    return(counts_df)
  }
  
  # prepare variables
  if (annotation_type == "none") {
    biomart_id_type <- "external_gene_name"
    annodbi_keytype <- "SYMBOL"
  } else if (annotation_type == "entrez") {
    biomart_id_type <- "entrezgene_id"
    annodbi_keytype <- "ENTREZID"
  }
  
  # annotation
  if (method == "biomart") {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attributes <- c(biomart_id_type, "ensembl_gene_id")
    
    gene_ids <- getBM(
      attributes = attributes,
      filters = biomart_id_type,
      values = ids,
      mart = ensembl
    )
    
    df_input <- data.frame(original = ids)
    df_merged <- merge(df_input, gene_ids, by.x = "original", by.y = biomart_id_type, all.x = TRUE)
    df_merged <- df_merged[!duplicated(df_merged$original), ]
    
  } else if (method == "annodbi") {
    annotations <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = ids,
      columns = c("ENSEMBL"),
      keytype = annodbi_keytype
    )
    
    df_input <- data.frame(original = ids)
    df_merged <- merge(df_input, annotations, by.x = "original", by.y = annodbi_keytype, all.x = TRUE)
    df_merged <- df_merged[!duplicated(df_merged$original), ]
  }
  
  # assign Ensembl as rownames
  counts_df$original <- ids
  final_df <- dplyr::left_join(counts_df, df_merged, by = "original")
  
  ensembl_col <- ifelse("ENSEMBL" %in% colnames(final_df), "ENSEMBL", "ensembl_gene_id")
  
  
  # Step 1: filter out missing Ensembl ID
  final_df <- final_df[!is.na(final_df[[ensembl_col]]), ]
  
  # Step 2: keep only the first occurrence of each Ensembl ID
  final_df <- final_df[!duplicated(final_df[[ensembl_col]]), ]
  
  # Step 3: assign Ensembl ID as rownames
  rownames(final_df) <- final_df[[ensembl_col]]
  
  # Step 4: drop helper columns
  final_df[[ensembl_col]] <- NULL
  final_df$original <- NULL
  return(final_df)
}

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

# ---- The Analysis Pipeline ----
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
    
    # get tissue for gsea combined analysis
    if (!"tissue" %in% colnames(colData(deres$dds))) tissue <- "unknown"
    else {
      tissues <- colData(deres$dds)$tissue
      if("pancreas" %in% tissues) tissue <- "pancreas"
      else tissue <- names(sort(table(tissues), decreasing = TRUE))[1]
    }
    deres$tissue <- tissue
    
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
  
  # sort res and take top 50 only
  res <- res[order(res$log2FoldChange),] %>% slice(1:50)
  
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
pathway_enrichment <- function(deres, view_all = TRUE, category = "H", subset = NULL) {
  species <- deres$species
  
  # get DESeq2 results
  res <- deres$res
  # remove NAs
  res <- res[!is.na(res$gene_name) & !is.na(res$stat), ]
  # remove duplicates (keep highest stat per gene)
  res <- res[!duplicated(res$gene_name), ]
  
  # create ranked gene list
  geneList <- res$stat
  names(geneList) <- res$gene_name
  geneList <- sort(geneList, decreasing = TRUE)
  
  # get MSigDB hallmark gene sets
  msigdb_sets <- msigdbr(species = species, collection = category)
  if (!is.null(subset)) msigdb_sets <- msigdb_sets[msigdb_sets$gs_name %in% subset, ]
  
  msigdbr_t2g <- dplyr::distinct(msigdb_sets, gs_name, gene_symbol)
  
  # set p-value cutoff depending on view_all
  pcut <- ifelse(view_all, 1, 0.05)
  
  # run GSEA
  gseaRes <- GSEA(geneList, TERM2GENE = msigdbr_t2g, pvalueCutoff = pcut)

  return(gseaRes)
}

custom_pathway_enrichment <- function(deres, view_all = TRUE) {
  species <- deres$species
  
  # get DESeq2 results
  res <- deres$res
  # remove NAs
  res <- res[!is.na(res$gene_name) & !is.na(res$stat), ]
  # remove duplicates (keep highest stat per gene)
  res <- res[!duplicated(res$gene_name), ]
  
  # create ranked gene list
  geneList <- res$stat
  names(geneList) <- res$gene_name
  geneList <- sort(geneList, decreasing = TRUE)
  
  # custom pathway
  matrisome_proteins <- read.csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
  if(species == "mouse") {
    custom_gene_set <- data.frame(
      term = matrisome_proteins$Category,
      gene = matrisome_proteins$mouse
    ) 
  } else {
    custom_gene_set <- data.frame(
      term = matrisome_proteins$Category,
      gene = matrisome_proteins$GeneSymbol
    ) 
  }
  
  # set p-value cutoff depending on view_all
  pcut <- ifelse(view_all, 1, 0.05)
  
  # run GSEA
  gseaRes <- GSEA(geneList, TERM2GENE = custom_gene_set, pvalueCutoff = pcut)
  
  return(gseaRes)
}

dotplot_pathways <- function(gseaRes, plot=TRUE, save=FALSE, save_name=NULL) { 
  pathways_dotplot <- clusterProfiler::dotplot(gseaRes, x = "GeneRatio",
                              color = "NES", size = "p.adjust")
  
  plot_or_save(pathways_dotplot, deres, plot, save, save_name, "pathways")
}

# ---- 7. Bulk Deconvolution  ----
prepare_deconvolution_matrix <- function(deres, annotation_method = "annodbi") {
  # VST-normalized expression matrix
  vsd <- vst(deres$dds, blind = TRUE)
  expr_matrix <- assay(vsd)
  
  # Annotate gene IDs to symbols
  annotated_expr <- annotate_df(as.data.frame(expr_matrix), method = annotation_method, species = "human")
  
  # Drop rows with missing symbols
  annotated_expr <- annotated_expr[!is.na(annotated_expr$gene_name), ]
  
  # Collapse to one row per gene symbol (xCell needs unique symbols)
  expr_by_symbol <- annotated_expr %>%
    group_by(gene_name) %>%
    summarise(across(where(is.numeric), max)) %>%
    column_to_rownames("gene_name") %>%
    as.matrix()
  
  if (nrow(expr_by_symbol) < 5000) warning("Low number of annotated genes; check annotation or input data.")
  
  return(expr_by_symbol)
}

xCell_deconvolution <- function(expr) {
  xcell_scores <- xCell::xCellAnalysis(expr)
  return(xcell_scores)
}

epic_deconvolution <- function(expr) {
  epic_result <- EPIC::EPIC(bulk = expr)
  return(epic_result$cellFractions)
}

mcp_deconvolution <- function(expr) {
  mcp_scores <- MCPcounter::MCPcounter.estimate(
    expr,
    featuresType = "HUGO_symbols"
  )
  
  return(mcp_scores)
}

decosus_deconvolution <- function(expr) {
  # Ensure no duplicated gene symbols
  expr <- expr[!duplicated(rownames(expr)), ]
  
  # Convert to data.frame if it's a matrix
  expr_df <- as.data.frame(expr)
  
  # Add required hgnc_symbols column
  expr_df$hgnc_symbols <- rownames(expr_df)
  
  decosus_result <- Decosus::cosDeco(x=expr_df, rnaseq = T, plot = TRUE, ext = FALSE)
  
  est_Decosus_cells <- decosus_result$main_cells
  
  return(est_Decosus_cells %>% as.matrix())
}




deconvolute_significant <- function(deres, method="xcell") {
  expr <- prepare_deconvolution_matrix(deres)
  
  # get scores and metadata
  if(method == "xcell") {
    scores <- xCell_deconvolution(expr) %>% as.data.frame()
  } else if(method == "epic") {
    scores <- epic_deconvolution(expr) %>% as.data.frame()
  } else if(method == "mcp") {
    scores <- mcp_deconvolution(expr) %>% as.data.frame()
  } else if(method == "decosus") {
    scores <- decosus_deconvolution(expr) %>% as.data.frame()
  }
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  
  # get responder status and filter
  responders <- metadata %>% filter(response == "responder") %>% rownames()
  scores_responders <- scores %>% select(responders)
  
  # get nonresponders
  non_responders <- metadata %>% filter(!response == "responder") %>% rownames()
  scores_non_responders <- scores %>% select(non_responders)
  
  # loop through cells
  cells <- c()
  pvalues <- c()
  for(i in seq_along(rownames(scores))) {
    test_stat <- t.test(as.numeric(scores_responders[i, ]), as.numeric(scores_non_responders[i, ]))
    p_value <- test_stat$p.value
    cell <- rownames(scores)[i]
    
    cells <- c(cells, cell)
    pvalues <- c(pvalues, p_value)
  }
  test_statistics <- data.frame(cell=cells, 
                                p_value=pvalues,
                                p_adj = p.adjust(pvalues, method = "BH"))
  
  mean_diff <- rowMeans(scores_responders) - rowMeans(scores_non_responders)
  test_statistics$mean_diff <- mean_diff
  
  log2_fc <- log2(rowMeans(scores_responders) + 1) - log2(rowMeans(scores_non_responders) + 1)
  test_statistics$log2_fc <- log2_fc
  
  return(test_statistics)
}



# ---- 8. Matrisome Index  ----
matrisome_index <- function(deres) {
  # get normalised counts
  norm_counts <- counts(deres$dds, normalized=TRUE)
  # sum of counts per sample
  lib_sizes <- colSums(norm_counts)
  # calculate counts per million
  cpm <- sweep(norm_counts, 2, lib_sizes, FUN = "/") * 1e6
  # log2 transform
  log2_cpm <- log2(cpm + 1)
  # annotated df
  annotated_log2_cpm_df <- annotate_df(as.data.frame(log2_cpm), "annodbi", deres$species)
  
  # genes
  positive_genes <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1")
  negative_genes <- c("AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")
  if(deres$species == "mouse") {
    positive_genes <- positive_genes %>% to_mouse_case(., "tomouse", "annodbi")
    negative_genes <- negative_genes %>% to_mouse_case(., "tomouse", "annodbi")
  }
  
  pos_ids <- rownames(annotated_log2_cpm_df)[annotated_log2_cpm_df$gene_name %in% positive_genes]
  neg_ids <- rownames(annotated_log2_cpm_df)[annotated_log2_cpm_df$gene_name %in% negative_genes]
  
  avg_pos <- colMeans(log2_cpm[pos_ids, , drop = FALSE])
  avg_neg <- colMeans(log2_cpm[neg_ids, , drop = FALSE])
  
  matrix_index <- avg_pos / avg_neg
  
  return(matrix_index)
}

sample_index <- function(deres, plot=TRUE, save=FALSE, save_name=NULL) {
  matrix_index <- matrisome_index(deres)
  
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  metadata$index <- matrix_index[rownames(metadata)]
  
  if (any(grepl("^group$", colnames(metadata)))) {
    metadata <- metadata[order(metadata$group, metadata$index), ]
    metadata$sample <- factor(rownames(metadata), levels = rownames(metadata))
    fill <- metadata$group
    
  }
  else {
    metadata <- metadata[order(metadata$response, metadata$index), ]
    metadata$sample <- factor(rownames(metadata), levels = rownames(metadata))
    fill <- metadata$response
  }
    
  sample_index_plot <- ggplot(metadata, aes(x = sample, y = index, color = fill)) +
    geom_point(size = 3) +
    labs(title = "Matrix Index per Sample",
         x = "Sample",
         y = "Matrix Index",
         color = "Response") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    scale_y_continuous(limits = c(0.25, 1.75)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plot_or_save(sample_index_plot, deres, plot, save, save_name, "sample_index")
}

boxplot_index <- function(deres, plot=TRUE, save=FALSE, save_name=NULL) {
  matrix_index <- matrisome_index(deres)
  
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  metadata$index <- matrix_index[rownames(metadata)]
  
  if (any(grepl("^group$", colnames(metadata)))) x <- metadata$group
  else x <- metadata$response
  
  index_boxplot <- ggplot(metadata, aes(x=x, y=index, fill=response)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    scale_y_continuous(limits = c(0.25, 1.75)) + 
    labs(title = "Matrix Index by Response Group") + 
    stat_compare_means()
  
  plot_or_save(index_boxplot, deres, plot, save, save_name, "boxplot_index")
}


# ---- 9. Printer Analysis  ----
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
  
  # naba pathways
  naba_pathways <- c(
    "NABA_MATRISOME", "NABA_ECM_REGULATORS", "NABA_ECM_GLYCOPROTEINS",
    "NABA_COLLAGENS", "NABA_PROTEOGLYCANS", "NABA_ECM_AFFILIATED",
    "NABA_CORE_MATRISOME", "NABA_MATRISOME_ASSOCIATED", "NABA_SECRETED_FACTORS"
  )
  
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
  
  
  ### Slide 10+: Pathways
  # GSEA hallmark
  gseaRes <- pathway_enrichment(deres, view_all = FALSE)
  dotplot_pathways(gseaRes, plot = FALSE, save = TRUE, save_name=paste0(file_path, "/dotplot_hallmark.png"))
  
  # GSEA naba
  gseaResNaba <- pathway_enrichment(deres, view_all = TRUE, category = "C2", subset = naba_pathways)
  dotplot_pathways(gseaResNaba, plot = FALSE, save = TRUE, save_name=paste0(file_path, "/dotplot_naba.png"))
  
  # GSEA custom matrisome
  gseaResMatrisome <- custom_pathway_enrichment(deres)
  dotplot_pathways(gseaResMatrisome, plot = FALSE, save = TRUE, save_name=paste0(file_path, "/dotplot_matrisome.png"))
  
  ### Slide 12+: Matrix Index
  sample_index(deres, plot = FALSE, save = TRUE, save_name = paste0(file_path, "/sample_index.png"))
  boxplot_index(deres, plot = FALSE, save = TRUE, save_name = paste0(file_path, "/boxplot_index.png"))
  
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
  slide_files <- c(slide_files, 
                   "dotplot_hallmark.png", "dotplot_naba.png", "dotplot_matrisome.png",
                   "sample_index.png", "boxplot_index.png")
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
  slide_titles <- c(slide_titles, 
                    "Dotplot: Hallmark Pathways", "Dotplot: NABA Pathways", "Dotplot: Matrisome Pathways", 
                    "Matrisome Index: Between Samples", "Matrisome Index: Between Response")
  
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

# combine pdfs
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

# perform printer analyses for all datasets in list
recursive_printer <- function(deres_list) {
  for (i in seq_along(deres_list)) printer_analysis(deres_list[[i]])
}

# ---- 10. Combination Analysis  ----
list_analyse_print <- function(datasets_list) {
  deres_list <- lapply(datasets_list, function(dataset) {
    deres <- response_analysis(dataset$name, dataset$counts, dataset$metadata, dataset$species, "annodbi")
    recursive_printer(deres)
    return(deres)
  })
  combine_pdfs(output_file = "publicdata.pdf")
  return(deres_list)
}

combined_printer <- function(deres_list, pdf_name = "combined_analysis.pdf") {
  gsea_df <- prep_gsea_combined(deres_list)
  plot_gsea_combined(gsea_df, save = TRUE)
  
  genes_combined(deres_list, save=TRUE)
  genes_combined(deres_list, "all", pvalues_cutoff = 2, save=TRUE)
  index_combined(deres_list, save=TRUE)
  
  deconvolute_combined(deres_list, save=TRUE)
  # deconvolute_combined(deres_list, method = "epic", save=TRUE, plot=FALSE)
  deconvolute_combined(deres_list, method = "decosus", save=TRUE, plot=FALSE)
  
  pdf(file = pdf_name, width = 8, height = 6)
  png_files <- list.files(path = "combined_analyses", pattern = "\\.png$", full.names = TRUE)
  for (i in seq_along(png_files)) {
    img <- readPNG(png_files[i])
    grid.newpage()
    grid.draw(rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")))
  }
  dev.off()
  
  print("------ UP ------")
  print(common_genes(deres_list, cutoff = 2))
  print("------ DOWN ------")
  print(common_genes(deres_list, "down", 2))
}

# ---- 10.1 Matrisome Combination Analysis  ----
genes_combined <- function(deres_list, genes = "index", pvalues_cutoff = 0, plot=TRUE, save=FALSE) {
  merged_res_list <- list()
  
  # required files
  matrisome_index <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")
  matrisome_file <- read_csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
  matrisome_file$...4 <- NULL
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      # required metadata
      res <- entry$res
      
      # species check and filter res
      if(entry$species == "human") {
        res <- res %>% filter(gene_name %in% matrisome_file$GeneSymbol)
        res <- merge(res, matrisome_file, 
                     by.x = "gene_name", by.y = "GeneSymbol")
        res$GeneSymbol <- res$gene_name
      }
      else {
        res <- res %>% filter(gene_name %in% matrisome_file$mouse)
        res <- merge(res, matrisome_file,
                     by.x = "gene_name", by.y = "mouse")
        res$mouse <- res$gene_name
      }
      
      # metadata
      res$treatment_status <- sub(".*_(.*)$", "\\1", entry$name)
      res$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
      res$species <- entry$species
      res$tissue <- entry$tissue
      
      merged_res_list[[entry$name]] <- res
    }
  }
  
  combined_res_df <- do.call(rbind, merged_res_list) %>% filter(!is.na(gene_name)) %>%
    mutate(significance = ifelse(is.na(padj), "NA", ifelse(padj < 0.05, "< 0.05", "≥ 0.05")))
  
  if(genes == "index") combined_res_df <- combined_res_df[combined_res_df$GeneSymbol %in% matrisome_index, ]
  
  significant_genes <- c()
  for(i in seq_along(matrisome_file$GeneSymbol)) {
    gene <- matrisome_file$GeneSymbol[i]
    filter_gene <- combined_res_df %>% filter(GeneSymbol == gene & treatment_status == "response")
    table_significance <- table(filter_gene$significance)
    if ( pvalues_cutoff == 0 && 
        ("≥ 0.05" %in% names(table_significance) || "< 0.05" %in% names(table_significance))) {
      significant_genes <- c(significant_genes, gene) 
    }
    else if ("< 0.05" %in% names(table_significance) && table_significance["< 0.05"] >= pvalues_cutoff) {
      significant_genes <- c(significant_genes, gene) 
    }
  }
  filtered_res_df <- combined_res_df[combined_res_df$GeneSymbol %in% significant_genes,]
  
  # colored labels for species
  filtered_res_df$dataset_label <- case_when(
    filtered_res_df$species == "human" ~ paste0("<span style='color:#2ca25f;'>", filtered_res_df$dataset, "</span>"),
    filtered_res_df$species == "mouse" ~ paste0("<span style='color:#de77ae;'>", filtered_res_df$dataset, "</span>"),
    TRUE ~ filtered_res_df$dataset
  )
  
  # order dataset by tissue
  filtered_res_df$dataset_label <- factor(
    filtered_res_df$dataset_label,
    levels = filtered_res_df %>%
      distinct(dataset_label, tissue) %>%
      arrange(tissue, dataset_label) %>%
      pull(dataset_label)
  ) %>% paste0(filtered_res_df$tissue, " | ", .)
  
  #plot
  gene_plot <- ggplot(filtered_res_df, aes(x = GeneSymbol, y = dataset_label)) +
    # dotplot
    geom_point(aes(size = significance, color = log2FoldChange)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5, "NA" = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # nested facet plot
    ggh4x::facet_nested(rows = vars(treatment_status), cols = vars(Division, Category), scales = "free", space = "free") +
    
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_markdown(size = 8),
      strip.text = element_text(size = 10)
    )
  
  if(plot) plot(gene_plot)
  if(save) {
    if(genes == "index") ggsave("combined_analyses/index_genes_combined.png", gene_plot, scale = 3)
    else ggsave("combined_analyses/matrisome_genes_combined.png", gene_plot, scale=3)
  }
}

index_combined <- function(deres_list, plot=TRUE, save=FALSE) {
  matrisome_index_list <- list()
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      if(grepl("_response$", entry$name) | grepl("_merged$", entry$name)) {
        matrix_index <- matrisome_index(entry)
        
        metadata <- entry$dds %>% colData() %>% as.data.frame()
        metadata$index <- matrix_index[rownames(metadata)]
        metadata$response <- tolower(metadata$response)
        metadata$tissue <- entry$tissue
        metadata$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
        
        metadata <- metadata %>% select(index, dataset, response, tissue)
        
        matrisome_index_list[[entry$name]] <- metadata
      }
    }
  }
  
  matrisome_index_df <- do.call(rbind, matrisome_index_list)
  
  # order dataset by tissue
  matrisome_index_df$dataset <- factor(
    matrisome_index_df$dataset,
    levels = matrisome_index_df %>%
      distinct(dataset, tissue) %>%
      arrange(tissue, dataset) %>%
      pull(dataset)
  ) %>% paste0(matrisome_index_df$tissue, " | ", .)
  
  matrisome_index_plot <- ggplot(matrisome_index_df, aes(x=dataset, y=index, fill=response)) +
    geom_boxplot() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    # scale_y_continuous(limits = c(0.7, 1.75)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = NULL, y = "Matrisome Index")
  
  
  
  if(plot) plot(matrisome_index_plot)
  if(save) ggsave("combined_analyses/matrisome_index_combined.png", matrisome_index_plot, scale = 3)
}


# ---- 10.2 Pathways Combination Analysis  ----

# function to prep gsea (to avoid running gsea analyses multiple times)
prep_gsea_combined <- function(deres_list, type = "hallmark") {
  merged_gsea_list <- list()
  
  # required pathways
  naba_pathways <- c(
    "NABA_MATRISOME", "NABA_ECM_REGULATORS", "NABA_ECM_GLYCOPROTEINS",
    "NABA_COLLAGENS", "NABA_PROTEOGLYCANS", "NABA_ECM_AFFILIATED",
    "NABA_CORE_MATRISOME", "NABA_MATRISOME_ASSOCIATED", "NABA_SECRETED_FACTORS"
  )
  
  hallmark_pathways <- data.frame(
    ID = c(
      "E2F_TARGETS", "G2M_CHECKPOINT", "MITOTIC_SPINDLE", "MYC_TARGETS_V1", "MYC_TARGETS_V2",
      "ALLOGRAFT_REJECTION", "COMPLEMENT", "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "TNFA_SIGNALING_VIA_NFKB",
      "HEDGEHOG_SIGNALING", "KRAS_SIGNALING_DN", "KRAS_SIGNALING_UP", "MTORC1_SIGNALING", "PI3K_AKT_MTOR_SIGNALING", "TGF_BETA_SIGNALING", "WNT_BETA_CATENIN_SIGNALING",
      "FATTY_ACID_METABOLISM", "GLYCOLYSIS", "HYPOXIA", "OXIDATIVE_PHOSPHORYLATION", "REACTIVE_OXYGEN_SPECIES_PATHWAY",
      "APOPTOSIS", "DNA_REPAIR", "P53_PATHWAY",
      "ANGIOGENESIS", "EPITHELIAL_MESENCHYMAL_TRANSITION"
    ),
    category = c(
      rep("proliferation", 5),
      rep("immune", 8),
      rep("signaling", 7),
      rep("metabolism", 5),
      rep("cell_death", 3),
      rep("other", 2)
    )
  )
  
  matrisome_pathways <- data.frame(
    ID = c("ECM Glycoproteins", "Collagens", "Proteoglycans",
            "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"),
    category = c(rep("Core", 3), rep("Assoc", 3))
  )
  
  # loop through deres
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
        entry <- deres_list[[i]][[j]]
        
        # prep hallmark pathways
        hallmark_enrich <- pathway_enrichment(
          entry, subset = paste0("HALLMARK_", hallmark_pathways$ID)
          )@result
        hallmark_enrich$facet <- "Hallmark"
        hallmark_enrich$ID <- sub("HALLMARK_", "", hallmark_enrich$ID)
        hallmark_enrich <- left_join(hallmark_enrich, hallmark_pathways, by="ID")
        
        # prep naba pathways
        naba_enrich <- pathway_enrichment(entry, view_all = TRUE, category = "C2", subset = naba_pathways)@result
        naba_enrich$facet <- "NABA"
        naba_enrich$ID <- sub("NABA_", "", naba_enrich$ID)
        naba_enrich$category <- "pathways"
        
        # prep matrisome pathways
        custom_enrich <- custom_pathway_enrichment(entry, view_all = TRUE)@result
        custom_enrich$facet <- "Matrisome"
        custom_enrich <- left_join(custom_enrich, matrisome_pathways, by="ID")
        
        # combine
        gsea <- rbind(hallmark_enrich, naba_enrich, custom_enrich)
        
        # metadata
        gsea$treatment_status <- sub(".*_(.*)$", "\\1", entry$name)
        gsea$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
        gsea$species <- entry$species
        gsea$tissue <- entry$tissue
        
        merged_gsea_list[[entry$name]] <- gsea
      }
    }
  
  # combine all
  combined_gsea_df <- do.call(rbind, merged_gsea_list) %>% 
    mutate(significance = ifelse(p.adjust < 0.05, "< 0.05", "≥ 0.05"))
  
  # colored labels for species
  combined_gsea_df$dataset_label <- case_when(
    combined_gsea_df$species == "human" ~ paste0("<span style='color:#2ca25f;'>", combined_gsea_df$dataset, "</span>"),
    combined_gsea_df$species == "mouse" ~ paste0("<span style='color:#de77ae;'>", combined_gsea_df$dataset, "</span>"),
    TRUE ~ combined_gsea_df$dataset
  )
  
  # order dataset by tissue
  combined_gsea_df$dataset_label <- factor(
    combined_gsea_df$dataset_label,
    levels = combined_gsea_df %>%
      distinct(dataset_label, tissue) %>%
      arrange(tissue, dataset_label) %>%
      pull(dataset_label)
  ) %>% paste0(combined_gsea_df$tissue, " | ", .)
  
  # order ID by facet/category
  combined_gsea_df$ID <- factor(combined_gsea_df$ID, 
                                levels = unique(combined_gsea_df$ID))

  return(combined_gsea_df)
}

plot_gsea_combined <- function(combined_gsea_df, plot=TRUE, save=FALSE) {
  gsea_combined <- ggplot(combined_gsea_df, aes(x = ID, y = dataset_label)) +
    # dotplot
    geom_point(aes(size = significance, color = NES)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # nested facet plot
    ggh4x::facet_nested(rows = vars(treatment_status), cols = vars(facet, category), scales = "free", space = "free") +
    
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_markdown(size = 8),
      strip.text = element_text(size = 10)
    ) + labs(x = NULL, y = NULL)

  
  if(plot) plot(gsea_combined)
  if(save) ggsave("combined_analyses/gsea_combined.png", gsea_combined, scale = 3)
}

# ---- 10.3 Common Genes  ----
# find commonly upregulated and downregulated genes
common_genes <- function(deres_list, direction="up", cutoff=5) {
  genes <- c("s")
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      if(grepl("_response$", entry$name)) {
        res_genes <- view_res_table(entry, direction)$gene_name
        if(entry$species == "mouse") res_genes <- to_mouse_case(res_genes, "tohuman", "annodbi")
        res_genes <- res_genes[!is.na(res_genes)]
        genes <- c(genes, res_genes)
      }
    }
  }
  gene_counts <- table(genes)
  common <- gene_counts[gene_counts >= cutoff]
  return(sort(common, decreasing = TRUE))
}

# ---- 10.3 Deconvolution Combined  ----
deconvolute_combined <- function(deres_list, method = "xcell", plot=TRUE, save=FALSE) {
  merged_deconvolution_list <- list()
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      treatment_status <- sub(".*_(.*)$", "\\1", entry$name)
      
      if (entry$species == "human" & (treatment_status == "response" || treatment_status == "merged")) {
        # Try to run deconvolute_significant, and skip on error
        tryCatch({
          deconv <- deconvolute_significant(entry, method = method)
          deconv$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
          deconv$tissue <- entry$tissue
          merged_deconvolution_list[[entry$name]] <- deconv
        }, error = function(e) {
          message(sprintf("Skipping %s due to error: %s", entry$name, e$message))
        })
      }
    }
  }
  
  # Combine all data and continue as normal
  combined_deconv_df <- do.call(rbind, merged_deconvolution_list) %>%
    mutate(log_p = -log10(p_adj)) %>%
    mutate(significance = ifelse(p_adj < 0.05, "< 0.05", "≥ 0.05"))
  
  # Order dataset by tissue
  combined_deconv_df$dataset <- factor(
    combined_deconv_df$dataset,
    levels = combined_deconv_df %>%
      distinct(dataset, tissue) %>%
      arrange(tissue, dataset) %>%
      pull(dataset)
  ) %>% paste0(combined_deconv_df$tissue, " | ", .)
  
  # Plot
  combined_deconv_plot <- ggplot(combined_deconv_df, aes(x = dataset, y = cell)) +
    geom_point(aes(size = significance, color=log2_fc)) + 
    scale_size_manual(values = c("≥ 0.05" = 2, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  if(plot) plot(combined_deconv_plot)
  if(save) ggsave( paste0("combined_analyses/", method, "_deconv_combined.png"), combined_deconv_plot, scale = 3)
}

# ---- 11. Merger Analysis ----

merge_datasets <- function(deres_list, tissue="pancreas") {
  required_deres <- list()
  # get required datasets
  for (i in seq_along(deres_list)) {
    deres <- deres_list[[i]][[1]]
    
    if (deres$species == "human") {
      if (tissue == "all" || deres$tissue == tissue) required_deres[[deres$name]] <- deres
    }
  }
  
  # prep dataframes
  metadata <- data.frame()
  counts_matrix <- data.frame()
  # loop through required and merge
  for (i in seq_along(required_deres)) {
    deres <- required_deres[[i]]
    required_met <- deres$dds %>% colData() %>% as.data.frame() %>% select(response, batch)
    metadata <- rbind(metadata, required_met)
    
    required_counts <- deres$dds %>% assay() %>% reverse_human_annotate()
    required_counts$gene_id <- rownames(required_counts)
    
    if (nrow(counts_matrix) == 0) {
      counts_matrix <- required_counts
    } else {
      counts_matrix <- full_join(counts_matrix, required_counts, by = "gene_id")
    }
  }
  
  # set rownames, set na to zero
  rownames(counts_matrix) <- counts_matrix$gene_id
  counts_matrix$gene_id <- NULL
  counts_matrix[is.na(counts_matrix)] <- 0
  
  dds <- DESeqDataSetFromMatrix(counts_matrix, metadata, ~batch+response)
  dds$response <- relevel(dds$response, ref = "nonresponder") # relevel dds
  
  return(list(metadata=metadata, counts=counts_matrix, dds=dds))
}

# ---- 11.1 Merged Differential Expression  ----
merged_de <- function(deres_list, tissue="pancreas") {
  merged_dataset <- merge_datasets(deres_list, tissue)
  deres <- run_deseq(merged_dataset$dds, "response_responder_vs_nonresponder")
  
  # annotate
  deres$res <- deres$res %>% as.data.frame() %>% annotate_df("annodbi", "human")
  
  # assign metadata
  deres$name <- paste0(tissue,"_merged")
  deres$tissue <- tissue
  deres$species <- "human"
  
  return(deres)
}

# ---- 11.2 Merged Machine Learning  ----

## Functions to prepare datasets 


# batch correct vsd
combat_correction <- function(deres_list, tissue="pancreas") {
  merged_dataset <- merge_datasets(deres_list, tissue)
  
  # Normalise and Log transform
  dds <- estimateSizeFactors(merged_dataset$dds)
  vsd <- vst(dds, blind = FALSE)
  log_expr <- assay(vsd)
  
  # Batch Correction
  mod <- model.matrix(~ response, data = merged_dataset$metadata)
  combat_expr <- ComBat(dat = log_expr, batch = merged_dataset$metadata$batch, mod = mod)
  
  return(list(expr=combat_expr, metadata=merged_dataset$metadata))
}

# annotate vsd
merged_filter_matrisome <- function(expr) {
  matrisome_file <- read_csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
  matrisome_file$...4 <- NULL
  
  expr <- expr %>% as.data.frame() %>% annotate_df("annodbi", "human") %>%
    filter(gene_name %in% matrisome_file$GeneSymbol)
  expr <- expr[!duplicated(expr$gene_name), ]
  
  rownames(expr) <- expr$gene_name
  expr$gene_name <- NULL
  
  return(expr)
}

# transpose vsd for ml packages
prepare_ml_data <- function(deres_list, use_matrisome=FALSE, tissue="pancreas") {
  cor_data <- combat_correction(deres_list, tissue)
  
  if(isTRUE(use_matrisome)) {cor_data$expr <- cor_data$expr %>% merged_filter_matrisome()}
  
  df <- as.data.frame(t(cor_data$expr))  # samples as rows
  df$response <- as.factor(cor_data$metadata$response)
  return(df)
}

pca_merged <- function(deres_list, tissue="pancreas") {
  # Perform PCA
  ml_df <- prepare_ml_data(deres_list, tissue)
  pca_res <- PCA(ml_df[ , !colnames(ml_df) %in% "response"], graph = FALSE)
  
  # Visualise
  fviz_pca_ind(pca_res,
               geom.ind = "point",
               habillage = ml_df$response,
               addEllipses = TRUE,
               palette = "jco")
}

# ---- 11.2.1 ML Models ----
run_ml_models <- function(df, method = "rf", tuneLength = 5, seed = 42) {
  set.seed(seed)
  
  # Split data
  trainIndex <- createDataPartition(df$response, p = 0.7, list = FALSE)
  train <- df[trainIndex, ]
  test <- df[-trainIndex, ]
  
  # Train model
  trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                         classProbs = TRUE, summaryFunction = twoClassSummary)
  
  model <- train(response ~ ., data = train, method = method,
                 metric = "ROC", trControl = trctrl, tuneLength = tuneLength)
  
  # Predict on test set
  predictions <- predict(model, newdata = test)
  prob <- predict(model, newdata = test, type = "prob")
  
  # Evaluate
  confusion <- confusionMatrix(predictions, test$response)
  roc_obj <- roc(response = test$response, predictor = prob[, "responder"])
  
  list(model = model,
       confusion = confusion,
       roc = roc_obj,
       auc = auc(roc_obj))
}

compare_ml_models <- function(df, models = c("rf", "svmRadial", "glmnet"),
                              tuneLength = 5, seed = 42) {
  set.seed(seed)
  
  trainIndex <- createDataPartition(df$response, p = 0.7, list = FALSE)
  train <- df[trainIndex, ]
  test <- df[-trainIndex, ]
  
  trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                         classProbs = TRUE, summaryFunction = twoClassSummary,
                         savePredictions = "final")
  
  results <- list()
  
  for (model in models) {
    message("Training model: ", model)
    
    fit <- train(response ~ ., data = train, method = model,
                 trControl = trctrl, metric = "ROC", tuneLength = tuneLength)
    
    pred <- predict(fit, newdata = test)
    prob <- predict(fit, newdata = test, type = "prob")
    
    roc_obj <- roc(response = test$response, predictor = prob[, "responder"])
    
    results[[model]] <- list(model = fit,
                             confusion = confusionMatrix(pred, test$response),
                             auc = auc(roc_obj),
                             roc = roc_obj)
  }
  
  return(results)
}


plot_model_aucs <- function(results) {
  aucs <- sapply(results, function(x) x$auc)
  barplot(aucs, beside = TRUE, col = "skyblue",
          main = "AUCs of ML Models", ylab = "AUC", ylim = c(0, 1))
}

run_rfe <- function(df, method = "rf", seed = 42) {
  set.seed(seed)
  
  ctrl <- rfeControl(functions = rfFuncs,
                     method = "repeatedcv",
                     number = 5, repeats = 3)
  
  predictors <- df[ , !names(df) %in% "response"]
  response <- df$response
  
  rfe_result <- rfe(predictors, response,
                    sizes = c(5, 10, 20, 50, 100),
                    rfeControl = ctrl)
  
  return(rfe_result)
}

extract_model_features <- function(model_results, top_n = 20, plot = TRUE) {
  all_features <- list()
  
  for (model_name in names(model_results)) {
    model_obj <- model_results[[model_name]]$model
    
    # Try extracting importance safely
    importance_obj <- try(varImp(model_obj), silent = TRUE)
    if (inherits(importance_obj, "try-error")) {
      message("Skipping model ", model_name, ": varImp not available.")
      next
    }
    
    importance_df <- importance_obj$importance
    importance_df$Feature <- rownames(importance_df)
    
    # Handle multi-class or multiple metrics
    if (!"Overall" %in% colnames(importance_df)) {
      # Take mean importance across columns
      importance_df$Overall <- rowMeans(importance_df[ , sapply(importance_df, is.numeric)], na.rm = TRUE)
    }
    
    top_df <- importance_df %>%
      arrange(desc(Overall)) %>%
      head(top_n)
    
    top_df$Model <- model_name
    all_features[[model_name]] <- top_df
    
    if (plot) {
      ggplot(top_df, aes(x = reorder(Feature, Overall), y = Overall)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        labs(title = paste("Top", top_n, "Features -", model_name),
             x = "Gene", y = "Importance") +
        theme_minimal() -> p
      print(p)
    }
  }
  
  return(all_features)
}

# ---- 11.2.2 ML Gene Signatures ----

build_gene_signature <- function(df, seed = 42) {
  set.seed(seed)
  
  x <- as.matrix(df[ , !colnames(df) %in% "response"])
  y <- as.factor(df$response)
  y_bin <- ifelse(y == "responder", 1, 0)
  
  # Fit LASSO model (logistic regression with L1 penalty)
  cvfit <- cv.glmnet(x, y_bin, family = "binomial", alpha = 1)
  
  # Extract coefficients at lambda.1se (simpler model)
  coefs <- coef(cvfit, s = "lambda.min")
  coef_df <- as.data.frame(as.matrix(coefs))
  coef_df$Gene <- rownames(coef_df)
  colnames(coef_df)[1] <- "Weight"
  
  # Filter non-zero coefficients (excluding intercept)
  gene_signature <- coef_df %>%
    filter(Weight != 0 & Gene != "(Intercept)")
  
  return(list(signature = gene_signature,
              model = cvfit))
}

prepare_counts_ml_df <- function(counts) {
  col_data <- data.frame(row.names = colnames(counts),
                         dummy = factor(rep("A", ncol(counts))))
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts ,
                                colData = col_data,
                                design = ~1)
  
  # Normalize + transform (VST recommended)
  dds <- estimateSizeFactors(dds)
  vst_expr <- assay(vst(dds, blind = TRUE)) %>% merged_filter_matrisome()
  
  # Transpose to match expected input (samples as rows)
  vst_df <- as.data.frame(t(vst_expr))
  
  return(vst_df)
}

apply_gene_signature <- function(expr_df, gene_signature) {
  genes <- gene_signature$Gene
  weights <- gene_signature$Weight
  
  expr_sub <- expr_df[ , genes, drop = FALSE]
  
  # If gene is missing in new data, fill with 0
  missing_genes <- setdiff(genes, colnames(expr_df))
  if (length(missing_genes) > 0) {
    expr_sub[ , missing_genes] <- 0
  }
  
  # Align columns
  expr_sub <- expr_sub[ , genes]
  
  # Compute score = weighted sum
  scores <- as.matrix(expr_sub) %*% as.matrix(weights)
  
  return(data.frame(Sample = rownames(expr_df), Score = as.numeric(scores)))
}

evaluate_signature_performance <- function(score_df, label_col = "response", score_col = "Score", plot = TRUE) {
  library(pROC)
  library(caret)
  library(ggplot2)
  library(dplyr)
  
  # Drop rows with missing response or score
  score_df <- score_df %>% filter(!is.na(.data[[label_col]]), !is.na(.data[[score_col]]))
  
  response <- factor(score_df[[label_col]], levels = c("nonresponder", "responder"))
  score <- score_df[[score_col]]
  
  # ROC + AUC
  roc_obj <- roc(response = response, predictor = score)
  auc_val <- auc(roc_obj)
  
  # Optimal threshold
  coords_opt <- coords(roc_obj, "best", best.method = "youden", transpose = FALSE)
  opt_threshold <- coords_opt["threshold"]
  
  # Predictions based on threshold
  predicted <- ifelse(score > opt_threshold, "responder", "nonresponder")
  predicted <- factor(predicted, levels = c("nonresponder", "responder"))
  
  # Ensure matching length
  valid_idx <- complete.cases(response, predicted)
  response <- response[valid_idx]
  predicted <- predicted[valid_idx]
  
  cm <- confusionMatrix(predicted, response)
  
  metrics <- data.frame(
    Accuracy = cm$overall["Accuracy"],
    Precision = cm$byClass["Pos Pred Value"],
    Recall = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    F1 = 2 * ((cm$byClass["Sensitivity"] * cm$byClass["Pos Pred Value"]) /
                (cm$byClass["Sensitivity"] + cm$byClass["Pos Pred Value"]))
  )
  
  if (plot) {
    # ROC plot
    plot(roc_obj, col = "blue", lwd = 2, main = "ROC Curve - Gene Signature")
    legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "blue", lwd = 2)
    
    # Confusion plot
    score_df$Predicted <- predicted
    ggplot(score_df, aes(x = .data[[label_col]], fill = Predicted)) +
      geom_bar(position = "dodge") +
      labs(title = "True vs Predicted Classes",
           x = "True Label", y = "Count") +
      theme_minimal() +
      scale_fill_manual(values = c("firebrick", "steelblue")) -> p
    print(p)
  }
  
  return(list(
    auc = auc_val,
    threshold = opt_threshold,
    confusion_matrix = cm,
    metrics = metrics,
    roc = roc_obj,
    predicted = predicted
  ))
}

# ---- 12. Other Useful Functions  ----


check_valid_data <- function(counts, metadata) {
  print(all(colnames(counts) == rownames(metadata)))
}

save_geo_data <- function(gse, counts, metadata) {
  write.csv(counts, paste0("bulk/", gse, "/counts.csv"), row.names = TRUE)
  write.csv(metadata, paste0("bulk/", gse, "/metadata.csv"), row.names = TRUE)
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

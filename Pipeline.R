# RNA Seq Response Analysis Pipeline

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
library(GSVA)

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
library(ComplexHeatmap)

# printer analysis 
library(grid)
library(gridExtra)
library(png)
library(qpdf)

# ML packages
library(caret)
library(pROC)
library(FactoMineR)
library(factoextra)
library(glmnet)

library(BiocParallel)


# ----- VARIABLES ALTER AS NEEDED  ----
matrisome_file <- read.csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
immuno_file <- read.csv("~/Repos/RNASeq Pipeline/InnateDB_genes.csv")
matrisome_index_genes <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")

positive_genes <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1")
negative_genes <- c("AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")

naba_pathways <- c(
  "NABA_MATRISOME", "NABA_ECM_REGULATORS", "NABA_ECM_GLYCOPROTEINS",
  "NABA_COLLAGENS", "NABA_PROTEOGLYCANS", "NABA_ECM_AFFILIATED",
  "NABA_CORE_MATRISOME", "NABA_MATRISOME_ASSOCIATED", "NABA_SECRETED_FACTORS"
)

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
    
    response_treatment_validity <- metadata %>%
      group_by(treatment = tolower(treatment)) %>%
      summarise(all_same = length(unique(tolower(response))) != 1) %>%
      deframe()
    
    # dds for UNTREATED response vs nonresponse
    if(response_treatment_validity["untreated"]) dds_list$untreated <- prepare_dds(counts, metadata, ~group, "nonresponderuntreated")
    
    # dds for TREATED response vs nonresponse
    if(response_treatment_validity["treated"]) dds_list$treated <- prepare_dds(counts, metadata, ~group, "nonrespondertreated")
    
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
annotate_df <- function(df, method="annodbi", species="human") {
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

annotate_df_remove_dup <- function(df, method="annodbi", species="human") {
  # annoate df using existing function
  df <- df %>% as.data.frame()
  annotated_df <- annotate_df(df, method, species)
  annotated_df <- annotated_df[!is.na(annotated_df$gene_name), ] # remove na
  
  expr_cols <- sapply(annotated_df, is.numeric)
  expr_cols["gene_name"] <- FALSE
  
  # take max of row sums
  annotated_df$expr_sum <- rowSums(annotated_df[, expr_cols], na.rm = TRUE) # row sums of numeric columns
  annotated_df <- annotated_df[order(-annotated_df$expr_sum), ] 
  annotated_df <- annotated_df[!duplicated(annotated_df$gene_name), ]
  
  # assign rownames
  rownames(annotated_df) <- annotated_df$gene_name
  annotated_df$gene_name <- NULL
  annotated_df$expr_sum <- NULL
  
  return(annotated_df)
}

reverse_human_annotate <- function(counts, method="annodbi") {
  counts_df <- as.data.frame(counts)
  ids <- rownames(counts_df)
  annotation_type <- annotation_status(ids)
  
  # handle Ensembl case - already good
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
  
  
  counts_df$original <- ids
  final_df <- dplyr::left_join(counts_df, df_merged, by = "original")
  
  ensembl_col <- ifelse("ENSEMBL" %in% colnames(final_df), "ENSEMBL", "ensembl_gene_id")
  
  
  # filter out missing Ensembl ID
  final_df <- final_df[!is.na(final_df[[ensembl_col]]), ]
  # keep only the first occurrence of each Ensembl ID
  final_df <- final_df[!duplicated(final_df[[ensembl_col]]), ]
  # assign Ensembl ID as rownames
  rownames(final_df) <- final_df[[ensembl_col]]
  # drop helper columns
  final_df[[ensembl_col]] <- NULL
  final_df$original <- NULL
  
  return(final_df)
}

# To mouse case
# Required if using homologene mapping
to_mouse_case <- function(symbols, method = "tomouse", source = "biomart") {
  
  # cheating method
  if (method == "case") {
    return(sapply(symbols, function(g) {
      if (is.na(g) || nchar(g) == 0) return("")
      first <- substr(g, 1, 1)
      rest <- tolower(substr(g, 2, nchar(g)))
      paste0(first, rest)
    }))
  }
  
  # filter out incorrect method
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
      
      # keep first valid mapping per symbol
      mouseg <- sapply(mouseg_list, function(x) if (length(x) > 0) x[1] else NA)
      
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
    if("untreated" %in% names(dds_list)) deres_list$untreated <- run_deseq(dds_list$untreated, "group_responderuntreated_vs_nonresponderuntreated")
    if("treated" %in% names(dds_list)) deres_list$treated <- run_deseq(dds_list$treated, "group_respondertreated_vs_nonrespondertreated")
  }
  message("Finished running DESeq2")
  
  names_deres_list <- names(deres_list)
  
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
      # assign pancreas if pancreas included despite not being the majority
      tissues <- colData(deres$dds)$tissue
      if("pancreas" %in% tissues) tissue <- "pancreas"
      else tissue <- names(sort(table(tissues), decreasing = TRUE))[1]
    }
    deres$tissue <- tissue
    
    return(deres)
  })
  names(deres_list) <- names_deres_list
  
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
  
  heatmap <- pheatmap::pheatmap(scale.dat,
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
prep_pathway_enrichment <- function(deres) {
  # get DESeq2 results
  res <- deres$res
  # remove NAs
  res <- res[!is.na(res$gene_name) & !is.na(res$stat), ]
  # remove duplicates (keep highest stat (Wald statistic) per gene)
  res <- res[!duplicated(res$gene_name), ]
  
  # create ranked gene list
  geneList <- res$stat
  names(geneList) <- res$gene_name
  geneList <- sort(geneList, decreasing = TRUE)
  
  return(geneList)
}

pathway_enrichment <- function(deres, view_all = TRUE, category = "H", subset = NULL) {
  species <- deres$species
  
  geneList <- prep_pathway_enrichment(deres)
  
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

matrisome_pathway_enrichment <- function(deres, view_all = TRUE) {
  species <- deres$species
  
  geneList <- prep_pathway_enrichment(deres)
  
  # custom pathway
  matrisome_proteins <- matrisome_file
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

# custom_human_pathway = df w/ 1 column term, 1 column gene
custom_pathway_enrichment <- function(deres, custom_human_pathway, view_all=TRUE) {
  species <- deres$species
  
  geneList <- prep_pathway_enrichment(deres)
  
  # set p-value cutoff depending on view_all
  pcut <- ifelse(view_all, 2, 0.05)
  
  # run GSEA
  gseaRes <- GSEA(geneList, TERM2GENE = custom_human_pathway, pvalueCutoff = pcut)
  
  return(gseaRes)
}

dotplot_pathways <- function(gseaRes, plot=TRUE, save=FALSE, save_name=NULL) { 
  pathways_dotplot <- clusterProfiler::dotplot(gseaRes, x = "GeneRatio",
                              color = "NES", size = "p.adjust")
  
  plot_or_save(pathways_dotplot, deres, plot, save, save_name, "pathways")
}

ssGSEA_score <- function(counts, metadata, gene_list, plot=TRUE) {
  gene_sets <- gene_list %>% group_by(term) %>% 
    summarise(genes = list(unique(gene))) %>% deframe()
  
  expr_mat <- counts %>% prepare_new_data() %>% as.matrix()
  
  ssgsea <- ssgseaParam(
    expr_mat, gene_sets,
    minSize = 1, maxSize = Inf,
    alpha = 0.25, normalize = TRUE
  )
  
  ssgsea_scores <- gsva(ssgsea, verbose = TRUE)
  
  ssgsea_df <- as.data.frame(t(ssgsea_scores)) %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "model", values_to = "ssgsea_score")
  
  metadata$sample <- rownames(metadata)
  ssgsea_df <- ssgsea_df %>%
    left_join(metadata, by = "sample")
  
  score_plot <- ggplot(ssgsea_df, aes(x = response, y = ssgsea_score, fill = response)) +
    geom_boxplot() +
    facet_wrap(~ model, scales = "free_y") +
    theme_minimal() +
    labs(title = "ssGSEA scores") + 
    stat_compare_means()
  
  if(plot) plot(score_plot)
  
  return(ssgsea_scores)
}

bagaev_clustering <- function(merged_dataset) {
  set.seed(42)
  
  # Load and Format Gene Sets
  bagaev_sets <- read.csv("~/Repos/RNASeq Pipeline/bagaev.csv", skip = 1) %>%
    rename(term = Gene.signature, gene = Gene)
  
  # Batch correct
  combat_corr <- combat_correction(merged_dataset)
  
  # ssGSEA scoring
  ssgsea_bagaev_score <- ssGSEA_score(
    combat_corr$expr,
    combat_corr$metadata,
    bagaev_sets,
    plot = FALSE
  )
  
  # prepare plot
  ## row groups
  row_groups <- setNames(bagaev_sets$Signature.Cluster, bagaev_sets$term)
  matched_groups <- row_groups[rownames(ssgsea_bagaev_score)]
  
  desired_row_group_order <- c(
    "Angiogenesis/Fibroblasts",
    "Pro-tumor/Immune infiltrate",
    "Anti-tumor/Immune infiltrate",
    "EMT signature/proliferation rate"
  )
  
  row_group_factor <- factor(
    matched_groups,
    levels = desired_row_group_order
  )
  
  ## colgroups
  cluster_metadata <- combat_corr$metadata
  cluster_metadata$response <- factor(cluster_metadata$response)
  
  ann_colors <- list(
    response = c("responder" = "forestgreen", "nonresponder" = "firebrick")
  )
  
  top_anno <- HeatmapAnnotation(
    response = cluster_metadata$response,
    batch = cluster_metadata$batch,
    col = ann_colors,
    annotation_legend_param = list(
      response = list(title = "Response")
    )
  )
  
  ## heatmap
  ssgsea_scaled <- t(scale(t(ssgsea_bagaev_score)))
  
  ssgsea_scaled_clean <- ssgsea_scaled
  for (i in 1:nrow(ssgsea_scaled_clean)) {
    if (anyNA(ssgsea_scaled_clean[i, ])) {
      ssq_row <- ssgsea_scaled_clean[i, ]
      ssq_row[is.na(ssq_row)] <- mean(ssq_row, na.rm = TRUE)
      ssgsea_scaled_clean[i, ] <- ssq_row
    }
  }
  
  Heatmap(
    ssgsea_scaled_clean,
    name = "Z-score",
    top_annotation = top_anno,
    row_split = row_group_factor,
    cluster_columns = TRUE,
    column_km = 4,
    cluster_rows = FALSE,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = TRUE,
    heatmap_legend_param = list(title = "Z-score", color_bar = "continuous")
  )
}

# ---- 7. Bulk Deconvolution  ----
raw_to_tpm <- function(raw_counts) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gene_lengths <- getBM(attributes = c("hgnc_symbol", "transcript_length"),
                        mart = mart)
  
  # get max length per gene
  gene_lengths_summary <- gene_lengths %>%
    filter(hgnc_symbol != "") %>%
    group_by(hgnc_symbol) %>%
    summarise(length_kb = max(transcript_length, na.rm = TRUE) / 1000)
  
  # get annotated raw counts
  raw_annotated_counts <- raw_counts %>% as.data.frame() %>% annotate_df_remove_dup()
  raw_counts_df <- raw_annotated_counts %>%
    as.data.frame() %>%
    rownames_to_column("hgnc_symbol")
  
  counts_with_lengths <- raw_counts_df %>%
    inner_join(gene_lengths_summary, by = "hgnc_symbol")
  
  # extract length and counts
  gene_lengths <- counts_with_lengths$length_kb
  count_matrix <- counts_with_lengths %>% select(-hgnc_symbol, -length_kb)
  
  # calculate RPK
  rpk <- sweep(count_matrix, 1, gene_lengths, FUN = "/")
  
  # scaling factors per sample
  scaling_factors <- colSums(rpk)
  
  # TPM = RPK / scaling_factor * 1e6
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  
  # reattach gene_names
  tpm <- tpm %>%
    as.data.frame() %>%
    mutate(hgnc_symbol = counts_with_lengths$hgnc_symbol) %>%
    column_to_rownames("hgnc_symbol")
  return(tpm)
}

deconvolute <- function(expr, method="xcell") {
  # get scores and metadata
  if(method == "xcell") {
    scores <- xCell::xCellAnalysis(expr) %>% as.data.frame()
  } else if(method == "epic") {
    scores <- EPIC::EPIC(bulk = expr) %>% .$cellFractions %>% t() %>% as.data.frame()
  } else if(method == "mcp") {
    scores <- MCPcounter::MCPcounter.estimate(
      expr, featuresType = "HUGO_symbols"
    ) %>% as.data.frame()
  } else if(method == "decosus") {
    expr_df <- as.data.frame(expr)
    # add required hgnc_symbols column
    expr_df$hgnc_symbols <- rownames(expr_df)
    decosus_result <- Decosus::cosDeco(x=expr_df, rnaseq = T, plot = TRUE, ext = FALSE)
    est_Decosus_cells <- decosus_result$main_cells
    
    scores <- est_Decosus_cells %>% as.data.frame() %>%
      mutate(across(everything(), ~ ifelse(is.nan(.), 0, .)))
  }
  
  return(scores)
}


deconvolute_significant <- function(deres, method="xcell") {
  expr <- raw_to_tpm(deres$dds %>% assay())
  
  # get scores and metadata
  scores <- deconvolute(expr, method)
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  
  # get responder status and filter
  responders <- metadata %>% filter(response == "responder") %>% rownames()
  scores_responders <- scores %>% select(responders)
  
  # get nonresponders
  non_responders <- metadata %>% filter(response != "responder") %>% rownames()
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
  
  test_statistics$overall_mean <- rowMeans(scores)
  
  return(test_statistics)
}

volcano_deconvolution <- function(deres, method="xcell", plot=TRUE, save=FALSE, save_name=NULL) {
  test_stats <- deconvolute_significant(deres, method)
  
  volcano_plot <- EnhancedVolcano(test_stats, lab=rownames(test_stats),
                                  x='log2_fc', y='p_adj',
                                  title=paste0(method, " ", deres$title), pCutoff=0.05)
  
  plot_or_save(volcano_plot, deres, plot, save, save_name, "volcano")
}


cluster_deconvolution_individual <- function(deres, method="xcell", plot=TRUE, save=FALSE, save_name=NULL) {
  set.seed(42)
  
  treatment_status <- sub(".*_(.*)$", "\\1", deres$name)
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  # filter untreated
  if (treatment_status != "response" & treatment_status != "merged") metadata <- metadata %>% filter(tolower(treatment) == treatment_status)

  counts <- deres$dds %>% assay() %>% as.data.frame() %>% select(rownames(metadata))
  expr <- counts %>% raw_to_tpm() %>% as.matrix()

  deconvolution <- deconvolute(expr, method)
  
  meta_cols <- c("response")
  if ("treatment" %in% colnames(metadata)) {
    meta_cols <- c(meta_cols, "treatment")
  }
  
  df <- as.data.frame(metadata[, meta_cols, drop = FALSE])
  
  heatmap <- pheatmap::pheatmap(deconvolution %>% as.matrix,
                                cutree_cols = 3, 
                                cluster_rows = FALSE,
                                show_rownames = TRUE,
                                show_colnames = FALSE,
                                annotation_col = df,
                                silent=!plot)
  
  if(is.null(save_name)) save_name <- paste0(deres$name, "_heatmap", ".png")
  if(save) save_pheatmap(heatmap, save_name)

  return(heatmap)
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
  if(deres$species == "mouse") {
    positive_genes <- positive_genes %>% to_mouse_case(., "tomouse", "annodbi")
    negative_genes <- negative_genes %>% to_mouse_case(., "tomouse", "annodbi")
  }
  
  # get genes available
  pos_ids <- rownames(annotated_log2_cpm_df)[annotated_log2_cpm_df$gene_name %in% positive_genes]
  neg_ids <- rownames(annotated_log2_cpm_df)[annotated_log2_cpm_df$gene_name %in% negative_genes]
  
  # calc average
  avg_pos <- colMeans(log2_cpm[pos_ids, , drop = FALSE])
  avg_neg <- colMeans(log2_cpm[neg_ids, , drop = FALSE])
  
  # index = ratio of pos and neg
  matrix_index <- avg_pos / avg_neg
  
  return(matrix_index)
}

sample_index <- function(deres, plot=TRUE, save=FALSE, save_name=NULL) {
  # calculate index
  matrix_index <- matrisome_index(deres)
  
  metadata <- deres$dds %>% colData() %>% as.data.frame()
  metadata$index <- matrix_index[rownames(metadata)]
  
  # check for treatment status
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
  
  # check for treatment status
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
printer_analysis <- function(deres, output_dir = "printer_analyses", plot_all=FALSE) {
  # get title
  title <- deres$name
  
  # create dir to store all printer_analyses files
  base_dir <- paste0(output_dir, "/")
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  file_path <- paste0(base_dir, title)
  if (!dir.exists(file_path)) dir.create(file_path, recursive = TRUE)
  
  # get metadata from deres object
  metadata <- deres$dds@colData %>% as.data.frame()
  
  # matrisome_genes
  if(deres$species == "mouse") {
    matrisome_index <- to_mouse_case(matrisome_index, method = "tomouse", source = "annodbi")
    immuno_genes <- to_mouse_case(immuno_file$name, method = "tomouse", source = "annodbi")
    matrisome_genes <- matrisome_file$mouse
  } else {
    matrisome_genes <- matrisome_file$GeneSymbol
    immuno_genes <- immuno_file$name
  }
  
  ### Slide 1: Title
  title_slide <- function(title) {
    grid.newpage()
    grid.text(title, gp = gpar(fontsize = 24, fontface = "bold"))
  }
  
  ### Slide 2: Heading Slide
  add_heading_slide <- function(main_title, subtitle) {
    grid.newpage()
    grid.text(main_title, y = 0.65, gp = gpar(fontsize = 20, fontface = "bold"))
    grid.text(subtitle, y = 0.45, gp = gpar(fontsize = 14))
  }
  
  ### Slide 3: metadata
  df_png(metadata, paste0(file_path, "/metadata.png"), col_width = 200)
  
  ### Slide 4: PCA
  visualise_pca(deres, plot=plot_all, save=TRUE, save_name=paste0(file_path, "/pca.png"))
  
  ### Slide 5-6: Volcano Plots
  visualise_volcano(deres, plot=plot_all, save=TRUE, save_name=paste0(file_path, "/volcano.png"))
  visualise_volcano(deres, plot=plot_all, save=TRUE, save_name=paste0(file_path, "/matrisome_volcano.png"), proteins = matrisome_genes)
  visualise_volcano(deres, plot=plot_all, save=TRUE, save_name=paste0(file_path, "/immuno_volcano.png"), proteins = immuno_genes)
  
  ### Slide 6-7: Heatmaps
  visualise_heatmap(deres, save=TRUE, proteins=matrisome_index_genes, save_name=paste0(file_path, "/heatmap_matrisome_index.png"), plot=plot_all)
  visualise_heatmap(deres, save=TRUE, proteins=matrisome_genes, save_name=paste0(file_path, "/heatmap_matrisome.png"), plot=plot_all)
  visualise_heatmap(deres, save=TRUE, proteins=immuno_genes, save_name=paste0(file_path, "/heatmap_immuno.png"), plot=plot_all)
  
  ### Slide 8-11: Pathway Analysis
  # GSEA hallmark
  gseaRes <- pathway_enrichment(deres, view_all = FALSE)
  dotplot_pathways(gseaRes, plot = plot_all, save = TRUE, save_name=paste0(file_path, "/dotplot_hallmark.png"))
  
  # GSEA naba
  gseaResNaba <- pathway_enrichment(deres, view_all = TRUE, category = "C2", subset = naba_pathways)
  dotplot_pathways(gseaResNaba, plot = plot_all, save = TRUE, save_name=paste0(file_path, "/dotplot_naba.png"))
  
  # GSEA custom matrisome
  gseaResMatrisome <- matrisome_pathway_enrichment(deres)
  dotplot_pathways(gseaResMatrisome, plot = plot_all, save = TRUE, save_name=paste0(file_path, "/dotplot_matrisome.png"))
  
  ### Slide 12-13: Matrix Index
  sample_index(deres, plot = plot_all, save = TRUE, save_name = paste0(file_path, "/sample_index.png"))
  boxplot_index(deres, plot = plot_all, save = TRUE, save_name = paste0(file_path, "/boxplot_index.png"))
  
  ### Slide 14-17: Deconvolute Significant
  volcano_deconvolution(deres, plot=plot_all, method="xcell", save=TRUE, save_name=paste0(file_path, "/xcell_volcano.png"))
  volcano_deconvolution(deres, plot=plot_all, method="mcp", save=TRUE, save_name=paste0(file_path, "/mcp_volcano.png"))
  volcano_deconvolution(deres, plot=plot_all, method="epic", save=TRUE, save_name=paste0(file_path, "/epic_volcano.png"))
  volcano_deconvolution(deres, plot=plot_all, method="decosus", save=TRUE, save_name=paste0(file_path, "/decosus_volcano.png"))
  
  ### Slide 18-21: Cluster Deconvolution
  cluster_deconvolution_individual(deres, method="xcell", plot=plot_all, save=TRUE, save_name=paste0(file_path, "/xcell_heatmap.png"))
  cluster_deconvolution_individual(deres, method="mcp", plot=plot_all, save=TRUE, save_name=paste0(file_path, "/mcp_heatmap.png"))
  cluster_deconvolution_individual(deres, method="epic", plot=plot_all, save=TRUE, save_name=paste0(file_path, "/epic_heatmap.png"))
  cluster_deconvolution_individual(deres, method="decosus", plot=plot_all, save=TRUE, save_name=paste0(file_path, "/decosus_heatmap.png"))
  
  
  ### Slide 22-23: Significant Genes (only if more than 1)
  significant_genes <- view_res_table(deres) %>% dplyr::slice(1:15)
  if(nrow(significant_genes) > 1) {
    ### Slide 22: Significant Genes
    df_png(significant_genes, save_name=paste0(file_path, "/significant_genes.png"), col_width = 110)
    
    ### Slide 23: Heatmap significant
    visualise_heatmap(deres, save=TRUE, save_name=paste0(file_path, "/heatmap_significant.png"), plot=plot_all)
    
    plot_significant <- TRUE
  } else plot_significant <- FALSE
  
  # Load in data
  slide_files <- c(
    "metadata.png", "pca.png",
    "volcano.png", "matrisome_volcano.png", "immuno_volcano.png",
    "heatmap_matrisome_index.png", "heatmap_matrisome.png", "heatmap_immuno.png",
    "dotplot_hallmark.png", "dotplot_naba.png", "dotplot_matrisome.png",
    "sample_index.png", "boxplot_index.png",
    "xcell_volcano.png", "mcp_volcano.png", "epic_volcano.png", "decosus_volcano.png",
    "xcell_heatmap.png", "mcp_heatmap.png", "epic_heatmap.png", "decosus_heatmap.png"
  )
  if (plot_significant) slide_files <- c(slide_files, "significant_genes.png",
                                         "heatmap_significant.png")
  slide_files <- file.path(file_path, slide_files)
  
  # Create the PDF!
  pdf(file = paste0(file_path, ".pdf"), width = 8, height = 6)
  
  ## Slide 1: Title
  title_slide(title)
  
  
  # Slide 2+: Annotated image slides
  slide_titles <- c(
    "Metadata Table", "PCA Plot",
    "Volcano Plot", "Matrisome Volcano Plot", "Immuno Volcano Plot",
    "Heatmap: Matrisome Index Genes", "Heatmap: Matrisome Genes", "Heatmap: Immuno Genes",
    "Dotplot: Hallmark Pathways", "Dotplot: NABA Pathways", "Dotplot: Matrisome Pathways", 
    "Matrisome Index: Between Samples", "Matrisome Index: Between Response",
    "Deconvolute Significant: xCell", "Deconvolute Significant: MCPcounter", "Deconvolute Significant: EPIC", "Deconvolute Significant: Decosus",
    "Cluster: xCell", "Cluster: MCPcounter", "Cluster: EPIC", "Cluster: Decosus"
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

# analyse datasets
analyse_datasets <- function(datasets_list) {
  deres_list <- list()
  for (i in seq_along(datasets_list)) {
    dataset <- datasets_list[[i]]
    deres <- response_analysis(dataset$name, dataset$counts, dataset$metadata, dataset$species, "annodbi")
    deres_list[[dataset$name]] <- deres
  }
  
  return(deres_list)
}

# print datasets
print_deres_list <- function(deres_list, output_dir="printer_analyses",
                             output_file = "publicdata.pdf") {
  for (i in seq_along(deres_list)) {
    deres <- deres_list[[i]]
    for (i in seq_along(deres)) printer_analysis(deres[[i]], output_dir=output_dir)
  }
  
  combine_pdfs(pdfs_dir = output_dir, output_file = output_file)
}


# ---- 10. Combination Analysis  ----

# combination analysis printer function
combined_printer <- function(deres_list, pdf_name = "combined_analysis.pdf") {
  gsea_df <- prep_gsea_combined(deres_list)
  plot_gsea_combined(gsea_df, plot_type = 0, save = TRUE)
  plot_gsea_combined(gsea_df, plot_type = 1,save = TRUE)
  plot_gsea_combined(gsea_df, plot_type = 2,save = TRUE)
  
  genes_combined(deres_list, save=TRUE)
  genes_combined(deres_list, "all", pvalues_cutoff = 1, save=TRUE)
  other_genes_combined(deres_list, immuno_file$name, cutoff = 1, save=TRUE)
  index_combined(deres_list, save=TRUE)
  
  deconvolute_combined(deres_list, save=TRUE)
  deconvolute_combined(deres_list, method = "mcp", save=TRUE)
  deconvolute_combined(deres_list, method = "epic", save=TRUE)
  deconvolute_combined(deres_list, method = "decosus", save=TRUE)
  
  cluster_deconvolution_combined(deres_list, method = "xcell", save=TRUE)
  cluster_deconvolution_combined(deres_list, method = "mcp", save=TRUE)
  cluster_deconvolution_combined(deres_list, method = "epic", save=TRUE)
  cluster_deconvolution_combined(deres_list, method = "decosus", save=TRUE)
  
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
  
  if(genes == "index") combined_res_df <- combined_res_df[combined_res_df$GeneSymbol %in% matrisome_index_genes, ]
  
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
      strip.text = element_text(size = 10, angle = 0)
    )
  
  if(plot) plot(gene_plot)
  if(save) {
    if(genes == "index") ggsave("combined_analyses/index_genes_combined.png", gene_plot, scale = 3)
    else ggsave("combined_analyses/matrisome_genes_combined.png", gene_plot, scale=3)
  }
  
  return(significant_genes)
}

other_genes_combined <- function(deres_list, gene_list, cutoff=0, plot=TRUE, save=FALSE) {
  merged_res_list <- list()
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
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
  combined_res_df <- combined_res_df[combined_res_df$GeneSymbol %in% gene_list, ]
  
  significant_genes <- c()
  for(i in seq_along(matrisome_file$GeneSymbol)) {
    gene <- matrisome_file$GeneSymbol[i]
    filter_gene <- combined_res_df %>% filter(GeneSymbol == gene & treatment_status == "response")
    table_significance <- table(filter_gene$significance)
    if ( cutoff == 0 && 
         ("≥ 0.05" %in% names(table_significance) || "< 0.05" %in% names(table_significance))) {
      significant_genes <- c(significant_genes, gene) 
    }
    else if ("< 0.05" %in% names(table_significance) && table_significance["< 0.05"] >= cutoff) {
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
  
  gene_plot <- ggplot(filtered_res_df, aes(x = GeneSymbol, y = dataset_label)) +
    # dotplot
    geom_point(aes(size = significance, color = log2FoldChange)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5, "NA" = 1)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggh4x::facet_nested(rows = vars(treatment_status), scales = "free", space = "free") +
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_markdown(size = 8),
      strip.text = element_text(size = 10)
    )
  
  if(plot) plot(gene_plot)
  
  if(save) {
    ggsave("combined_analyses/other_genes_combined.png", gene_plot, scale = 3)
  }
  
  return(significant_genes)
}

index_combined <- function(deres_list, plot=TRUE, save=FALSE) {
  matrisome_index_list <- list()
  
  for (i in seq_along(deres_list)) {
    if("untreated" %in% names(deres_list[[i]])) entry <- deres_list[[i]]$untreated
    else if ("merged" %in% names(deres_list[[i]])) entry <- deres_list[[i]]$merged
    else entry <- deres_list[[i]]$response

    matrix_index <- matrisome_index(entry)
    
    metadata <- entry$dds %>% colData() %>% as.data.frame()
    metadata$index <- matrix_index[rownames(metadata)]
    metadata$response <- tolower(metadata$response)
    metadata$tissue <- entry$tissue
    metadata$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
    
    metadata <- metadata %>% select(index, dataset, response, tissue)
    
    matrisome_index_list[[entry$name]] <- metadata
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
prep_gsea_combined <- function(deres_list, species_prep=FALSE, tissue_prep=FALSE) {
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
        custom_enrich <- matrisome_pathway_enrichment(entry, view_all = TRUE)@result
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
  
  
  
  if(species_prep) {
    # colored labels for species
    combined_gsea_df$dataset_label <- case_when(
      combined_gsea_df$species == "human" ~ paste0("<span style='color:#2ca25f;'>", combined_gsea_df$dataset, "</span>"),
      combined_gsea_df$species == "mouse" ~ paste0("<span style='color:#de77ae;'>", combined_gsea_df$dataset, "</span>"),
      TRUE ~ combined_gsea_df$dataset
    )
  } else combined_gsea_df$dataset_label <- combined_gsea_df$dataset
  
  # order dataset by tissue
  if(tissue_prep) {
    combined_gsea_df$dataset_label <- factor(
      combined_gsea_df$dataset_label,
      levels = combined_gsea_df %>%
        distinct(dataset_label, tissue) %>%
        arrange(tissue, dataset_label) %>%
        pull(dataset_label)
    ) %>% paste0(combined_gsea_df$tissue, " | ", .)
  }
  
  # order ID by facet/category
  combined_gsea_df$ID <- factor(combined_gsea_df$ID, 
                                levels = unique(combined_gsea_df$ID))

  return(combined_gsea_df)
}

plot_gsea_combined <- function(combined_gsea_df, plot=TRUE, save=FALSE, plot_type=0) {
  # plot type: 0 = all, 1 = matrisome only, 2 = hallmark only
  if(plot_type == 1) combined_gsea_df <- combined_gsea_df %>% filter(facet != "Hallmark")
  else if(plot_type == 2) combined_gsea_df <- combined_gsea_df %>% filter(facet == "Hallmark")
  
  gsea_combined <- ggplot(combined_gsea_df, aes(x = ID, y = dataset_label)) +
    # dotplot
    geom_point(aes(size = significance, color = NES)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # nested facet plot
    ggh4x::facet_nested(rows = vars(treatment_status), cols = vars(facet, category), scales = "free", space = "free") + theme_pubr() + 
    
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=8),
      axis.text.y = element_markdown(size = 8),
      strip.text.x = element_text(size = 8),
      strip.text.y = element_text(angle=0)
    ) + labs(x = NULL, y = NULL) 

  
  if(plot) plot(gsea_combined)
  if(save) ggsave("combined_analyses/gsea_combined.png", gsea_combined, scale = 3)
}

plot_custom_gsea_combined <- function(deres_list, custom_pathway) {
  merged_gsea_list <- list()
  
  for (i in seq_along(deres_list)) {
    for (j in seq_along(deres_list[[i]])) {
      entry <- deres_list[[i]][[j]]
      treatment_status <- sub(".*_(.*)$", "\\1", entry$name)
      if (treatment_status == "response") {
        gsea <- custom_pathway_enrichment(entry, custom_pathway)@result
        if (nrow(gsea) > 0) {
          gsea$category <- "ML"
          
          gsea$treatment_status <- sub(".*_(.*)$", "\\1", entry$name)
          gsea$dataset <- sub("_(?=[^_]*$).*", "", entry$name, perl = TRUE)
          gsea$species <- entry$species
          gsea$tissue <- entry$tissue
          
          merged_gsea_list[[entry$name]] <- gsea
        }
      }
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
  
  gsea_combined <- ggplot(combined_gsea_df, aes(x = ID, y = dataset_label)) +
    # dotplot
    geom_point(aes(size = significance, color = NES)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_markdown(size = 8),
      strip.text = element_text(size = 10)
    ) + labs(x = NULL, y = NULL)
  
  plot(gsea_combined)
}

# ---- 10.3 Common Genes  ----
# find commonly upregulated and downregulated genes
common_genes <- function(deres_list, direction="up", cutoff=2) {
  genes <- c()
  for (i in seq_along(deres_list)) {
    if("untreated" %in% names(deres_list[[i]])) entry <- deres_list[[i]]$untreated
    else if("merged" %in% names(deres_list[[i]])) entry <- NULL
    else entry <- deres_list[[i]]$response
    
    if(!is.null(entry)) {
      res_genes <- view_res_table(entry, direction)$gene_name
      if(entry$species == "mouse") res_genes <- to_mouse_case(res_genes, "tohuman", "annodbi")
      res_genes <- res_genes[!is.na(res_genes)]
      genes <- c(genes, res_genes)
    }
  }
  gene_counts <- table(genes)
  common <- gene_counts[gene_counts >= cutoff]
  return(sort(common, decreasing = TRUE))
}

# ---- 10.4 Deconvolution Combined  ----
deconvolute_combined <- function(deres_list, method = "xcell", plot=TRUE, save=FALSE, tissue_prep=FALSE) {
  merged_deconvolution_list <- list()
  
  for (i in seq_along(deres_list)) {
    if("untreated" %in% names(deres_list[[i]])) entry <- deres_list[[i]]$untreated
    if("merged" %in% names(deres_list[[i]])) entry <- deres_list[[i]]$merged
    else entry <- deres_list[[i]]$response
      
    if (entry$species == "human") {
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
  
  # Combine all data and continue as normal
  combined_deconv_df <- do.call(rbind, merged_deconvolution_list) %>%
    mutate(log_p = -log10(p_adj))
  
  combined_deconv_df$significance <- apply(combined_deconv_df, 1, function(row) {
    if(is.na(row[["p_adj"]])) return("≥ 0.05")
    else if(row[["p_adj"]] < 0.05) return("< 0.05")
    else return("≥ 0.05")
  })
  
  # Order dataset by tissue
  if(tissue_prep){
    combined_deconv_df$dataset <- factor(
      combined_deconv_df$dataset,
      levels = combined_deconv_df %>%
        distinct(dataset, tissue) %>%
        arrange(tissue, dataset) %>%
        pull(dataset)
    ) %>% paste0(combined_deconv_df$tissue, " | ", .)
  }
  
  
  # Plot
  combined_deconv_plot <- ggplot(combined_deconv_df, aes(x = cell, y = dataset)) +
    geom_point(aes(size = significance, color=log2_fc)) + 
    scale_size_manual(values = c("≥ 0.05" = 2, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if(plot) plot(combined_deconv_plot)
  if(save) ggsave( paste0("combined_analyses/", method, "_deconv_combined.png"), combined_deconv_plot, scale = 3)
}

cluster_deconvolution_combined <- function(deres_list, method="xcell", clusters = 3, plot=TRUE, save=FALSE) {
  set.seed(42)
  
  # get merged
  merged_dataset <- merge_datasets(deres_list)
  merged_dataset$counts <- merged_dataset$counts %>% raw_to_tpm() 
  combat_corr <- combat_correction(merged_dataset, FALSE)
  
  expr <- combat_corr$expr
  metadata <- combat_corr$metadata
  
  deconvolution <- deconvolute(expr, method)
  
  # heatmap
  ann_colors <- list(
    response = c("responder" = "forestgreen", "nonresponder" = "firebrick")
  )
  
  top_anno <- HeatmapAnnotation(
    response = metadata$response,
    batch = metadata$batch,
    col = ann_colors,
    annotation_legend_param = list(
      response = list(title = "Response")
    )
  )
  
  
  heatmap <- Heatmap(
    deconvolution %>% as.matrix(),
    top_annotation = top_anno,
    cluster_columns = TRUE,
    column_km = clusters,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = TRUE,
    heatmap_legend_param = list(color_bar = "continuous")
  )
  
  if(plot) plot(heatmap)

  if(save){
    png(paste0("combined_analyses/", method, "_deconv_cluster.png"),
        width=10,height=10,units="in",res=800)
    print(heatmap)
    dev.off()
  }  
}

# ---- 11. Merger Analysis ----

merge_datasets <- function(deres_list, tissue="pancreas") {
  required_deres <- list()
  # get required datasets
  for (i in seq_along(deres_list)) {
    deres <- deres_list[[i]]$response
  
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
    deres_metadata <- deres$dds %>% colData() %>% as.data.frame()
    if ("treatment" %in% colnames(deres_metadata)) required_met <- deres_metadata %>% 
      filter(tolower(treatment) == "untreated") %>% select(response, batch)
    else required_met <- deres_metadata  %>% select(response, batch)
    
    metadata <- rbind(metadata, required_met)
    
    required_counts <- deres$dds %>% assay() %>% reverse_human_annotate()
    required_counts$gene_id <- rownames(required_counts)
    required_counts <- required_counts %>% select(c(rownames(required_met), gene_id))
    
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

# ---- 12 Machine Learning  ----

## Functions to prepare datasets
# batch correct vsd
combat_correction <- function(merged_dataset, normalise=TRUE) {
  # annotation
  raw_counts <- merged_dataset$counts %>% annotate_df()
  raw_counts <- raw_counts[!is.na(raw_counts$gene_name), ]
  rownames(raw_counts) <- make.unique(as.character(raw_counts$gene_name))
  raw_counts$gene_name <- NULL
  
  # normalisation
  if(normalise) {
    # compute cpm
    lib_sizes <- colSums(raw_counts)
    cpm <- sweep(raw_counts, 2, lib_sizes, FUN = "/") * 1e6
    # log transform
    transformed_cpm <- log2(cpm + 1)
    transformed_cpm <- as.matrix(transformed_cpm)
  } else transformed_cpm <- as.matrix(raw_counts)
  
  
  # Batch Correction
  BiocParallel::register(SerialParam())
  mod <- model.matrix(~ response, data = merged_dataset$metadata)
  combat_expr <- ComBat(dat = transformed_cpm, batch = merged_dataset$metadata$batch, mod = mod)
  
  return(list(expr=combat_expr, metadata=merged_dataset$metadata))
}

# prepare unseen data counts
prepare_new_data <- function(counts) {
  lib_sizes <- colSums(counts)
  cpm <- sweep(counts, 2, lib_sizes, FUN = "/") * 1e6
  log2_cpm <- log2(cpm + 1) %>% as.data.frame() %>% annotate_df_remove_dup()
  
  return(log2_cpm)
}

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

prepare_ml_data <- function(merged_dataset, use_matrisome=FALSE, top_n = 1000, use_mads=TRUE) {
  cor_data <- combat_correction(merged_dataset)
  
  if(isTRUE(use_matrisome)) {cor_data$expr <- cor_data$expr %>% merged_filter_matrisome()}
  else {
    # Select top N most variable genes using MAD or variance
    if (use_mads) {
      gene_mads <- apply(cor_data$expr, 1, mad)
      top_genes <- names(sort(gene_mads, decreasing = TRUE))[1:top_n]
    } else {
      gene_vars <- apply(cor_data$expr, 1, var)
      top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:top_n]
    }
    cor_data$expr <- cor_data$expr[top_genes, ]
    
    cor_data$expr <- cor_data$expr %>% as.data.frame() %>% annotate_df("annodbi", "human")
    cor_data$expr <- cor_data$expr[!duplicated(cor_data$expr$gene_name), ]
    
    rownames(cor_data$expr) <- cor_data$expr$gene_name
    cor_data$expr$gene_name <- NULL
  }
  
  df <- as.data.frame(t(cor_data$expr))  # samples as rows
  df$response <- as.factor(cor_data$metadata$response)
  return(df)
}

pca_merged <- function(merged_dataset) {
  # Perform PCA
  ml_df <- prepare_ml_data(merged_dataset)
  pca_res <- PCA(ml_df[ , !colnames(ml_df) %in% "response"], graph = FALSE)
  
  # Visualise
  fviz_pca_ind(pca_res,
               geom.ind = "point",
               habillage = ml_df$response,
               addEllipses = TRUE,
               palette = "jco")
}

# ---- 12.1 ML Models ----
run_ml_model <- function(df, method = "rf", tuneLength = 5, seed = 42) {
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
  
  df$response <- factor(df$response, levels = c("nonresponder", "responder"))
  
  trainIndex <- createDataPartition(df$response, p = 0.7, list = FALSE)
  train <- df[trainIndex, ]
  test <- df[-trainIndex, ]
  
  trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                         classProbs = TRUE, summaryFunction = twoClassSummary,
                         savePredictions = "final")
  
  results <- list()
  
  for (model in models) {
    message("Training model: ", model)
    
    if(model == "nnet") fit <- train(response ~ ., data = train, method = model,
                 trControl = trctrl, metric = "ROC", tuneLength = tuneLength)
    else fit <- train(response ~ ., data = train, method = model,
                      trControl = trctrl, metric = "ROC", tuneLength = tuneLength, preProcess = c("center", "scale"))
    
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
      head(top_n) %>%
      filter(Overall > 0)
    
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

view_compare_ml <- function(deres_list, use_matrisome=TRUE, tissue="pancreas",
                            no_features=60) {
  merged_dataset <- merge_datasets(deres_list, tissue=tissue)
  ml_df <- prepare_ml_data(merged_dataset, use_matrisome = use_matrisome)
  
  # model comparison
  model_results <- compare_ml_models(ml_df)
  
  # visualise
  plot_model_aucs(model_results)
  
  # get top features
  top_features_list <- extract_model_features(model_results, top_n = no_features, plot=TRUE)
  
  merge_features_df <- lapply(top_features_list, function(df) {
    df <- df[, c("Model", "Feature")]  # Keep only necessary columns (optional)
    names(df)[names(df) == "Model"] <- "term"
    names(df)[names(df) == "Feature"] <- "gene"
    return(df)
  })
  merged_features <- do.call(rbind, merge_features_df)
  
  plot_custom_gsea_combined(deres_list, merged_features)
  return(merged_features)
}

# ---- 12.2 ML Gene Signatures ----
# create gene signature using caret
caret_gene_signature <- function(df, seed=42) {
  set.seed(seed)
  
  # 1. Prepare training data
  df$response <- factor(df$response, levels = c("nonresponder", "responder"))
  
  # 2. Define training control
  # model evaluation CV
  train_control <- caret::trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  # includes CV for hyperparameter (alpha, lambda) tuning
  model_glmnet <- train(
    response ~ .,
    data = df,
    method = "glmnet",
    family = "binomial",
    trControl = train_control,
    metric = "ROC",
    tuneLength = 20, 
    preProcess = c("center", "scale")
  )
  
  # 3. Extract Gene Signature
  
  # Best lambda value
  best_lambda <- model_glmnet$bestTune$lambda
  
  # Coefficients at best lambda
  coefs <- coef(model_glmnet$finalModel, s = best_lambda)
  
  # Convert to dataframe and filter
  coef_df <- as.data.frame(as.matrix(coefs))
  coef_df$Gene <- rownames(coef_df)
  colnames(coef_df)[1] <- "Weight"
  
  # Keep non-zero (excluding intercept)
  gene_signature <- coef_df %>%
    filter(Weight != 0 & Gene != "(Intercept)")
  
  return(gene_signature)
}

signature_finder <- function(deres_list, tissue="all", use_matrisome = TRUE) {
  merged_dataset <- merge_datasets(deres_list, tissue = tissue)
  ml_df <- prepare_ml_data(merged_dataset, use_matrisome)
  signature <- caret_gene_signature(ml_df)
  
  return(signature)
}

# apply various calculated signatures
signature_applier <- function(signature, raw_counts, metadata=NULL) {
  prepped_counts <- prepare_new_data(raw_counts)
  positive_genes <- signature %>% filter(Weight > 0)
  negative_genes <- signature %>% filter(Weight < 0)
  
  # coef method
  coef_weights <- signature %>% filter(Gene %in% rownames(prepped_counts))
  coef_scores <- colSums(prepped_counts[coef_weights$Gene, , drop = FALSE] * coef_weights$Weight)
  
  # mean diff
  pos_means <- colMeans(prepped_counts[positive_genes$Gene, , drop = FALSE], na.rm = TRUE)
  neg_means <- colMeans(prepped_counts[negative_genes$Gene, , drop = FALSE], na.rm = TRUE)
  mean_diff <- pos_means - neg_means
  
  # mean ratio
  mean_ratio <- pos_means / neg_means
  
  # log-ratio
  log_ratio <- log2(pos_means + 1) - log2(neg_means + 1)
  
  # ssgsea
  gene_sets <- list(sig_pos=positive_genes$Gene, sig_neg=negative_genes$Gene)
  ssgsea <- ssgseaParam(prepped_counts %>% as.matrix(), gene_sets, 
                        minSize = 1, maxSize = Inf,
                        alpha = 0.25, normalize = TRUE)
  ssgsea_scores <- gsva(ssgsea, verbose = TRUE) %>% t() %>% as.data.frame() %>%
    rownames_to_column("sample") %>% pivot_longer(-sample, names_to = "direction", values_to = "ssgsea_score")
  ssgsea_pos <- ssgsea_scores %>% filter(direction == "sig_pos")
  ssgsea_neg <- ssgsea_scores %>% filter(direction == "sig_neg")
  ssgsea_sig_score <- ssgsea_pos$ssgsea_score - ssgsea_neg$ssgsea_score
  
  score_df <- data.frame(
    sample = colnames(prepped_counts),
    coef = coef_scores,
    mean_diff = mean_diff,
    mean_ratio = mean_ratio,
    log_ratio = log_ratio,
    ssgsea = ssgsea_sig_score
  )
  if(!is.null(metadata)) score_df$response <- metadata$response
  
  return(score_df)
}

# view clusterprofiler plots of signature datasets
signature_dotplot <- function(signature, deres_list) {
  # include direction
  direction_features <- signature
  direction_features$term <- ifelse(direction_features$Weight > 0, "Signature_Positive", "Signature_Negative")
  direction_features$gene <- direction_features$Gene
  direction_features <- direction_features %>% select(term, gene)
  
  # include all
  all_features <- signature
  all_features$term <- "Signature_All"
  all_features$gene <- all_features$Gene
  all_features <- all_features %>% select(term, gene)
  
  # combine
  required_features <- rbind(direction_features, all_features)
  
  # dotplot
  plot_custom_gsea_combined(deres_list, required_features)
}

# view all types of signatures calculated
view_signatures <- function(signature, counts, metadata = NULL) {
  signs_plot <- signature_applier(signature, counts, metadata) %>%
    pivot_longer(cols = -c(sample, response), names_to = "method", values_to = "score") %>% 
    ggplot(aes(x = score, fill = response)) + 
    facet_wrap(vars(method)) + 
    geom_density(alpha = 0.5)
  
  return(signs_plot)
}

# ---- 12.3 Immune-Response-related Risk Score  ----
get_intersection_genes <- function(deres_list, selected_deres, intersection_list) {
  deres <- deres_list[[selected_deres]]
  if("untreated" %in% names(deres)) res <- deres$untreated$res
  else res <- deres$response$res
  
  matrisome_genes <- read.csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv") %>% .$GeneSymbol
  
  intersection_list <- res %>% filter(abs(log2FoldChange) > 1 & padj < 0.05) %>% filter(gene_name %in% intersection_list) %>% .$gene_name
  
  return(intersection_list)
}


# ---- 12.4 Depreciated ML methods ----
compare_sig_pipeline <- function(deres_list, use_matrisome = TRUE, 
                               test_counts, test_metadata) {
  # 1. prepare data
  merged_dataset <- merge_datasets(deres_list)
  ml_df <- prepare_ml_data(merged_dataset, use_matrisome)
  
  # 2. build signature
  if (method=="caret") sig_res <- caret_gene_signature(ml_df)
  else sig_res <- glmnet_gene_signature(ml_df)
  
  # 3. apply to train data
  train_scores <- signature_applier(sig_res, merged_dataset$counts, merged_dataset$metadata) %>% select(sample, coef, response)
  
  train_scores$Z <- scale(train_scores$coef)
  
  mean_train <- mean(train_scores$coef)
  sd_train   <- sd(train_scores$coef)
  
  train_roc_obj <- roc(train_scores$response, train_scores$coef)
  train_auc_val <- auc(train_roc_obj)
  
  # if (method=="batch") test_scores <- prepare_new_data_batch(test_counts, test_metadata, signature_data$train_counts, signature_data$train_metadata)

  test_scores <- signature_applier(signature_data$signature, test_counts, test_metadata) %>% select(sample, coef, response)
  
  test_scores$Z <- (test_scores$coef - mean_train) / sd_train
  
  train_plot <- ggplot(train_scores, aes(x = Z, fill = response)) +
    geom_vline(xintercept = mean(train_scores$Z), linetype = "dashed") +
    geom_density(alpha = 0.5) +
    ggtitle("Train data")
  
  test_plot <- ggplot(test_scores, aes(x = Z, fill = response)) +
    geom_vline(xintercept = mean(train_scores$Z), linetype = "dashed") +
    geom_density(alpha = 0.5) +
    ggtitle("Test Data")
  
  test_roc_obj <- roc(test_scores$response, test_scores$coef)
  test_auc_val <- auc(test_roc_obj)
  
  plot(ggarrange(train_plot, test_plot))

  return(train_auc_val, test_auc_val)
}

signature_lobo_validate <- function(deres_list, tissue="all", use_matrisome=TRUE) {
  set.seed(42)
  
  merged_dataset <- merge_datasets(deres_list, tissue = tissue)
  ml_df <- prepare_ml_data(merged_dataset, use_matrisome = use_matrisome)
  
  metadata <- merged_dataset$metadata
  batches <- unique(metadata$batch)
  
  results <- list()
  
  for (b in batches) {
    message("Leaving out batch: ", b)
    
    # Identify samples
    test_ids <- rownames(metadata[metadata$batch == b, ])
    train_ids <- rownames(metadata[metadata$batch != b, ])
    
    train_df <- ml_df[train_ids, ]
    test_df  <- ml_df[test_ids, ]
    
    # Train signature
    sig <- caret_gene_signature(train_df)
    
    # Score training
    train_scores <- apply_gene_signature(train_df[, -which(names(train_df) == "response")], sig)
    train_scores$Z <- scale(train_scores$Score)
    train_scores$response <- train_df$response
    
    # Score testing (Z-normalize using train stats)
    mean_train <- mean(train_scores$Score)
    sd_train   <- sd(train_scores$Score)
    
    test_scores <- apply_gene_signature(test_df[, -which(names(test_df) == "response")], sig)
    test_scores$Z <- (test_scores$Score - mean_train) / sd_train
    test_scores$response <- test_df$response
    
    # ROC
    roc_obj <- roc(response = test_scores$response, predictor = test_scores$Score)
    auc_val <- auc(roc_obj)
    
    # Save
    results[[b]] <- list(
      batch = b,
      train_scores = train_scores,
      test_scores = test_scores,
      auc = auc_val,
      roc = roc_obj,
      signature = sig
    )
  }
  
  return(results)
}

plot_lobo_aucs <- function(lobo_results) {
  aucs <- sapply(lobo_results, function(x) x$auc)
  
  # Barplot AUCs
  barplot(aucs, beside = TRUE, col = "steelblue",
          main = "LOBO AUCs per Batch", ylab = "AUC", ylim = c(0, 1))
}

plot_lobo_results <- function(lobo_results) {
  plots <- list()
  
  for (batch_name in names(lobo_results)) {
    res <- lobo_results[[batch_name]]
    train_df <- res$train_scores
    test_df <- res$test_scores
    
    # Train plot
    p_train <- ggplot(train_df, aes(x = Z, fill = response)) +
      geom_density(alpha = 0.5) +
      ggtitle(paste("Train (excluding", batch_name, ")")) +
      theme_minimal()
    
    # Test plot
    p_test <- ggplot(test_df, aes(x = Z, fill = response)) +
      geom_density(alpha = 0.5) +
      ggtitle(paste("Test (batch =", batch_name, ")\nAUC =", round(res$auc, 3))) +
      theme_minimal()
    
    # Arrange side by side
    plots[[batch_name]] <- ggarrange(p_train, p_test, ncol = 2)
  }
  
  # Combine all arranged plots into a grid
  final_plot <- ggarrange(plotlist = plots, ncol = 1)
  return(final_plot)
}

# ---- 13. Other Useful Functions  ----
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

# ----- Load all packages  ----

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)
library(ggpubr)

library(NMF)
library(limma)
library(glue)
library(fgsea)
library(nnls)

# ----- 1. Download TCGA PAAD data  ----

download_data <- function() {
  paad_query <- GDCquery(
    project = "TCGA-PAAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts",
    access = "open",
  )
  GDCdownload(query = paad_query, method = 'client')
  data <- GDCprepare(query = paad_query)
  saveRDS(data, "paad_data.rds")
  
  return(data)
}

# load data (if exists otherwise download)
load_data <- function() {
  if (file.exists("paad_data.rds")) {
    data <- readRDS("paad_data.rds")
  } else {
    print("downloading data")
    data <- download_data()
  }
  
  pdac_labels <- c(
    "infiltrating duct carcinoma, nos",
    "adenocarcinoma, nos",
    "adenocarcinoma with mixed subtypes"
  )
  
  data <- data[, (colData(data)$tissue_type %>% tolower()) != "normal"]
  data <- data[, (colData(data)$classification_of_tumor %>% tolower()) != "metastasis"]
  data <- data[, (colData(data)$primary_diagnosis %>% tolower()) ==  "infiltrating duct carcinoma, nos"]
  
  
  # Filter the data
  
  return(data)
}

# annotate df
get_annotated_data <- function() {
  data <- load_data()
  
  # extract expression + survival
  expr <- assay(data)
  clinical <- colData(data)
  
  # get gene ensembl_id
  annotated <- expr %>% as.data.frame() %>% annotate_df()
  
  return(list(expr=expr, annotated=annotated, clinical=clinical))
}

# ----- 2. Survival Analyses  ----

# Gene data survival
plot_survival <- function(gene_selected, annotated_data = get_annotated_data()) {
  clinical <- annotated_data$clinical
  expr <- annotated_data$expr
  annotated <- annotated_data$annotated
  
  gene_row <- annotated %>%
    as.data.frame() %>%
    filter(gene_name == gene_selected) %>%
    slice(1)
  ensembl_id <- rownames(gene_row)
  
  # extract expression 
  expr_vector <- expr[ensembl_id, ]
  expr_group <- ifelse(expr_vector > median(expr_vector, na.rm = TRUE), "High", "Low")
  
  # get clincial data
  surv_df <- clinical %>% as.data.frame() %>%
    mutate(group = expr_group[colnames(expr)])
  
  # get plot
  fit <- survfit(Surv(days_to_death, vital_status == "Dead") ~ group, data = surv_df)
  plot <- ggsurvplot(fit, data = surv_df, pval = TRUE,
             title = paste0("Survival by ", gene_selected,  " Expression"))$plot
  
  return(plot)
}

# multiple genes survival
plot_survivals <- function(gene_list) {
  annotated_data <- get_annotated_data()
  
  km_plots <- list()
  for (gene in gene_list) {
    km_plots[[gene]] <- plot_survival(gene, annotated_data)
  }
  
  ggarrange(plotlist = km_plots)
}

find_sig_survivals <- function(gene_list, annotated_data = get_annotated_data()) {
  clinical <- annotated_data$clinical
  expr <- annotated_data$expr
  annotated <- annotated_data$annotated
  
  sig_genes <- c()
  
  for (gene in gene_list) {
    gene_row <- annotated %>%
      as.data.frame() %>%
      filter(gene_name == gene) %>%
      slice(1)
    ensembl_id <- rownames(gene_row)
  
    # extract expression 
    expr_vector <- expr[ensembl_id, ]
    expr_group <- ifelse(expr_vector > median(expr_vector, na.rm = TRUE), "High", "Low")
    
    # get clincial data
    surv_df <- clinical %>% as.data.frame() %>%
      mutate(group = expr_group[colnames(expr)])
    
    # get survival_data
    fit <- tryCatch(
      {survfit(Surv(days_to_death, vital_status == "Dead") ~ group, data = surv_df)
      }, error = function(error) {
        print(paste0("Failed ", gene))
        return(NULL)}
    )
    
    if(!is.null(fit)) {
      p_val_data <- surv_pvalue(fit, surv_df)
      p_val <- p_val_data$pval[[1]]
      if(!is.na(p_val) & p_val < 0.05) sig_genes <- c(gene, sig_genes) 
    }
  }
  
  return(sig_genes)
}

plot_gene_set_survivals <- function(gene_set, annotated_data = get_annotated_data()) {
  clinical <- annotated_data$clinical
  expr <- annotated_data$expr %>% prepare_new_data() %>% as.matrix()
  annotated <- annotated_data$annotated
  
  ssgsea <- ssgseaParam(
    expr, list(set = gene_set),
    minSize = 1, maxSize = Inf,
    alpha = 0.25, normalize = TRUE
  )
  
  # Run GSVA with the ssGSEA parameter object
  ssgsea_scores <- gsva(ssgsea, verbose = TRUE)
  
  # Extract the scores
  score_vector <- as.numeric(ssgsea_scores["set", ])
  names(score_vector) <- colnames(expr)
  
  # Group samples
  expr_group <- ifelse(score_vector > median(score_vector, na.rm = TRUE), "High", "Low")
  
  # Match group with clinical data
  surv_df <- clinical %>%
    as.data.frame() %>%
    mutate(group = expr_group[rownames(clinical)])
  
  # Kaplan-Meier survival fit
  fit <- survfit(Surv(days_to_death, vital_status == "Dead") ~ group, data = surv_df)
  # Plot survival
  plot <- ggsurvplot(fit, data = surv_df, pval = TRUE,
                     title = paste0("Survival by ssGSEA Score of Gene Set"))$plot
  
  return(plot)
}

# signature score survival
plot_survival_signature <- function(signature, method="coef", title=FALSE) {
  expr_data <- load_data() %>% assay() %>% as.data.frame()
  
  # get score
  score <- signature_applier(signature, expr_data)
  
  clinical <- load_data() %>% colData()
  
  # set expression group
  expr_group <- ifelse(score[[method]] > median(score[[method]], na.rm=TRUE), "High", "Low")
  
  # get survival data
  surv_df <- clinical %>% as.data.frame()
  surv_df$group <- expr_group
  fit <- survfit(Surv(days_to_death, vital_status == "Dead") ~ group, data = surv_df)
  
  if(title) plot_title <- paste0("Survival by Signature Score")
  else plot_title <- NULL
  
  # get plot
  plot <- ggsurvplot(fit, data = surv_df, pval = TRUE,
                     title = plot_title)$plot
  
  return(plot)
}

# ----- 3. NMF Analysis  ----


# ----- 3.1 NMF prep  ----
prepare_d_exp <- function() {
  exp_data <- load_data() %>% assay() %>% as.data.frame() %>% 
    prepare_new_data()
  
  d_exp <- exp_data
  mads <- apply(d_exp, 1, mad)
  d_exp <- d_exp[rev(order(mads))[1:1500],]
  
  return(d_exp)
}

matrisome_prepare_d_exp <- function() {
  exp_data <- load_data() %>% assay() %>% prepare_new_data()
  
  matrisome_file <- read_csv("~/Repos/RNASeq Pipeline/Matrisome_Hs_MasterList_SHORT.csv")
  
  d_exp <- exp_data %>% filter(rownames(.) %in% matrisome_file$GeneSymbol)
  
  zero_rows <- apply(d_exp, 1, function(row) all(row == 0))
  d_exp <- d_exp[!zero_rows, ]
  
  return(d_exp)
}
 

# ----- 3.2 Core NMF  ----
nmf_check <- function(use_matrisome=FALSE) {
  if(!use_matrisome) d_exp <- prepare_d_exp()
  else d_exp <- matrisome_prepare_d_exp()
  
  estim.r <- nmf(d_exp, 2:6, nrun=10, seed=42)
  plot(estim.r)
  consensusmap(estim.r)
  
  return(estim.r)
}

nmf_run <- function(num_clusters, use_matrisome=FALSE) {
  if(!use_matrisome) d_exp <- prepare_d_exp()
  else d_exp <- matrisome_prepare_d_exp()
  
  res <- nmf(d_exp, num_clusters, "brunet", seed=42)
  w <- basis(res)
  h <- coef(res)
  
  s <- extractFeatures(res)
  gene_lists <- list()
  for (i in 1:num_clusters) {
    gene_lists[[i]] <- rownames(d_exp)[s[[i]]]
  }
  sample_membership <- predict(res)
  
  limma_res <- limma_cluster_analysis(d_exp, sample_membership)
  gsea_res <- nmf_pathway_analysis(limma_res)
  
  return(list(original=d_exp, res=res, gene_lists=gene_lists, limma_res=limma_res, gsea_res=gsea_res))
}

nmf_survival <- function(nmf_data) {
  clinical <- load_data() %>% colData %>% as.data.frame()
  
  sample_membership <- predict(nmf_data$res)
  sample_membership_df <- data.frame(
    sample_id = names(sample_membership),
    cluster = factor(sample_membership)
  )
  clinical$sample_id <- clinical$barcode
  
  
  merged_df <- merge(clinical, sample_membership_df, by = "sample_id")
  
  merged_df$OS_time <- as.numeric(merged_df$days_to_death)
  merged_df$OS_event <- ifelse(is.na(merged_df$OS_time), 0, 1)
  merged_df$OS_time[is.na(merged_df$OS_time)] <- as.numeric(merged_df$days_to_last_follow_up[is.na(merged_df$OS_time)])
  
  fit <- survfit(Surv(OS_time, OS_event) ~ cluster, data = merged_df)
  
  ggsurvplot(fit,
             data = merged_df,
             pval = TRUE,
             ggtheme = theme_minimal())
}

get_cluster_defining_genes <- function(nmf_data, matrisome=TRUE, no_genes=40) {
  cluster_defining_genes <- list()
  for (i in 1:length(nmf_data$limma_res)) {
    limma_res <- nmf_data$limma_res[[i]]
    
    if(matrisome) {
      pos_genes <- limma_res %>% filter(rownames(.) %in% matrisome_file$GeneSymbol) %>% filter(logFC > 1) %>% arrange(desc(logFC)) %>% slice(1:no_genes) %>% rownames()
      neg_genes <- limma_res%>% filter(rownames(.) %in% matrisome_file$GeneSymbol) %>% filter(logFC < -1) %>% arrange(logFC) %>% slice(1:no_genes) %>% rownames()
    } else {
      pos_genes <- limma_res %>% filter(logFC > 1) %>% arrange(desc(logFC)) %>% slice(1:no_genes) %>% rownames()
      neg_genes <- limma_res %>% filter(logFC < -1) %>% arrange(logFC) %>% slice(1:no_genes) %>% rownames()
    }
    
    
    cluster_defining_genes[[paste0("Positive_Cluster_",i)]] <- pos_genes
    cluster_defining_genes[[paste0("Negative_Cluster_",i)]] <- neg_genes
  }
  
  return(cluster_defining_genes)
}

predict_clusters <- function(nmf_data, deres_list, matrisome=TRUE, no_genes=40) {
  merged_dataset <- merge_datasets(deres_list)
  prepped_counts <- prepare_new_data(merged_dataset$counts)
  n_clusters <- length(nmf_data$limma_res)
  
  cluster_defining_genes <- get_cluster_defining_genes(nmf_data, matrisome, no_genes)
  
  ssgsea <- ssgseaParam(prepped_counts %>% as.matrix(), cluster_defining_genes, 
                        minSize = 1, maxSize = Inf,
                        alpha = 0.25, normalize = TRUE)
  ssgsea_scores <- gsva(ssgsea, verbose = TRUE) %>% t() %>% as.data.frame()
  
  cluster_scores <- sapply(1:n_clusters, function(i) {
    ssgsea_scores[[paste0("Positive_Cluster_", i)]] - ssgsea_scores[[paste0("Negative_Cluster_", i) ]]
  })

  cluster_scores_df <- cluster_scores %>% as.data.frame()
  colnames(cluster_scores_df) <- paste0("Cluster_", 1:length(nmf_data$limma_res))
  rownames(cluster_scores_df) <- rownames(ssgsea_scores)
  
  predicted_clusters <- apply(cluster_scores_df, 1, function(x) {
    names(x)[which.max(x)]
  })
  
  return(data.frame(sample_id=rownames(cluster_scores_df), cluster=factor(predicted_clusters)))
}

predict_clusters_nnls <- function(nmf_data, deres_list) {
  # Step 1: Prepare external expression data 
  merged_dataset <- merge_datasets(deres_list)
  
  merged_counts <- merged_dataset$counts %>%
    prepare_new_data()
  
  # Step 2: Match genes with original NMF data
  d_exp <- nmf_data$original
  
  common_genes <- intersect(rownames(d_exp), rownames(merged_counts))
  merged_counts_matched <- merged_counts[common_genes, ]
  d_exp_matched <- d_exp[common_genes, ]
  
  # Step 3: Match genes with NMF basis (W) 
  res <- nmf_data$res
  W <- basis(res)
  
  common_genes <- intersect(rownames(W), rownames(merged_counts_matched))
  W_matched <- as.matrix(W[common_genes, ])
  new_expr <- as.matrix(merged_counts_matched[common_genes, ])  # genes x samples
  
  # Step 4: Assign clusters using NNLS
  assign_clusters_nnls <- function(expr_mat, W_basis) {
    expr_mat <- as.matrix(expr_mat)
    W_basis <- as.matrix(W_basis)
    
    cluster_scores <- apply(expr_mat, 2, function(sample_vec) {
      fit <- nnls(W_basis, as.numeric(sample_vec))
      return(fit$x)
    })
    
    cluster_scores <- t(cluster_scores)  # samples x clusters
    assignments <- max.col(cluster_scores)
    
    return(list(scores = cluster_scores, assignments = assignments))
  }
  
  nnls_result <- assign_clusters_nnls(new_expr, W_matched)
  assigned_clusters <- nnls_result$assignments
  names(assigned_clusters) <- colnames(new_expr)
  
  # Step 5: Merge with external response metadata
  # (Assumes response_meta exists with rownames = sample IDs
  
  assign_df <- data.frame(
    sample_id = names(assigned_clusters),
    cluster = factor(assigned_clusters)
  )
  
  return(assign_df)
}


visualise_clusters_new_data <- function(assign_df, deres_list, title=NULL) {
  merged_dataset <- merge_datasets(deres_list)
  response_meta <- merged_dataset$metadata
  response_meta$sample_id <- rownames(response_meta)
  
  merged_response <- merge(response_meta, assign_df, by = "sample_id")
  
  print(table(merged_response$cluster, merged_response$response))
  print(fisher.test(table(merged_response$cluster, merged_response$response)))
  
  df_counts <- merged_response %>%
    count(cluster, response) %>%
    group_by(cluster) %>%
    mutate(
      prop = n / sum(n),
      ypos = cumsum(prop) - prop / 2  # position at center of bar segment
    ) %>%
    ungroup()
  
  ggplot(df_counts, aes(x = cluster, y = prop, fill = response)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), 
              size = 4, color = "black") +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    labs(title=title, y="Proportion of Response", x=NULL)
}

# ----- 3.3 NMF Downstream  ----

limma_cluster_analysis <- function(original_data, sample_clustering) {
  # Get number of clusters
  numeric_vector <- as.numeric(as.character(sample_clustering))
  no_of_clusters <- max(numeric_vector)
  
  limma_res<-list()
  for (i in 1:no_of_clusters) {
    coi <- original_data[,sample_clustering==i]
    rc <- original_data[,sample_clustering!=i]
    
    factors<- factor(rep(c("coi", "others"), c(length(which(sample_clustering==i)),length(which(sample_clustering!=i)))))
    
    mm<- model.matrix(~factors-0)
    colnames(mm) <- levels(factors)
    row.names(mm)<-c(colnames(coi),colnames(rc))
    
    mm[(ncol(coi)+1):(ncol(original_data)),1]<-0
    
    contrast.mat <- makeContrasts(
      Diff = coi - others,
      levels = mm)
    
    colnames(mm)[1]<-paste0("k",i)
    rownames(contrast.mat)[1]<-paste0("k",i)
    
    fit <- lmFit(cbind(coi,rc), mm)
    
    contrast <- contrasts.fit(fit, contrast.mat)
    contrast2<-eBayes(contrast)
    
    output<-topTable(contrast2,number = nrow(original_data))
    output<-list(output)
    
    limma_res<-append(limma_res,output)
  }
  
  # Display top 20 genes in each cluster according to adj.P.val (already adjusted)
  for (i in 1:no_of_clusters) {
    print(glue("\n\nTop 20 DE genes in cluster = {i}"))
    print(rownames(limma_res[[i]][1:20,]))
  }
  
  return(limma_res)
}

limma_print_volcano <- function(limma_res) {
  plot_list <- lapply(seq_along(limma_res[1:4]), function(i) {
    EnhancedVolcano(limma_res[[i]],
                    lab = rownames(limma_res[[i]]),
                    x = 'logFC',
                    y = 'adj.P.Val', pCutoff = 0.05, FCcutoff = 1,
                    subtitle=NULL, title=paste0("Cluster ",i))
  })
  ggarrange(plotlist = plot_list, ncol = 2, nrow = 2) 
}

# Function to perform GSEA analysis
nmf_pathway_analysis <- function(limma_res){
  no_of_clusters <- length(limma_res)
  
  # Process data for GSEA analysis
  gene_lists<-list()
  for (i in 1:no_of_clusters) {
    gene_list<-limma_res[[i]]$logFC
    #Limma outputs log10FC, but clusterProfiler requires log2FC to the values are transformed:
    # gene_list<-log2(10^gene_list)
    names(gene_list)<-rownames(limma_res[[i]])
    
    #The lists are sorted in decreasing FC for clusterProfiler to work
    gene_list<-sort(gene_list, decreasing = TRUE)
    gene_list<-list(gene_list)
    gene_lists<-append(gene_lists,gene_list)
  }
  
  # Perform GSEA analysis
  GSEA_res<-list()
  msigdb_sets <- msigdbr(species = "human", collection = "H")
  msigdbr_t2g <- dplyr::distinct(msigdb_sets, gs_name, gene_symbol)
  
  for (i in 1:no_of_clusters) {
    set.seed(55)
    gse <- gseGO(geneList=gene_lists[[i]], 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 nPerm = 1000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Hs.eg.db", 
                 pAdjustMethod = "none")
    
    # gseaRes <- GSEA(gene_lists[[i]], TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05)
    
    GSEA_res<-append(GSEA_res, gse)
  }
  
  return(GSEA_res)
}


prep_GSEA_NMF <- function(limma_res) {
  no_of_clusters <- length(limma_res)
  
  # Process data for GSEA analysis
  gene_lists<-list()
  for (i in 1:no_of_clusters) {
    gene_list<-limma_res[[i]]$logFC
    names(gene_list)<-rownames(limma_res[[i]])
    
    #The lists are sorted in decreasing FC for clusterProfiler to work
    gene_list<-sort(gene_list, decreasing = TRUE)
    gene_list<-list(gene_list)
    gene_lists<-append(gene_lists,gene_list)
  }
  
  # Perform GSEA analysis
  GSEA_res<-list()
  for (i in 1:no_of_clusters) {
    # prep hallmark pathways
    msigdb_sets_hallmark <- msigdbr(species = "human", collection = "H")
    msigdb_sets_hallmark <- msigdb_sets_hallmark[msigdb_sets_hallmark$gs_name %in% paste0("HALLMARK_", hallmark_pathways$ID), ]
    msigdbr_t2g_hallmark <- dplyr::distinct(msigdb_sets_hallmark, gs_name, gene_symbol)
    
    hallmark_enrich <- GSEA(gene_lists[[i]], TERM2GENE = msigdbr_t2g_hallmark, pvalueCutoff = 2)
    hallmark_enrich <- hallmark_enrich@result
    hallmark_enrich$facet <- "Hallmark"
    hallmark_enrich$ID <- sub("HALLMARK_", "", hallmark_enrich$ID)
    hallmark_enrich <- left_join(hallmark_enrich, hallmark_pathways, by="ID")
    
    # custom enrich
    custom_gene_set <- data.frame(term = matrisome_file$Category, gene = matrisome_file$GeneSymbol)
    custom_enrich <- GSEA(gene_lists[[i]],  TERM2GENE = custom_gene_set, pvalueCutoff = 2)
    custom_enrich <- custom_enrich@result
    custom_enrich$facet <- "Matrisome"
    custom_enrich <- left_join(custom_enrich, matrisome_pathways, by="ID")
  
    
    # combine
    gsea <- rbind(hallmark_enrich, custom_enrich)
    gsea$cluster <- paste0("Cluster ", i)
    
    GSEA_res[[paste0("Cluster_", i)]] <- gsea
  }
  
  # combine all
  combined_gsea_df <- do.call(rbind, GSEA_res) %>% 
    mutate(significance = ifelse(p.adjust < 0.05, "< 0.05", "≥ 0.05"))
  
  gsea_combined <- ggplot(combined_gsea_df, aes(x = ID, y = cluster)) +
    theme_classic() + 
    
    # dotplot
    geom_point(aes(size = significance, color = NES)) +
    scale_size_manual(values = c("≥ 0.05" = 1, "< 0.05" = 5)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    
    # nested facet plot
    ggh4x::facet_nested(cols = vars(facet, category), scales = "free", space = "free") + 
    
    # text theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=7),
      axis.text.y = element_markdown(size = 8),
      strip.text.x = element_text(size = 8)
    ) + labs(x = NULL, y = NULL) 
  
  return(gsea_combined)
}



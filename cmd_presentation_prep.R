source("./Pipeline.R", local = knitr::knit_global())

matrisome_genes <- c("COL11A1", "COMP", "FN1", "VCAN", "CTSB", "COL1A1", "AGT", "ANXA5", "ANXA6", "LAMB1", "FBLN2", "LAMC1", "LGALS3", "CTSG", "HSPG2", "COL15A1", "ANXA1", "LAMA4", "COL6A6", "VWF", "ABI3BP", "TNXB")

# Get command line args
args <- commandArgs(trailingOnly = TRUE)

counts_matrix <- args[1]
metadata      <- args[2]
design        <- args[3]
species       <- if (length(args) >= 4) args[4] else "human"
txi_status    <- if (length(args) >= 5) as.logical(args[5]) else FALSE
annotation    <- if (length(args) >= 6) args[6] else "biomart"
treatment     <- if (length(args) >= 7) args[7] else "untreated"
do_pathway    <- if (length(args) >= 8) as.logical(args[8]) else FALSE
extra_info    <- if (length(args) >= 9) args[9] else NULL

# Run differential expression analysis
de_res <- rna_seq_analysis(counts_matrix, metadata, design, species,
                 txi_status, annotation, treatment)
title <- de_res$title

### Slide 1: Title

### Slide 2: extra info (if any)

### Slide 3: metadata
df_png(metadata, "metadata.png", col_width = 200)

### Slide 3: PCA
visualise_pca(de_res, save=TRUE) # saved as pca.png

### Slide 4: Volcano
visualise_volcano(de_res, save=TRUE) # saved as volcano.png

### Slide 5: Significant Genes
significant_genes <- view_res_table(de_res) %>% dplyr::slice(1:15)
df_png(significant_genes, "significant_genes.png", col_width = 110)

### Slide 6: Heatmap matrisome
visualise_heatmap(de_res, save=TRUE, proteins=matrisome_genes, save_name="heatmap_matrisome.png")

### Slide 7: Heatmap significant
visualise_heatmap(de_res, save=TRUE, save_name="heatmap_significant.png")

### Slide 8: Pathway (as required)
if (do_pathway) df_png(pathway_enrichment(de_res), "pathway_analysis.png", col_width = 200)


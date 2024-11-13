# required packages ####
library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggbreak)
library(ggfortify)
library(mosaic)
library(RColorBrewer)
library(dplyr)
library(EnhancedVolcano)
library(fgsea)
library(qusage)
library(data.table)
library(pheatmap)
library(tidyverse)
library(clustree)
library(celldex)
library(SingleR)
library(DoubletFinder)
library(harmony)
library(anndata)
library(chisq.posthoc.test)
library(presto)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
library(reticulate)
library(ComplexHeatmap)
library(monocle3)
library(SeuratWrappers)
library(magrittr)
library(infercnv)
library(sctransform)
library(glmGamPoi)
library(scCustomize)
library(parallel)
library(BiocParallel)
library(fpc)
library(genomicInstability)
library(ineq)
library(plyr)
library(cowplot)
library(sceasy)
library(dendsort)
library(liana)
library(lmerTest)
library(emmeans)
library(rstatix)
library(expm)
library(Regmex)
library(decoupleR)
library(gprofiler2)

# database of color vectors ####
treatment_colors = c("dodgerblue2", "orange")
sample_colors = c("#E41A1C", "#f58231", "#ffe119", "#3cb44b", "#42d4f4", "#911eb4")

all_broad_colors_alpha = c("#911eb4", "#3cb44b", "#f58231", "#E41A1C", "#ffe119", "#42d4f4")
all_broad_order_alpha = c("B", "DC", "Macrophage", "Monocyte", "Neutrophil", "T/NK")

all_broad_colors_rainbow = c("#E41A1C", "#f58231", "#ffe119", "#3cb44b", "#42d4f4", "#911eb4")
all_broad_order_rainbow = c("Monocyte", "Macrophage", "Neutrophil", "DC", "T/NK", "B")

myeloid_fine_colors_alpha = c("#6a822e", "#024d0d", 
                              "#ffc296", "#f58231", "#7d4722",
                              "#E41A1C", "#ffe119", 
                              "#bfef45", "#3cb44b")
myeloid_fine_order_alpha = c("cDC1", "cDC2", 
                             "M0", "M1-like", "M2-like", 
                             "Monocyte",  "Neutrophil", 
                             "pDC", "pre-cDC")

myeloid_fine_colors_rainbow = c("#E41A1C", "#f58231", "#ffc296", "#ffe119", 
                                "#bfef45", "#3cb44b", "#024d0d", "#6a822e", "#7d4722")
myeloid_fine_order_rainbow = c("Monocyte", "M1-like", "M0", "Neutrophil", 
                               "pDC", "pre-cDC", "cDC2", "cDC1", "M2-like")


lymphoid_fine_colors_alpha = c("#911eb4", "#000075", "#f032e6", "#dcbeff", 
                               "#afebfa", "#86b1f7", "#42d4f4", "#4363d8")
lymphoid_fine_order_alpha = c("B", "CTL", "NK", "Plasma", 
                              "Tcyc", "Tex", "Tn", "Treg")

lymphoid_fine_colors_rainbow = c("#f032e6", "#911eb4", "#dcbeff", "#4363d8", 
                                 "#000075", "#86b1f7", "#42d4f4", "#afebfa")
lymphoid_fine_order_rainbow = c("NK", "B", "Plasma", "Treg", 
                                "CTL", "Tex", "Tn", "Tcyc")


all_fine_colors_alpha = c("#911eb4", "#6a822e", "#024d0d", "#000075",
                          "#ffc296", "#f58231", "#7d4722",
                          "#E41A1C", "#ffe119", "#f032e6", 
                          "#bfef45", "#dcbeff", "#3cb44b", 
                          "#afebfa", "#86b1f7", "#42d4f4", "#4363d8")
all_fine_order_alpha = c("B", "cDC1", "cDC2", "CTL",
                         "M0", "M1-like", "M2-like", 
                         "Monocyte",  "Neutrophil", "NK", 
                         "pDC", "Plasma", "pre-cDC", 
                         "Tcyc", "Tex", "Tn", "Treg")

all_fine_colors_rainbow = c("#E41A1C", "#f58231", "#ffc296", "#ffe119", 
                            "#bfef45", "#3cb44b", "#024d0d", "#6a822e",
                            "#7d4722", "#f032e6", "#911eb4", "#dcbeff", "#4363d8", 
                            "#000075", "#86b1f7", "#42d4f4", "#afebfa")

all_fine_order_rainbow = c("Monocyte", "M1-like", "M0", "Neutrophil", 
                           "pDC", "pre-cDC", "cDC2", "cDC1", 
                           "M2-like", "NK", "B", "Plasma", "Treg",
                           "CTL", "Tex", "Tn", "Tcyc")

# preprocessing functions ####
process_and_doublet_finder <- function(obj, nExp_multiplier) {
  # objects need to be processed through standard workflow before calculating doublets
  obj_dd <- subset(obj, subset = nFeature_RNA > 50 & nFeature_RNA < 10000)
  obj_dd <- NormalizeData(obj_dd)
  obj_dd <- FindVariableFeatures(obj_dd, selection.method = "vst", nfeatures = 2000)
  obj_dd <- ScaleData(obj_dd)
  obj_dd <- RunPCA(obj_dd, npcs = 25)
  obj_dd <- FindNeighbors(obj_dd, dims = 1:15)
  obj_dd <- FindClusters(obj_dd, resolution = 0.5)
  obj_dd <- RunUMAP(obj_dd, dims = 1:15)
  
  # calculate pK
  sweep.res.list <- paramSweep(obj_dd, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # determine the optimal pK value
  optimal_pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
  print(paste("Optimal pK value for current object:", optimal_pK))
  
  # calculate homotypic proportion and expected number of doublets
  annotations <- obj_dd@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(nExp_multiplier * nrow(obj_dd@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # apply DoubletFinder
  obj_dd <- doubletFinder(obj_dd, 
                          PCs = 1:15, 
                          pN = 0.25, 
                          pK = as.numeric(optimal_pK),
                          nExp = nExp_poi.adj, 
                          reuse.pANN = FALSE, 
                          sct = FALSE)
  
  return(obj_dd)
}
doublet_detected_objs <- lapply(seq_along(raw_objs), function(x) {
  obj <- raw_objs[[x]]$object
  expected_rate <- raw_objs[[x]]$expected_rate
  cat(paste("Processing object", x, "of", length(raw_objs), "\n"))
  processed_obj <- process_and_doublet_finder(obj, expected_rate)
  cat(paste("Finished processing object", x, "with optimal pK value\n\n"))
  return(processed_obj)
})
doublet_dimplot_check <- function(doublet_detected_objs) {
  plot_list <- list()
  
  for (obj_name in names(doublet_detected_objs)) {
    obj <- doublet_detected_objs[[obj_name]]
    
    # find the doublet column
    group_by_column <- grep("^DF.classifications_", colnames(obj@meta.data), value = TRUE)
    
    # generate DimPlot if the column is found
    if (length(group_by_column) > 0) {
      plot <- scCustomize::DimPlot_scCustom(obj, reduction = "umap", group.by = group_by_column[1], colors_use = c("red", "grey80"))
      plot_list[[obj_name]] <- plot
    } else {
      stop("No column starting with 'DF.classifications_' found in object: ", obj_name)
    }
  }
  
  # combine all plots using cowplot
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 3)
  
  return(combined_plot)
}
print_metadata_table <- function(seurat_obj, column_prefix) {
  # find the column name that starts with the given prefix
  column_name <- grep(column_prefix, colnames(seurat_obj@meta.data), value = TRUE)
  
  if (length(column_name) > 0) {
    for (col in column_name) {
      cat("Table for column:", col, "\n")
      print(table(seurat_obj@meta.data[[col]]))
      cat("\n")
    }
  } else {
    print(paste("No columns starting with", column_prefix, "found in metadata"))
  }
}
update_doublets <- function(raw, doublet) {
  for (obj_name in names(raw)) {
    raw_obj <- raw[[obj_name]]$object
    doublet_obj <- doublet[[obj_name]]
    
    if (!is.null(raw_obj) && !is.null(doublet_obj)) {
      # find the correct DF.classifications column in doublet_obj
      classification_col <- grep("^DF.classifications", colnames(doublet_obj@meta.data), value = TRUE)
      
      if (length(classification_col) > 0) {
        classification_col <- classification_col[1]
        
        # subset the doublet detected metadata
        doublet_cells <- rownames(subset(doublet_obj@meta.data, doublet_obj@meta.data[[classification_col]] == "Doublet"))
        
        # update the doublet column in the raw object metadata
        raw_obj@meta.data$doublet = paste("single")
        raw_obj@meta.data$doublet[rownames(raw_obj@meta.data) %in% doublet_cells] <- "doublet"
        
        # assign the updated raw object back to the list
        raw_objs[[obj_name]]$object <- raw_obj
      }
    }
  }
  
  return(raw_objs)
}
qc_metrics <- function(obj){
  obj <- scCustomize::Add_Mito_Ribo(object = obj, species = "Mouse")
  obj <- scCustomize::Add_Cell_Complexity(object = obj)
  
  return(obj)
}
qc_plots <- function(obj, sample) {
  # create plots
  p1 <- scCustomize::QC_Plots_Genes(seurat_object = obj, low_cutoff = 100, high_cutoff = 8500)
  p2 <- scCustomize::QC_Plots_UMIs(seurat_object = obj, low_cutoff = 500, high_cutoff = 50000)
  p3 <- scCustomize::QC_Plots_Mito(seurat_object = obj, high_cutoff = 20)
  p4 <- scCustomize::QC_Plots_Complexity(seurat_object = obj, high_cutoff = 0.75)
  
  qc_plots <- wrap_plots(p1, p2, p3, p4, ncol = 4)
  
  p1 <- scCustomize::QC_Plot_UMIvsGene(seurat_object = obj, 
                                       low_cutoff_gene = 100, high_cutoff_gene = 8500, 
                                       low_cutoff_UMI = 500, high_cutoff_UMI = 50000)
  p2 <- scCustomize::QC_Plot_GenevsFeature(seurat_object = obj, feature1 = "percent_mito", 
                                           low_cutoff_gene = 100, high_cutoff_gene = 8500, 
                                           high_cutoff_feature = 20)
  
  umi_gene_plots <- wrap_plots(p1, p2)
  
  # combine all plots 
  combined_plots <- wrap_plots(qc_plots, umi_gene_plots)
  
  # display plots in the RStudio plots window
  print(combined_plots + plot_annotation(title = paste("QC Plots:", sample)), ncol = 1)
}
# annotation functions ####
do_clustree <- function(obj, res_vector = NULL, res_seq = NULL){
  # set resolutions to use
  if(is.null(res_vector)){
    if(!is.null(res_seq)){
      res_vector = seq(res_seq[1], res_seq[2], res_seq[3])
    } else stop("Please input res_vector or res_seq parameters")
  }
  
  # perform multiple clusterings at resolutions specified
  for (i in res_vector){
    obj <- FindClusters(obj, resolution = i)
  }
  
  # plot clustree graph
  clustree::clustree(obj@meta.data, prefix = "SCT_snn_res.")
}
# analysis functions ####
combined_volcano_plot <- function(data_list, x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, FCcutoff = 1, col = c('black', 'black', 'black', 'red3'), ncol = 1) {
  plots <- lapply(names(data_list), function(name) {
    data <- data_list[[name]]
    EnhancedVolcano(data, lab = data$X, x = x, y = y, pCutoff = pCutoff, FCcutoff = FCcutoff, col = col, title = name, legendPosition = 'none', 
                    labSize = 4, axisLabSize = 10, titleLabSize = 12, subtitleLabSize = 10, captionLabSize = 6)
  })
  
  combined_plot <- plot_grid(plotlist = plots, ncol = ncol)
  
  return(combined_plot)
}
create_gsea_vectors <- function(data_list) {
  gsea_data_list <- lapply(data_list, function(data) {
    gsea_data <- setNames(data$avg_log2FC, data$X)
    gsea_data <- sort(gsea_data, decreasing = TRUE)
    return(gsea_data)
  })
  
  return(gsea_data_list)
}
run_through_fgsea <- function(msigdb_list, named_degs_list, minSize = 10, maxSize = 500, nPermSimple = 50000, gseaParam = 0.5) {
  # initialize empty results list
  fgsea_results <- list()
  
  # iterate over msigdb and degs
  for (i in seq_along(msigdb_list)) {
    msigdb_name <- names(msigdb_list)[i]
    for (j in seq_along(named_degs_list)) {
      degs_name <- names(named_degs_list)[j]
      
      # print progress message since this takes a little bit of time
      message(paste("Processing:", degs_name, "through", msigdb_name))
      
      # running fgsea
      fgsea_res <- fgsea(pathways = msigdb_list[[i]], 
                         stats    = named_degs_list[[j]],
                         minSize  = minSize,
                         maxSize  = maxSize, 
                         nPermSimple = nPermSimple)
      fgsea_res$Enrichment = ifelse(fgsea_res$NES > 0, "Up-regulated", "Down-regulated")
      fgsea_results[[paste0(msigdb_name, "_", degs_name)]] <- fgsea_res
    }
  }
  return(fgsea_results)
}
# figure functions ####
plotEnrichment_customcolor <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2, custom_color) {
  pd <- plotEnrichmentData(pathway = pathway, stats = stats, 
                           gseaParam = gseaParam)
  with(pd, ggplot(data = curve) + geom_line(aes(x = rank, y = ES), linewidth = 2,
                                            color = custom_color) + geom_segment(data = ticks, mapping = aes(x = rank, 
                                                                                                             y = -spreadES/16, xend = rank, yend = spreadES/16), linewidth = ticksSize) + 
         geom_hline(yintercept = posES, colour = "red", linetype = "dashed") + 
         geom_hline(yintercept = negES, colour = "red", linetype = "dashed") + 
         geom_hline(yintercept = 0, colour = "black") + theme(panel.background = element_blank(), 
                                                              panel.grid.major = element_line(color = "grey92")) + 
         labs(x = "rank", y = "enrichment score"))
} 
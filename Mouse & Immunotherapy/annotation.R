# 3.0 load objects ####
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
options(Seurat.object.assay.version = "v5")
basefile_path <- "/path/to/files" # please download files and update this line to point to your directory
data_obj <- readRDS(file=file.path(save_path, "corin_pd1_obj.rds"))

# 3.1 initial preprocessing ####
data_obj <- SCTransform(data_obj)
data_obj <- RunPCA(data_obj)
data_obj <- FindNeighbors(data_obj, reduction = "pca", dims = 1:30)
do_clustree(obj = data_obj, res_seq = c(0.1, 1, 0.1))
data_obj <- FindClusters(data_obj, resolution = 0.6, cluster.name = "allcell_clusters")
data_obj <- RunUMAP(data_obj, dims = 1:30, reduction = "pca", reduction.name = "umap")

# cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mmus_s = gprofiler2::gorth(s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gprofiler2::gorth(g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

data_obj <- CellCycleScoring(data_obj,
                             s.features = mmus_s, g2m.features = mmus_g2m, 
                             set.ident = TRUE, assay = 'SCT')

p1 <- DimPlot_scCustom(data_obj, reduction = "umap", group.by = "allcell_clusters")+NoLegend()
p2 <- DimPlot_scCustom(data_obj, reduction = "umap", group.by = "sample", colors_use = sample_colors)
p3 <- DimPlot_scCustom(data_obj, reduction = "umap", group.by = "treatment", colors_use = treatment_colors)
p4 <- FeaturePlot_scCustom(data_obj, features = "S.Score", reduction = "umap")
p5 <- FeaturePlot_scCustom(data_obj, features = "G2M.Score", reduction = "umap")
p6 <- FeaturePlot_scCustom(data_obj, features = "percent_mito", reduction = "umap")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

# checking identities
FeaturePlot_scCustom(data_obj, features = c("Sox10", "Dct", "Mlana", # melanoma
                                            "Pecam1", "Vwf", "Cdh5", # endothelial 
                                            "Ptprc", "Cd33", "Csf1r", # general immune/myeloid
                                            "Fap", "Col1a1", "Pdpn"), # fibroblast
                     reduction = "umap")

# removing everything that is not immune for this analysis
imm_obj <- subset(x = data_obj, subset = allcell_clusters %in% c("6", "22", "13", "24"), invert = T)

# 3.2 reprocessing only immune cells ====
imm_obj <- SCTransform(imm_obj)
imm_obj <- RunPCA(imm_obj)
imm_obj <- FindNeighbors(imm_obj, dims = 1:20)

# finding right resolution
do_clustree(obj = imm_obj, res_seq = c(0.1, 1, 0.1)) # function in 3.1
do_clustree(obj = imm_obj, res_seq = c(0.6, 0.7, 0.02)) # function in 3.1

imm_obj <- FindClusters(imm_obj, resolution = 0.66, cluster.name = "immune_clusters") 
imm_obj <- RunUMAP(imm_obj, dims = 1:20, reduction.name = "umap.imm")

p1 <- DimPlot_scCustom(imm_obj, reduction = "umap.imm", group.by = "immune_clusters")+NoLegend()
p2 <- DimPlot_scCustom(imm_obj, reduction = "umap.imm", group.by = "sample", colors_use = sample_colors)
p3 <- DimPlot_scCustom(imm_obj, reduction = "umap.imm", group.by = "treatment", colors_use = treatment_colors)
p4 <- FeaturePlot_scCustom(imm_obj, features = "S.Score", reduction = "umap.imm")
p5 <- FeaturePlot_scCustom(imm_obj, features = "G2M.Score", reduction = "umap.imm")
p6 <- FeaturePlot_scCustom(imm_obj, features = "percent_mito", reduction = "umap.imm")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

# 3.3 coarse annotations ####
FeaturePlot_scCustom(imm_obj, features = c("Sox10", "Dct", "Mlana", # melanoma
                                           "Pecam1", "Vwf", "Cdh5", # endothelial 
                                           "Ptprc", "Cd33", "Csf1r", # general immune/myeloid
                                           "Fap", "Col1a1", "Pdpn"), # fibroblast
                     reduction = "umap.imm") # looks great

FeaturePlot_scCustom(imm_obj, features = c("Itgam", # myeloid
                                           "Cd14", #monocyte
                                           "Cd68", "Adgre1", "Mrc1", #macrophage
                                           "Itgax", "H2-Ab1", "Cd86", # dendritic
                                           "Ly6g", "S100a8"), #neutrophil
                     reduction = "umap.imm")

FeaturePlot_scCustom(imm_obj, features = c("Cd3e", "Cd4", "Cd8a", # t
                                           "Cd19", "Ms4a1", "Cd79a", # b
                                           "Ncr1", "Klrb1c", # nk
                                           "Sdc1", "Tnfrsf17"), # plasma
                     reduction = "umap.imm")

# helpful throughout annotations: https://www.cellsignal.com/pathways/immune-cell-markers-mouse

imm_obj@meta.data$coarse_anno[imm_obj@meta.data$immune_clusters %in% c("0", "1", "4", "8", "5",
                                                                       "21", "14", "9", "6", "15",
                                                                       "22", "3")] <- "myeloid"

imm_obj@meta.data$coarse_anno[imm_obj@meta.data$immune_clusters %in% c("20", "12", "7", "16", "19",
                                                                       "17", "11", "10", "13", "23",
                                                                       "2", "18")] <- "lymphoid"

lymphoid_obj <- subset(x = imm_obj, subset = coarse_anno =="lymphoid")
myeloid_obj <- subset(x = imm_obj, subset = coarse_anno =="myeloid")

# 3.4 lymphoid fine annotation ====
lymphoid_obj <- SCTransform(lymphoid_obj)
lymphoid_obj <- RunPCA(lymphoid_obj)
lymphoid_obj <- FindNeighbors(lymphoid_obj, dims = 1:20)
do_clustree(obj = lymphoid_obj, res_seq = c(0.1, 1, 0.1)) # function in 3.1
lymphoid_obj <- FindClusters(lymphoid_obj, resolution = 0.4, cluster.name = "lymphoid_clusters") 
lymphoid_obj <- RunUMAP(lymphoid_obj, dims = 1:20, reduction.name = "umap.lymph")

p1 <- DimPlot_scCustom(lymphoid_obj, reduction = "umap.lymph", group.by = "lymphoid_clusters")+NoLegend()
p2 <- DimPlot_scCustom(lymphoid_obj, reduction = "umap.lymph", group.by = "sample", colors_use = sample_colors)
p3 <- DimPlot_scCustom(lymphoid_obj, reduction = "umap.lymph", group.by = "treatment", colors_use = treatment_colors)
p4 <- FeaturePlot_scCustom(lymphoid_obj, features = "S.Score", reduction = "umap.lymph")
p5 <- FeaturePlot_scCustom(lymphoid_obj, features = "G2M.Score", reduction = "umap.lymph")
p6 <- FeaturePlot_scCustom(lymphoid_obj, features = "percent_mito", reduction = "umap.lymph")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

lymphoid_obj <- PrepSCTFindMarkers(lymphoid_obj)
lymphoid.markers <- FindAllMarkers(lymphoid_obj, only.pos = TRUE,
                                   min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

write.csv(lymphoid.markers, 
          file=file.path(save_path, "lymphoid_markers_MAST.csv"), 
          row.names=TRUE)

lymphoid.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(lymphoid_obj, features = top5$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

lymphoid.markers %>% 
  dplyr::filter(cluster == 0) %>% # Change this to explore markers of other clusters
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n=20) %>%
  dplyr::pull(gene)

# LYMPHOID ANNOTATIONS
# 0 - b cell (Cd79a, Ms4a1, Cd19, Bank1)
# 1 -  possibly doublet contam
# 2 - cd8 t cell (gzmb, cd8a, cd8b1, cd3)
# 3 - nk cell (no cd3, klrk1, klrb1c, Il2rb)
# 4 - naive t cell (lef1, tcf7, Skap1, and Prkcq) 
# 5 - b cells (ebf1, bank1, cd79a, pax5. bach2)
# 6 - cd4 treg (ctla4, foxp3, hif1a, tnfrsf4, icos, cd4)
# 7 - b cell (ebf1, bank1, bach2, pax5)
# 8 - cycling t cells (top2a, mki67, tubb5)
# 9 - exhausted cd8 t cells (tox, havcr2)
# 10 - plasma cell (jchain, mzb1)
# 11 - possibly doublet contam
# 12 - stressed

lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("0", "5", "7")] <- "B"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("1")] <- "test_1"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("2")] <- "CTL"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("3")] <- "NK"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("4")] <- "Tn"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("6")] <- "Treg"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("8")] <- "Tcyc"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("9")] <- "Tex"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("10")] <- "Plasma"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("11")] <- "test_11"
lymphoid_obj[[]]$fine_anno[lymphoid_obj[[]]$lymphoid_clusters %in% c("12")] <- "test_12"

# adding labels to main object
for (category in names(table(lymphoid_obj@meta.data$fine_anno))) {
  imm_obj@meta.data$fine_anno[imm_obj@meta.data$barcode %in% subset(lymphoid_obj@meta.data, fine_anno == category)$barcode] <- category
}

DimPlot_scCustom(imm_obj, group.by = "fine_anno", 
                 reduction = "umap.imm", figure_plot = T)

# 3.5 myeloid fine annotation ====
myeloid_obj <- SCTransform(myeloid_obj)
myeloid_obj <- RunPCA(myeloid_obj)
myeloid_obj <- FindNeighbors(myeloid_obj, dims = 1:20)
do_clustree(obj = myeloid_obj, res_seq = c(0.1, 1, 0.1)) # function in 3.1
myeloid_obj <- FindClusters(myeloid_obj, resolution = 0.6, cluster.name = "myeloid_clusters") 
myeloid_obj <- RunUMAP(myeloid_obj, dims = 1:20, reduction.name = "umap.myel")

p1 <- DimPlot_scCustom(myeloid_obj, reduction = "umap.myel", group.by = "myeloid_clusters")+NoLegend()
p2 <- DimPlot_scCustom(myeloid_obj, reduction = "umap.myel", group.by = "sample", colors_use = sample_colors)
p3 <- DimPlot_scCustom(myeloid_obj, reduction = "umap.myel", group.by = "treatment", colors_use = treatment_colors)
p4 <- FeaturePlot_scCustom(myeloid_obj, features = "S.Score", reduction = "umap.myel")
p5 <- FeaturePlot_scCustom(myeloid_obj, features = "G2M.Score", reduction = "umap.myel")
p6 <- FeaturePlot_scCustom(myeloid_obj, features = "percent_mito", reduction = "umap.myel")

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

myeloid_obj <- PrepSCTFindMarkers(myeloid_obj)
myeloid.markers <- FindAllMarkers(myeloid_obj, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

write.csv(myeloid.markers, 
          file=file.path(save_path, "myeloid_markers_MAST.csv"), 
          row.names=TRUE)

myeloid.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(myeloid_obj, features = top5$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

myeloid.markers %>% 
  dplyr::filter(cluster == 0) %>% # Change this to explore markers of other clusters
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n=20) %>%
  dplyr::pull(gene)

FeaturePlot_scCustom(myeloid_obj, features = c("Cd14", #mono
                                               "Cd83", # activation
                                               "Adgre1", "Cd68", #mac
                                               "Cd86", "Cd80", "Nos2", #m1
                                               "Arg1", "Cd163", "Mrc1"), #m2
                     reduction = "umap.myel")

FeaturePlot_scCustom(myeloid_obj, features = c("Itgax", "H2-Ab1", #cDC
                                               "Clec9a", "Xcr1", #cDC1
                                               "Cd8a", #resident DC
                                               "Itgae", #migratory DC
                                               "Itgam", "Sirpa", #cDC2
                                               "Siglech", "Bst2"), #pDC
                     reduction = "umap.myel")

FeaturePlot_scCustom(myeloid_obj, features = c("Ly6c1", "Ly6g", "Arg1", #MDSC - combos
                                               "Ly6g", "Adgre1", #neutrophil
                                               "Ccr3", "Siglecf", #eosinophil
                                               "Fcer1a", "Kit", "Fcer2a"), #basophil and mast
                     reduction = "umap.myel")

FeaturePlot_scCustom(myeloid_obj, 
                     reduction = "umap.myel", 
                     features = c("Cd86", "Cd80", "Cd68", "H2-Ab1", "Il1r1", 
                                  "Tlr2", "Tlr4","Nos2", "Socs3")) # m1 like

FeaturePlot_scCustom(myeloid_obj, 
                     reduction = "umap.myel", 
                     features = c("Cd163", "H2-Ab1", "Msr1", "Mrc1", 
                                  "Chil3", "Retnla", "Arg1", "Vegfa", "Tlr8")) # m2 like

FeaturePlot_scCustom(myeloid_obj, 
                     reduction = "umap.myel", 
                     features = c("Ccl2", "Cd3e", "Pecam1", "Il10", "Pcna", "Vegfa")) # tam

FeaturePlot_scCustom(myeloid_obj, 
                     reduction = "umap.myel", 
                     features = c("Fyn", "Lat", "Lck", "Zap70")) # tcr

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6476302/

# MYELOID ANNOTATIONS
# 0 - M2c (gas6, mrc1, selenop, timp2, Lyve1)
# 1 - possibly doublet contam
# 2 - m0 (Ms4a7, C1qa, C1qb, C1qc, Cd68)
# 3 - m1 (Il31ra, Hexb high perc mito aka highly metabolic)
# 4 - neutrophil (Cxcr2, S100a8, S100a9, Il1b, and Retnlg)
# 5 -  dendritic cells (lots of HLA)
# 6 - monocyte (Ccr2, Lyz2, Nr4a1)
# 7 -  m1 (nos2, Zeb2, Ptprj, Gab2, Arid5b)
# 8 - m2a (Mrc1 (CD206), Arg1, Fn1, and Ccl24)
# 9 - m2d (Vegfa, Gpnmb, Lgals3, and Spp1)
# 10 - neutrophil (Ly6g, S100a8, S100a9, Cxcr2, and Cxcl3)
# 11 - pDC (Siglech, Irf8, Tcf4, Ccr9, and Bcl11a)
# 12 - cDC1 (Xcr1, Clec9a, Irf8, H2-Ab1, and Cd24a)
# 13 - cDC2 (Ccr7, Ccl22, Relb, and Ly75)

myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("0")] <- "M2c"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("1")] <- "test_mye_1"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("2")] <- "M0"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("3")] <- "M1"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("4")] <- "Neutrophil"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("5")] <- "DC"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("6")] <- "Monocyte"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("7")] <- "M1"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("8")] <- "M2a"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("9")] <- "M2d"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("10")] <- "Neutrophil"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("11")] <- "pDC"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("12")] <- "cDC1"
myeloid_obj[[]]$fine_anno[myeloid_obj[[]]$myeloid_clusters %in% c("13")] <- "cDC2"

# adding labels to main object
for (category in names(table(myeloid_obj@meta.data$fine_anno))) {
  imm_obj@meta.data$fine_anno[imm_obj@meta.data$barcode %in% subset(myeloid_obj@meta.data, fine_anno == category)$barcode] <- category
}

DimPlot_scCustom(imm_obj, group.by = "fine_anno", 
                 reduction = "umap.imm", figure_plot = T)

# 3.6 cleaning up data ====
# some cells are very clearly doublets that were not properly removed - going to subset and recluster
VlnPlot_scCustom(imm_obj, features = c("percent_mito", "S.Score", "G2M.Score", 
                                       "log10GenesPerUMI", "nCount_RNA", "nFeature_RNA"), group.by = "fine_anno", pt.size = 0)

imm_obj_clean <- subset(x = imm_obj, subset = fine_anno %in% c("test_1", "test_11", "test_mye_1", "test_12"), invert = T)

imm_obj_clean <- SCTransform(imm_obj_clean)
imm_obj_clean <- RunPCA(imm_obj_clean)
imm_obj_clean <- FindNeighbors(imm_obj_clean, dims = 1:20)

# finding right resolution
do_clustree(obj = imm_obj_clean, res_seq = c(0.1, 1, 0.1)) # function in 3.1

imm_obj_clean <- FindClusters(imm_obj_clean, resolution = 0.8, cluster.name = "immune_clean_clusters") 
imm_obj_clean <- RunUMAP(imm_obj_clean, dims = 1:20, reduction.name = "umap.imm.clean")

# seeing how well clusters overlap
table(imm_obj_clean@meta.data$immune_clean_clusters, imm_obj_clean@meta.data$fine_anno)
imm_obj_clean@meta.data$finalized_fine_anno = imm_obj_clean@meta.data$fine_anno
imm_obj_clean@meta.data$finalized_fine_anno[imm_obj_clean@meta.data$finalized_fine_anno %in% c("M1")] <- "M1-like"
imm_obj_clean@meta.data$finalized_fine_anno[imm_obj_clean@meta.data$finalized_fine_anno %in% c("M2a", "M2c", "M2d")] <- "M2-like"

# t cells only
t_obj = subset(x=imm_obj_clean, subset = finalized_fine_anno %in% c("Treg", "CTL", "Tex", "Tn", "Tcyc"))
t_obj <- SCTransform(t_obj)
t_obj <- RunPCA(t_obj)
t_obj <- FindNeighbors(t_obj, dims = 1:20)

# finding right resolution
t_obj = subset(imm_obj_clean, subset = fine_anno %in% c("Treg", "CTL", "Tex", "Tn", "Tcyc"))
do_clustree(obj = t_obj, res_seq = c(0.1, 1, 0.1)) # function in 3.1
t_obj <- FindClusters(t_obj, resolution = 0.15, cluster.name = "t_clusters") 
t_obj <- RunUMAP(t_obj, dims = 1:20, reduction.name = "umap.t")

DimPlot_scCustom(t_obj, split.by = "treatment", group.by = "finalized_fine_anno", reduction = "umap.t",
                 colors_use = c(lymphoid_fine_colors_alpha[2], lymphoid_fine_colors_alpha[5:8]))

#removing b cell contam
t_obj_clean = subset(t_obj, subset = t_clusters == 4, invert = T)
t_obj_clean <- SCTransform(t_obj_clean)
t_obj_clean <- RunPCA(t_obj_clean)
t_obj_clean <- FindNeighbors(t_obj_clean, dims = 1:20)
t_obj_clean <- FindClusters(t_obj_clean, resolution = 0.1, cluster.name = "t_clean_clusters") 
t_obj_clean <- RunUMAP(t_obj_clean, dims = 1:20, reduction.name = "umap.t.clean")

t_obj_clean <- FindClusters(t_obj_clean, resolution = 2, cluster.name = "t_clean_clusters_xtreme") 

t_obj_cleaned = subset(t_obj_clean, subset = t_clean_clusters_xtreme == 17, invert = T)
t_obj_cleaned <- SCTransform(t_obj_cleaned)
t_obj_cleaned <- RunPCA(t_obj_cleaned)
t_obj_cleaned <- FindNeighbors(t_obj_cleaned, dims = 1:20)
t_obj_cleaned <- FindClusters(t_obj_cleaned, resolution = 0.1, cluster.name = "t_cleaned_clusters") 
t_obj_cleaned <- RunUMAP(t_obj_cleaned, dims = 1:20, reduction.name = "umap.t.cleaned")

DimPlot_scCustom(t_obj_cleaned, split.by = "treatment", group.by = "finalized_fine_anno", reduction = "umap.t.cleaned",
                 colors_use = c(lymphoid_fine_colors_alpha[2], lymphoid_fine_colors_alpha[5:8]))

# removing the t/b from main object
imm_obj_cleaned <- subset(x = imm_obj_clean, subset = barcode %in% subset(t_obj@meta.data, t_clusters ==4)$barcode, invert = T)
imm_obj_cleaned <- subset(x = imm_obj_cleaned, subset = barcode %in% subset(t_obj_clean@meta.data, t_clean_clusters_xtreme ==17)$barcode, invert = T)

imm_obj_cleaned <- SCTransform(imm_obj_cleaned)
imm_obj_cleaned <- RunPCA(imm_obj_cleaned)
imm_obj_cleaned <- FindNeighbors(imm_obj_cleaned, dims = 1:20)

# finding right resolution
do_clustree(obj = imm_obj_cleaned, res_seq = c(0.1, 1, 0.1)) # function in 3.1

imm_obj_cleaned <- FindClusters(imm_obj_cleaned, resolution = 0.9, cluster.name = "immune_cleaned_clusters") 
imm_obj_cleaned <- RunUMAP(imm_obj_cleaned, dims = 1:20, reduction.name = "umap.imm.cleaned")

p1 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "immune_cleaned_clusters")+NoLegend()
p2 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "sample", colors_use = sample_colors)
p3 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "treatment", colors_use = treatment_colors, shuffle = T)
p4 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "coarse_anno")
p5 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "finalized_fine_anno")
p6 <- FeaturePlot_scCustom(imm_obj_cleaned, features = "S.Score", reduction = "umap.imm.cleaned")
p7 <- FeaturePlot_scCustom(imm_obj_cleaned, features = "G2M.Score", reduction = "umap.imm.cleaned")
p8 <- FeaturePlot_scCustom(imm_obj_cleaned, features = "percent_mito", reduction = "umap.imm.cleaned")
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4)

p1 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "broad_anno", colors_use = all_broad_colors_alpha)
p2 <- DimPlot_scCustom(imm_obj_cleaned, reduction = "umap.imm.cleaned", group.by = "finalized_fine_anno", colors_use = all_fine_colors_alpha)
plot_grid(p1, p2)

# 3.7 saving objects ####
saveRDS(data_obj, file = file.path(save_path, "corin_pd1_obj.rds"))
saveRDS(imm_obj, file = file.path(save_path, "corin_pd1_imm_obj.rds")) 
saveRDS(lymphoid_obj, file = file.path(save_path, "corin_pd1_lymphoid_obj.rds")) 
saveRDS(myeloid_obj, file = file.path(save_path, "corin_pd1_myeloid_obj.rds")) 
saveRDS(imm_obj_clean, file = file.path(save_path, "corin_pd1_imm_obj_clean.rds")) 
saveRDS(imm_obj_cleaned, file = file.path(save_path, "corin_pd1_imm_obj_cleaned.rds")) 
saveRDS(t_obj_cleaned, file = file.path(save_path, "corin_pd1_t_obj_cleaned.rds")) 

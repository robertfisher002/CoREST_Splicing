# 5.0 load data needed ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
options(Seurat.object.assay.version = "v5")
basefile_path <- "/path/to/files" # please download files and update this line to point to your directory
imm_obj_clean <- readRDS(file=file.path(save_path, "corin_pd1_imm_obj_cleaned.rds"))

imm_obj_clean <- readRDS(file=file.path(save_path, "corin_pd1_imm_obj_cleaned.rds"))
t_obj_clean <- readRDS(file=file.path(save_path, "corin_pd1_t_obj_cleaned.rds"))

tcell_degs = read.csv(file=file.path(save_path, "tcell_pd1corin_v_pd1_corrected.csv"))
fgsea_results = readRDS(file=file.path(save_path, "fgsea_results_tcells_pd1corin_v_pd1.rds"))

imm.markers = read.csv(file=file.path(save_path, "all_immune_markers_MAST.csv"))

# please export gmt files from msigdb: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
all_msigdb = read.gmt("path_to/msigdb.v2024.1.Mm.symbols.gmt")
hallmarks = read.gmt("path_to/mh.all.v2024.1.Mm.symbols.gmt")
gobp = read.gmt("path_to/m5.go.bp.v2024.1.Mm.symbols.gmt")
cp = read.gmt("path_to/m2.cp.v2024.1.Mm.symbols.gmt")

msigdb_of_interest = list(all = all_msigdb, hallmarks = hallmarks, 
                          gobp = gobp, cp = cp)
# 5.1 main umap figures ====
Idents(object = imm_obj_clean) <- "finalized_fine_anno"
imm_obj_clean@meta.data$finalized_fine_anno <- factor(imm_obj_clean@meta.data$finalized_fine_anno, levels = c(all_fine_order_rainbow))
DimPlot_scCustom(imm_obj_clean, colors_use = all_fine_colors_rainbow, reduction = "umap.imm.cleaned",
                 shuffle = T, pt.size = 0.01, group.by="finalized_fine_anno", figure_plot = T)
ggsave(filename = file.path(save_path, "mainUMAP_by_celltype.pdf"), height = 8, width = 11)
ggsave(filename = file.path(save_path, "mainUMAP_by_celltype_small.pdf"), height = 4, width = 5.5)

DimPlot_scCustom(imm_obj_clean, colors_use = c("grey30", "#c778ff"), reduction = "umap.imm.cleaned",
                 shuffle = T, pt.size = 0.01, group.by="treatment", figure_plot = T)
ggsave(filename = file.path(save_path, "mainUMAP_by_treatment.pdf"), height = 8, width = 11)
ggsave(filename = file.path(save_path, "mainUMAP_by_treatment_small.pdf"), height = 4, width = 5.5)

DimPlot_scCustom(imm_obj_clean, colors_use = brewer.pal(6, "Set1"), reduction = "umap.imm.cleaned",
                 shuffle = T, pt.size = 0.01, group.by="sample", figure_plot = T)
ggsave(filename = file.path(save_path, "mainUMAP_by_sample.pdf"), height = 8, width = 11)
ggsave(filename = file.path(save_path, "mainUMAP_by_sample_small.pdf"), height = 4, width = 5.5)

# 5.2 heatmaps ====
# version with top 3 markers per celltype
immune.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

imm_top_degs = top3$gene

imm_pseudo <- AverageExpression(imm_obj_clean, assay = "SCT", return.seurat = T)

imm_pseudo_df = imm_pseudo@assays$SCT$scale.data
imm_pseudo_df = as.data.frame(imm_pseudo_df)
imm_pseudo_df = subset(imm_pseudo_df, rownames(imm_pseudo_df) %in% imm_top_degs)
imm_pseudo_df = na.omit(imm_pseudo_df)
imm_pseudo_df <- imm_pseudo_df[!grepl("\\.1$", row.names(imm_pseudo_df)), ]

#set up for heatmap
imm_pseudo_mat = as.matrix(imm_pseudo_df)

row_dend = dendsort::dendsort(hclust(dist(imm_pseudo_mat)))
col_dend = dendsort::dendsort(hclust(dist(t(imm_pseudo_mat))))

pdf(file=file.path(save_path, "heatmap_top3markers.pdf"), height = 3, width = 11)
ht <- ComplexHeatmap::Heatmap(t(imm_pseudo_mat), 
                              col = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "red")),
                              border_gp = gpar(col = "black"),
                              cluster_rows = col_dend, cluster_columns = row_dend,
                              show_row_dend = FALSE, show_column_dend = FALSE,
                              column_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                              row_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                              heatmap_legend_param = list( #specifying legend 
                                title = "", at = c(-4, 0, 4)),
                              row_names_centered = F, row_names_rot = 0,
                              column_names_rot = 45,
                              row_names_side = "left"
) 
draw(ht)
dev.off()

# version with canonical markers
imm_pseudo <- AverageExpression(imm_obj_clean, assay = "SCT", return.seurat = T)

imm_pseudo_df = imm_pseudo@assays$SCT$scale.data
imm_pseudo_df = as.data.frame(imm_pseudo_df)

imm_canonical <- c(
  "Adgre1", "Cd68",   # m0
  "Arg1", "Cd163", "Mrc1",  # m2-like
  "Cd86", "Cd80", "Nos2",        # m1-like
  "Cd83", "Itgax", "H2-Ab1", #pre-cDC
  "Ly6g", # neutrophil
  "Cd14", # monocyte
  "Clec9a", "Xcr1", # cDC1
  "Irf8", "Itgam", "Sirpa", # cDC2
  "Siglech", # pDC
  "Cd3e", "Klrk1", "Klrb1c", # nk
  "Tcf7", "Lef1", # tn
  "Tox", "Havcr2", # tex
  "Gzmb", "Cd8a", "Cd8b1", # ctl
  "Ctla4", "Foxp3", "Cd4", "Icos", # treg
  "Top2a", "Mki67", # tcyc
  "Ms4a1", "Cd79a", "Bank1", # b
  "Jchain", "Mzb1" # plasma
)

imm_pseudo_df = subset(imm_pseudo_df, rownames(imm_pseudo_df) %in% imm_canonical)
imm_pseudo_df = na.omit(imm_pseudo_df)
imm_pseudo_df <- imm_pseudo_df[!grepl("\\.1$", row.names(imm_pseudo_df)), ]

#set up for heatmap
imm_pseudo_mat = as.matrix(imm_pseudo_df)

row_dend = dendsort::dendsort(hclust(dist(imm_pseudo_mat)))
col_dend = dendsort::dendsort(hclust(dist(t(imm_pseudo_mat))))

pdf(file=file.path(save_path, "heatmap_canonicalmarkers.pdf"), height = 3, width = 11)
ht <- ComplexHeatmap::Heatmap(t(imm_pseudo_mat), 
                              col = circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "red")), # specifying color palette
                              border_gp = gpar(col = "black"),
                              cluster_rows = col_dend, cluster_columns = row_dend,
                              show_row_dend = FALSE, show_column_dend = FALSE,
                              column_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                              row_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                              heatmap_legend_param = list( #specifying legend 
                                title = "", at = c(-4, 0, 4)),
                              row_names_centered = F, row_names_rot = 0,
                              column_names_rot = 45,
                              row_names_side = "left"
) 
draw(ht)
dev.off()

# 5.3 t cell umap figures ====
Idents(object = t_obj_clean) <- "finalized_fine_anno"
t_obj_clean@meta.data$finalized_fine_anno <- factor(t_obj_clean@meta.data$finalized_fine_anno, levels = lymphoid_fine_order_rainbow[4:8])
DimPlot_scCustom(t_obj_clean, colors_use = lymphoid_fine_colors_rainbow[4:8], 
                 reduction = "umap.t.cleaned", shuffle = T, pt.size = 0.01, split.by = "treatment",
                 group.by="finalized_fine_anno", figure_plot = T)
ggsave(filename = file.path(save_path, "t_UMAP_by_celltype.pdf"), height = 4, width = 11)
ggsave(filename = file.path(save_path, "t_UMAP_by_celltype_small.pdf"), height = 2, width = 5.5)

# 5.4 bar plots ====
prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::group_by(treatment, finalized_fine_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = all_fine_order_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = all_fine_colors_rainbow)
ggsave(filename = file.path(save_path, "allcell_proportion_stacked.pdf"), height = 5, width = 3)


prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::filter(fine_anno %in% c("Treg", "CTL", "Tex", "Tn", "Tcyc")) %>%
  dplyr::group_by(treatment, finalized_fine_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = lymphoid_fine_order_rainbow[4:8])

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = lymphoid_fine_colors_rainbow[4:8])
ggsave(filename = file.path(save_path, "tcell_proportion_stacked.pdf"), height = 5, width = 3)

total_per_treatment <- prop_calcs %>%
  dplyr::group_by(treatment, sample) %>%
  dplyr::summarise(total_count = sum(count))

prop_calcs <- prop_calcs %>%
  dplyr::left_join(total_per_treatment, by = c("treatment", "sample")) %>%
  dplyr::mutate(cells_not_in_celltype = total_count - count)

count_glm <- glm(
  formula = cbind(count, cells_not_in_celltype) ~ finalized_fine_anno*treatment+sample,
  family = binomial(link = 'logit'),
  data = prop_calcs
)

emm_Count <- emmeans(count_glm, specs = revpairwise ~ treatment|finalized_fine_anno)
emm_Count$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

prop_calcs$prop = prop_calcs$count/prop_calcs$total_count

pw_pval_vec = c_results$p.value
pw_pval_vec_corrected = length(pw_pval_vec)
pw_pval_vec_corrected[pw_pval_vec > 0.05] <- "ns"
pw_pval_vec_corrected[pw_pval_vec <= 0.05] <- "*"
pw_pval_vec_corrected[pw_pval_vec <= 0.01] <- "**"
pw_pval_vec_corrected[pw_pval_vec < 0.001] <- "***"

p1<- ggplot(subset(prop_calcs, finalized_fine_anno == "Treg"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none", plot.title=element_text(hjust=0.5))+
  ggtitle("Treg")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  ylim(0, 0.75)+
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")), y_position = 0.45,
              annotations = pw_pval_vec_corrected[1], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)

p2<- ggplot(subset(prop_calcs, finalized_fine_anno == "CTL"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+ 
  theme_bw()+ theme(legend.position = "none", plot.title=element_text(hjust=0.5))+
  ggtitle("CTL")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  ylim(0, 0.75)+
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")), y_position = 0.65,
              annotations = pw_pval_vec_corrected[2], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)

p3<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tex"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none", plot.title=element_text(hjust=0.5))+
  ggtitle("Tex")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  ylim(0, 0.75)+
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")), y_position = 0.3,
              annotations = pw_pval_vec_corrected[3], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)

p4<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tn"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none", plot.title=element_text(hjust=0.5))+
  ggtitle("Tn")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  ylim(0, 0.75)+
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")), y_position = 0.5,
              annotations = pw_pval_vec_corrected[4], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)

p5<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tcyc"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none", plot.title=element_text(hjust=0.5))+
  ggtitle("Tcyc")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  ylim(0, 0.75)+
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")), y_position = 0.25,
              annotations = pw_pval_vec_corrected[5], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)

cowplot::plot_grid(p1, 
                   p2 + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()), 
                   p3 + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()),
                   p4 + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()), 
                   p5 + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank()),
                   nrow = 1, align = "v")
ggsave(filename = file.path(save_path, "tcell_proportion_barplots.pdf"), height = 3, width = 10)

# 5.5 degs and gsea figures ====
tcell_degs$sig = rep("NS")
tcell_degs$sig[tcell_degs$p_val_adj<0.05 & tcell_degs$avg_log2FC>1] <- "*up"
tcell_degs$sig[tcell_degs$p_val_adj<0.05 & tcell_degs$avg_log2FC<(-1)] <- "*down"
tcell_degs$gene_name = tcell_degs$X

genes_to_include = c("Gzmb", "Ccl5", "Prf1", "Ifng", "Nkg7", "Cd8a",
                     "Tcf7", "Foxp1", "Sell", "Il7r", "Foxo1", "Lef1")

# volcano
ggplot(tcell_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = sig), size = 1.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Transformed threshold line
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +  # Transformed threshold line
  scale_y_continuous(
    name = "p-value (adjusted)",
    # Apply both a reverse log10 transformation and negative transformation
    labels = function(x) format(10^(-x), scientific = TRUE)) +
  labs(x = "log2(PD1+Corin/PD1)") +
  scale_color_manual(values = c("grey30", "#c778ff", "#b2b2b2"))+
  guides(color=guide_legend(title="Significance"))+
  theme_bw() +
  geom_label_repel(data = subset(tcell_degs, tcell_degs$gene_name %in% genes_to_include),
                   aes(label = gene_name), min.segment.length = 0,
                   max.overlaps = Inf, size = 3, segment.size =0.3, force = 20) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")
ggsave(filename = file.path(save_path, "tcell_degs_volcano.pdf"), height = 5, width = 5)

# violins
Idents(object = t_obj_clean) <- "treatment"
Stacked_VlnPlot(t_obj_clean, features = genes_to_include[1:6], 
                colors_use = c("grey30", "#c778ff"))
ggsave(filename = file.path(save_path, "tcell_degs_up_violins.pdf"), height = 8, width = 3)

Stacked_VlnPlot(t_obj_clean, features = genes_to_include[7:12], 
                colors_use = c("grey30", "#c778ff"))
ggsave(filename = file.path(save_path, "tcell_degs_down_violins.pdf"), height = 8, width = 3)

# bubble plot
tohightlight = c(
  "GOLDRATH_ANTIGEN_RESPONSE", 
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "REACTOME_SEPARATION_OF_SISTER_CHROMATIDS",
  "GOMF_CHEMOATTRACTANT_ACTIVITY",
  "GOMF_CHEMOKINE_ACTIVITY",
  "GOMF_CYTOKINE_ACTIVITY",
  "GOMF_CHEMOKINE_RECEPTOR_BINDING",
  "GOBP_LEUKOCYTE_CHEMOTAXIS",
  "GOBP_LYMPHOCYTE_CHEMOTAXIS",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "GOBP_APOPTOTIC_CELL_CLEARANCE",
  "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL",
  "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_CELL_KILLING",
  "WP_INFLAMMATORY_RESPONSE_PATHWAY",
  "GOBP_LEUKOCYTE_MIGRATION",
  "GOBP_INFLAMMATORY_RESPONSE"
)

colos = setNames(c("grey30", "#c778ff"), c("Down-regulated", "Up-regulated"))

cherry = subset(fgsea_results$all_t.cells, pathway %in% tohightlight)

cherry$Enrichment = ifelse(cherry$NES > 0, "Up-regulated", "Down-regulated")

ggplot(cherry, aes(reorder(pathway, NES), NES)) +
  geom_point( aes(fill = Enrichment, size = padj, col = Enrichment), shape=21) +
  scale_fill_manual(values = colos, guide = "none" ) +
  scale_color_manual(values = colos, guide = "none" ) +
  scale_size_continuous(range = c(2,8), trans = c("log10", "reverse")) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score", size = "Adjusted P Value") +
  theme_bw() +
  theme(axis.line.y.left =element_line(color="white"),
        axis.ticks.y=element_blank()) + 
  guides(size = guide_legend(override.aes = list(color = "black", fill = "white")))
ggsave(filename = file.path(save_path, "tcell_gsea_bubbleplot.pdf"), height = 5, width = 8)

# enrichment plots
t_degs <- setNames(tcell_degs$avg_log2FC, tcell_degs$X)
t_degs <- sort(t_degs, decreasing = TRUE)

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

plotEnrichment_customcolor(all_msigdb[["GOLDRATH_ANTIGEN_RESPONSE"]],
                           t_degs, ticksSize = 0.1, custom_color = '#c778ff') +
  ggtitle("GOLDRATH_ANTIGEN_RESPONSE") +
  annotate("text", x = 9000, y = 0.5, hjust = 1,
           label = paste0("NES = ", round(subset(fgsea_results$all_t.cells, pathway == "GOLDRATH_ANTIGEN_RESPONSE")$NES, 3))) +
  annotate("text", x = 9000, y = 0.45, hjust = 1,
           label = paste0("Padj = ", format(signif(subset(fgsea_results$all_t.cells, pathway == "GOLDRATH_ANTIGEN_RESPONSE")$padj, 3), scientific = TRUE))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(save_path, "goldrathantigenresponse_enrichmentplot.pdf"), height = 2.5, width = 4.5)

plotEnrichment_customcolor(all_msigdb[["GOMF_CYTOKINE_ACTIVITY"]],
                           t_degs, ticksSize = 0.1, custom_color = '#c778ff') +
  ggtitle("GOMF_CYTOKINE_ACTIVITY") +
  annotate("text", x = 9000, y = 0.5, hjust = 1,
           label = paste0("NES = ", round(subset(fgsea_results$all_t.cells, pathway == "GOMF_CYTOKINE_ACTIVITY")$NES, 3))) +
  annotate("text", x = 9000, y = 0.45, hjust = 1,
           label = paste0("Padj = ", format(signif(subset(fgsea_results$all_t.cells, pathway == "GOMF_CYTOKINE_ACTIVITY")$padj, 3), scientific = TRUE))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(save_path, "cytokineactivity_enrichmentplot.pdf"), height = 2.5, width = 4.5)

plotEnrichment_customcolor(all_msigdb[["GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL"]],
                           t_degs, ticksSize = 0.1, custom_color = '#c778ff') +
  ggtitle("GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL") +
  annotate("text", x = 9000, y = 0.5, hjust = 1,
           label = paste0("NES = ", round(subset(fgsea_results$all_t.cells, pathway == "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL")$NES, 3))) +
  annotate("text", x = 9000, y = 0.45, hjust = 1,
           label = paste0("Padj = ", format(signif(subset(fgsea_results$all_t.cells, pathway == "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL")$padj, 3), scientific = TRUE))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(save_path, "immuneresponsetotumorcell_enrichmentplot.pdf"), height = 2.5, width = 4.5)

plotEnrichment_customcolor(all_msigdb[["GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE"]],
                           t_degs, ticksSize = 0.1, custom_color = '#c778ff') +
  ggtitle("GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE") +
  annotate("text", x = 9000, y = 0.5, hjust = 1,
           label = paste0("NES = ", round(subset(fgsea_results$all_t.cells, pathway == "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$NES, 3))) +
  annotate("text", x = 9000, y = 0.45, hjust = 1,
           label = paste0("Padj = ", format(signif(subset(fgsea_results$all_t.cells, pathway == "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$padj, 3), scientific = TRUE))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(save_path, "leukocytemigration_enrichmentplot.pdf"), height = 2.5, width = 4.5)
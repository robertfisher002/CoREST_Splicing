# 4.0 load object ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
options(Seurat.object.assay.version = "v5")
basefile_path <- "/path/to/files" # please download files and update this line to point to your directory
imm_obj_clean <- readRDS(file=file.path(save_path, "corin_pd1_imm_obj_cleaned.rds"))

# please export gmt files from msigdb: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
all_msigdb = read.gmt("path_to/msigdb.v2024.1.Mm.symbols.gmt")
hallmarks = read.gmt("path_to/mh.all.v2024.1.Mm.symbols.gmt")
gobp = read.gmt("path_to/m5.go.bp.v2024.1.Mm.symbols.gmt")
cp = read.gmt("path_to/m2.cp.v2024.1.Mm.symbols.gmt")

msigdb_of_interest = list(all = all_msigdb, hallmarks = hallmarks, 
                          gobp = gobp, cp = cp)
# 4.1 heatmap of markers ====
p1 <- DimPlot_scCustom(imm_obj_clean, reduction = "umap.imm.clean", group.by = "broad_anno", colors_use = all_broad_colors_alpha)
p2 <- DimPlot_scCustom(imm_obj_clean, reduction = "umap.imm.clean", group.by = "finalized_fine_anno", colors_use = all_fine_colors_alpha)
plot_grid(p1, p2)

Idents(object = imm_obj_clean) <- "finalized_fine_anno"
imm_obj_clean <- PrepSCTFindMarkers(imm_obj_clean)
immune.markers <- FindAllMarkers(imm_obj_clean, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

write.csv(immune.markers, 
          file=file.path(save_path, "all_immune_markers_MAST.csv"), 
          row.names=TRUE)

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

extra_rdbu <- colorRampPalette(rev(brewer.pal(11,"RdBu")), bias = 2.4)

# plotting
# super helpful: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
ComplexHeatmap::Heatmap(t(imm_pseudo_mat), 
                        col = extra_rdbu(100), # specifying color palette
                        rect_gp = gpar(col = "white", lwd = 1), # adding borders
                        cluster_rows = col_dend, cluster_columns = row_dend,
                        show_row_dend = FALSE, show_column_dend = FALSE,
                        column_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                        row_names_gp = gpar(fontsize = 8), # setting fontsize of genes
                        heatmap_legend_param = list( #specifying legend 
                          title = "", at = c(-2, 0, 4)),
                        row_names_centered = F, row_names_rot = 0,
                        column_names_rot = 45,
                        row_names_side = "left"
) 

# 4.2 all immune broad proportion analysis ====
prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::group_by(treatment, broad_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$broad_anno = factor(prop_calcs$broad_anno, levels = all_broad_order_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=broad_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = all_broad_colors_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=broad_anno), stat = "identity")+theme_pubr(legend="right")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Count", fill = "")+
  scale_fill_manual(values = all_broad_colors_rainbow)

# finding significance
total_per_treatment <- prop_calcs %>%
  dplyr::group_by(treatment, sample) %>%
  dplyr::summarise(total_count = sum(count))

prop_calcs <- prop_calcs %>%
  dplyr::left_join(total_per_treatment, by = c("treatment", "sample")) %>%
  dplyr::mutate(cells_not_in_celltype = total_count - count)

count_glm <- glm(
  formula = cbind(count, cells_not_in_celltype) ~ broad_anno*treatment+sample,
  family = binomial(link = 'logit'),
  data = prop_calcs
)

emm_Count <- emmeans(count_glm, specs = revpairwise ~ treatment|broad_anno)
emm_Count$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results
c_results

# broad_anno = Monocyte:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.52437 0.030017 Inf   0.46872   0.58663    1 -11.277  <.0001
# 
# broad_anno = Macrophage:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    2.98231 0.088886 Inf   2.81309   3.16171    1  36.662  <.0001
# 
# broad_anno = Neutrophil:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   46.13203 7.304963 Inf  33.82327  62.92011    1  24.197  <.0001
# 
# broad_anno = DC:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.94795 0.087188 Inf   1.78434   2.12655    1  14.897  <.0001
# 
# broad_anno = T/NK:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.53205 0.018516 Inf   0.49697   0.56961    1 -18.132  <.0001
# 
# broad_anno = B:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.01796 0.001679 Inf   0.01495   0.02157    1 -42.991  <.0001

# 4.3 all immune fine proportion analysis ====
prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::group_by(treatment, finalized_fine_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = all_fine_order_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = all_fine_colors_rainbow)

# finding significance
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
c_results

# finalized_fine_anno = Monocyte:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.52437 0.030007 Inf   0.46873   0.58661    1 -11.281  <.0001
# 
# finalized_fine_anno = M1-like:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.82839 0.033474 Inf   0.76531   0.89666    1  -4.659  <.0001
# 
# finalized_fine_anno = M0:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.95413 0.110280 Inf   1.74951   2.18268    1  11.871  <.0001
# 
# finalized_fine_anno = Neutrophil:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   46.13203 7.304649 Inf  33.82372  62.91927    1  24.198  <.0001
# 
# finalized_fine_anno = pDC:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.60600 0.152123 Inf   1.33388   1.93363    1   5.001  <.0001
# 
# finalized_fine_anno = pre-cDC:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    2.09758 0.128246 Inf   1.86070   2.36462    1  12.116  <.0001
# 
# finalized_fine_anno = cDC2:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    2.02279 0.282962 Inf   1.53773   2.66087    1   5.036  <.0001
# 
# finalized_fine_anno = cDC1:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.47931 0.143200 Inf   1.22366   1.78837    1   4.045  0.0001
# 
# finalized_fine_anno = M2-like:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    6.62363 0.322852 Inf   6.02013   7.28762    1  38.788  <.0001
# 
# finalized_fine_anno = NK:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.02489 0.069060 Inf   0.89809   1.16959    1   0.365  0.7152
# 
# finalized_fine_anno = B:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.01914 0.001850 Inf   0.01584   0.02313    1 -40.941  <.0001
# 
# finalized_fine_anno = Plasma:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.02063 0.007402 Inf   0.01021   0.04168    1 -10.817  <.0001
# 
# finalized_fine_anno = Treg:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.65009 0.059574 Inf   0.54321   0.77800    1  -4.699  <.0001
# 
# finalized_fine_anno = CTL:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.82366 0.125128 Inf   1.59419   2.08616    1   8.757  <.0001
# 
# finalized_fine_anno = Tex:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.07513 0.012845 Inf   0.05374   0.10504    1 -15.141  <.0001
# 
# finalized_fine_anno = Tn:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    0.03795 0.005835 Inf   0.02808   0.05130    1 -21.277  <.0001
# 
# finalized_fine_anno = Tcyc:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1    1.20788 0.124478 Inf   0.98696   1.47823    1   1.833  0.0669

# 4.4 myeloid proportion analysis ====
prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::filter(coarse_anno =="myeloid") %>%
  dplyr::group_by(treatment, finalized_fine_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = myeloid_fine_order_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = myeloid_fine_colors_rainbow)

# finding significance
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
c_results

# finalized_fine_anno = Monocyte:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.223869 0.013262 Inf  0.199328  0.251431    1 -25.265  <.0001
# 
# finalized_fine_anno = M1-like:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.316070 0.013886 Inf  0.289992  0.344493    1 -26.216  <.0001
# 
# finalized_fine_anno = M0:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.892934 0.052047 Inf  0.796535  1.001001    1  -1.943  0.0520
# 
# finalized_fine_anno = Neutrophil:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1  23.036711 3.661446 Inf 16.870613 31.456479    1  19.738  <.0001
# 
# finalized_fine_anno = pDC:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.762077 0.072948 Inf  0.631712  0.919344    1  -2.838  0.0045
# 
# finalized_fine_anno = pre-cDC:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.970184 0.060934 Inf  0.857814  1.097274    1  -0.482  0.6298
# 
# finalized_fine_anno = cDC2:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.972089 0.136623 Inf  0.738029  1.280380    1  -0.201  0.8404
# 
# finalized_fine_anno = cDC1:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.701943 0.068636 Inf  0.579524  0.850223    1  -3.619  0.0003
# 
# finalized_fine_anno = M2-like:
#   contrast          odds.ratio       SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   3.161185 0.161348 Inf  2.860253  3.493779    1  22.550  <.0001

# 4.5 lymphoid proportion analysis ====
prop_calcs = imm_obj_clean@meta.data %>%
  dplyr::filter(coarse_anno =="lymphoid") %>%
  dplyr::group_by(treatment, finalized_fine_anno, sample) %>%
  dplyr::summarise(count = n())
prop_calcs$finalized_fine_anno = factor(prop_calcs$finalized_fine_anno, levels = lymphoid_fine_order_rainbow)

ggplot(prop_calcs, aes(x=treatment, y=count))+
  geom_bar(aes(fill=finalized_fine_anno), stat = "identity", position = "fill")+theme_pubr(legend="right")+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = lymphoid_fine_colors_rainbow)

# finding significance
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
c_results

# finalized_fine_anno = NK:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   5.035887 0.3776194 Inf  4.347585  5.833160    1  21.559  <.0001
# 
# finalized_fine_anno = B:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.051967 0.0053018 Inf  0.042549  0.063470    1 -28.985  <.0001
# 
# finalized_fine_anno = Plasma:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.077655 0.0279556 Inf  0.038348  0.157253    1  -7.099  <.0001
# 
# finalized_fine_anno = Treg:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   2.713155 0.2623819 Inf  2.244694  3.279382    1  10.321  <.0001
# 
# finalized_fine_anno = CTL:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1  10.328742 0.7934991 Inf  8.884938 12.007164    1  30.393  <.0001
# 
# finalized_fine_anno = Tex:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.284596 0.0494000 Inf  0.202524  0.399926    1  -7.240  <.0001
# 
# finalized_fine_anno = Tn:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.138709 0.0217363 Inf  0.102027  0.188578    1 -12.606  <.0001
# 
# finalized_fine_anno = Tcyc:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   5.135532 0.5526212 Inf  4.159009  6.341340    1  15.205  <.0001

# 4.6 t only proportion analysis ====
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

# finding significance
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
c_results

# finalized_fine_anno = Treg:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   1.343355 0.1402041 Inf  1.094844  1.648275    1   2.828  0.0047
# 
# finalized_fine_anno = CTL:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   6.756374 0.6071497 Inf  5.665288  8.057595    1  21.260  <.0001
# 
# finalized_fine_anno = Tex:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.126833 0.0225552 Inf  0.089508  0.179724    1 -11.611  <.0001
# 
# finalized_fine_anno = Tn:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   0.047976 0.0077803 Inf  0.034913  0.065927    1 -18.727  <.0001
# 
# finalized_fine_anno = Tcyc:
#   contrast          odds.ratio        SE  df asymp.LCL asymp.UCL null z.ratio p.value
# (PD1+Corin) / PD1   2.673369 0.3059679 Inf  2.136183  3.345639    1   8.592  <.0001

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
  theme_bw()+ theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # to ensure p values are not partially hidden
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")),
              annotations = pw_pval_vec_corrected[1], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)
p2<- ggplot(subset(prop_calcs, finalized_fine_anno == "CTL"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # to ensure p values are not partially hidden
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")),
              annotations = pw_pval_vec_corrected[2], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)
p3<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tex"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # to ensure p values are not partially hidden
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")),
              annotations = pw_pval_vec_corrected[3], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)
p4<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tn"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # to ensure p values are not partially hidden
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")),
              annotations = pw_pval_vec_corrected[4], 
              step_increase = 0.1, textsize = 8, tip_length = 0.01)
p5<- ggplot(subset(prop_calcs, finalized_fine_anno == "Tcyc"), aes(x=treatment, y=prop))+
  geom_bar(aes(fill=treatment), stat='summary', color = "black", size = 0.2) + 
  geom_jitter(color="black", width = 0.2, size = 2)+
  geom_errorbar(stat='summary', width=.2)+
  theme_bw()+ theme(legend.position = "none")+
  labs(x="", y = "Proportion", fill = "")+
  scale_fill_manual(values = c("grey30", "#c778ff"))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # to ensure p values are not partially hidden
  geom_signif(comparisons = list(c("PD1", "PD1+Corin")),
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
                   nrow = 1)

# 4.7 t only differential expression and GSEA ====
Idents(object = t_obj_cleaned) <- "treatment"
tcell_pd1_v_pd1corin <- FindMarkers(t_obj_cleaned, ident.1 = "PD1+Corin", ident.2 = "PD1", test.use="MAST")

write.csv(tcell_pd1_v_pd1corin, file=file.path(save_path, "tcell_pd1corin_v_pd1_corrected.csv"))

# running through gsea
tcell_degs = read.csv(file=file.path(save_path, "tcell_pd1corin_v_pd1_corrected.csv"))
markers_list = list(t.cells = tcell_degs)

# running setup functions
combined_volcano_plot(markers_list)
imm_for_gsea = create_gsea_vectors(markers_list)
fgsea_results = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                  named_degs_list = imm_for_gsea) # this takes a few minutes
saveRDS(fgsea_results, file=file.path(save_path, "fgsea_results_tcells_pd1corin_v_pd1.rds"))

fgsea_results = readRDS(file=file.path(save_path, "fgsea_results_tcells_pd1corin_v_pd1.rds"))

# making summary/nice gsea plot
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

colos = setNames(c("dodgerblue2", "orange"), c("Up-regulated", "Down-regulated"))

cherry = subset(fgsea_results$all_t.cells, pathway %in% tohightlight)

cherry$Enrichment = ifelse(cherry$NES > 0, "Up-regulated", "Down-regulated")

ggplot(cherry, aes(reorder(pathway, NES), NES)) +
  geom_vline(xintercept = c(1:18), alpha = 0.1, size = 0.5) +
  geom_hline(yintercept = seq(from = 0, to = 3, by = 1), alpha = 0.2) +
  geom_point( aes(fill = Enrichment, size = padj, col = Enrichment), shape=21) +
  scale_fill_manual(values = colos, guide = "none" ) +
  scale_color_manual(values = colos, guide = "none" ) +
  scale_size_continuous(range = c(2,8), trans = c("log10", "reverse")) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  labs(x="", y="Normalized Enrichment Score", size = "Adjusted P Value") +
  theme_pubr(legend = "right") +
  theme(axis.line.y.left =element_line(color="white"),
        axis.ticks.y=element_blank()) + 
  guides(size = guide_legend(override.aes = list(color = "black", fill = "white")))

# 4.8 t states differential expression and GSEA ====
# treg
treg_obj <- subset(x = imm_obj_clean, subset = finalized_fine_anno %in% c("Treg"))
treg_obj <- SCTransform(treg_obj)
treg_obj <- RunPCA(treg_obj)
treg_obj <- FindNeighbors(treg_obj, dims = 1:20)
Idents(object = treg_obj) <- "treatment"
treg_pd1_v_pd1corin <- FindMarkers(treg_obj, ident.1 = "PD1+Corin", ident.2 = "PD1", test.use="MAST")
write.csv(treg_pd1_v_pd1corin, file=file.path(save_path, "treg_pd1corin_v_pd1_corrected.csv"))

# tex
tex_obj <- subset(x = imm_obj_clean, subset = finalized_fine_anno %in% c("Tex"))
tex_obj <- SCTransform(tex_obj)
tex_obj <- RunPCA(tex_obj)
tex_obj <- FindNeighbors(tex_obj, dims = 1:20)
Idents(object = tex_obj) <- "treatment"
tex_pd1_v_pd1corin <- FindMarkers(tex_obj, ident.1 = "PD1+Corin", ident.2 = "PD1", test.use="MAST")
write.csv(tex_pd1_v_pd1corin, file=file.path(save_path, "tex_pd1corin_v_pd1_corrected.csv"))

# ctl
ctl_obj <- subset(x = imm_obj_clean, subset = finalized_fine_anno %in% c("CTL"))
ctl_obj <- SCTransform(ctl_obj)
ctl_obj <- RunPCA(ctl_obj)
ctl_obj <- FindNeighbors(ctl_obj, dims = 1:20)
Idents(object = ctl_obj) <- "treatment"
ctl_pd1_v_pd1corin <- FindMarkers(ctl_obj, ident.1 = "PD1+Corin", ident.2 = "PD1", test.use="MAST")
write.csv(ctl_pd1_v_pd1corin, file=file.path(save_path, "ctl_pd1corin_v_pd1_corrected.csv"))

# running through gsea
treg_degs = read.csv(file=file.path(save_path, "treg_pd1corin_v_pd1_corrected.csv"))
tex_degs = read.csv(file=file.path(save_path, "tex_pd1corin_v_pd1_corrected.csv"))
ctl_degs = read.csv(file=file.path(save_path, "ctl_pd1corin_v_pd1_corrected.csv"))
markers_list = list(treg = treg_degs, tex = tex_degs, ctl = ctl_degs)

# running setup functions
combined_volcano_plot(markers_list)
imm_for_gsea = create_gsea_vectors(markers_list)
fgsea_results = run_through_fgsea(msigdb_list = msigdb_of_interest, 
                                  named_degs_list = imm_for_gsea) # this takes a few minutes
saveRDS(fgsea_results, file=file.path(save_path, "fgsea_results_tcellstates_pd1corin_v_pd1.rds"))

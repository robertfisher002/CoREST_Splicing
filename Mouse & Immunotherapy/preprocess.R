# 1.0 load packages and initial data ====
source("../scripts/utils.R")
save_path <- "path/to/save/dir/" # please update file path to save in your desired location
options(Seurat.object.assay.version = "v5")
basefile_path <- "/path/to/files" # please download files and update this line to point to your directory

data_1A4 <- Read10X(data.dir = file.path(basefile_path, "1_A4"))
data_1B3 <- Read10X(data.dir = file.path(basefile_path, "1_B3"))
data_1B5 <- Read10X(data.dir = file.path(basefile_path, "1_B5"))
data_4A3 <- Read10X(data.dir = file.path(basefile_path, "4_A3"))
data_4A4 <- Read10X(data.dir = file.path(basefile_path, "4_A4"))
data_4B3 <- Read10X(data.dir = file.path(basefile_path, "4_B3"))

# initialize the Seurat object with the raw (non-normalized data).
d1A4 <- CreateSeuratObject(counts = data_1A4, project = "pdo", min.cells = 3, min.features = 200) # 6602 cells
d1B3 <- CreateSeuratObject(counts = data_1B3, project = "pdo", min.cells = 3, min.features = 200) # 5688 cells
d1B5 <- CreateSeuratObject(counts = data_1B5, project = "pdo", min.cells = 3, min.features = 200) # 6818 cells
d4A3 <- CreateSeuratObject(counts = data_4A3, project = "pdo", min.cells = 3, min.features = 200) # 12812 cells
d4A4 <- CreateSeuratObject(counts = data_4A4, project = "pdo", min.cells = 3, min.features = 200) # 8600 cells
d4B3 <- CreateSeuratObject(counts = data_4B3, project = "pdo", min.cells = 3, min.features = 200) # 8704 cells

# 1.1 doublet detection ====
# must be performed on processed single lanes - creating function so this isn't hundreds of lines long
# used 10X 5' v3 kit - https://www.10xgenomics.com/support/single-cell-immune-profiling/documentation/steps/library-prep/chromium-gem-x-single-cell-5-v3-gene-expression-user-guide
# define a list of objects with their parameters
raw_objs <- list(
  list(object = d1A4, expected_rate = 0.028), # expected doublet rate when recovering ~ 6600 cells is ~2.8%
  list(object = d1B3, expected_rate = 0.024), # expected doublet rate when recovering ~5688 cells is ~2.4%
  list(object = d1B5, expected_rate = 0.028), # expected doublet rate when recovering ~6818 cells is ~2.8%
  list(object = d4A3, expected_rate = 0.052), # expected doublet rate when recovering ~12,812 cells is ~5.2%
  list(object = d4A4, expected_rate = 0.036), # expected doublet rate when recovering ~8600 cells is ~3.6%
  list(object = d4B3, expected_rate = 0.036) # expected doublet rate when recovering ~8704 cells is ~3.6%
)

names(raw_objs) = c("d1A4", "d1B3", "d1B5", "d4A3", "d4A4", "d4B3")

# process, calculate pK, and detect doublets
doublet_detected_objs <- lapply(seq_along(raw_objs), function(x) {
  obj <- raw_objs[[x]]$object
  expected_rate <- raw_objs[[x]]$expected_rate
  cat(paste("Processing object", x, "of", length(raw_objs), "\n"))
  processed_obj <- process_and_doublet_finder(obj, expected_rate)
  cat(paste("Finished processing object", x, "with optimal pK value\n\n"))
  return(processed_obj)
})

names(doublet_detected_objs) = c("d1A4", "d1B3", "d1B5", "d4A3", "d4A4", "d4B3") # double checked the order and this is correct

# checking results
doublet_dimplot_check(doublet_detected_objs)
for (obj in doublet_detected_objs) {
  print_metadata_table(obj, "DF.classifications_")
}

# adding doublets as column in real data
raw_doublet_detected_objs <- update_doublets(raw = raw_objs, doublet = doublet_detected_objs)

# 1.2 QC visualizations ====
# iterate over each object in the list and generate plots
for (i in seq_along(raw_doublet_detected_objs)) {
  raw_doublet_detected_objs[[i]] <- qc_metrics(raw_doublet_detected_objs[[i]]$object)
}

for (i in seq_along(raw_doublet_detected_objs)) {
  qc_plots(raw_doublet_detected_objs[[i]], names(raw_doublet_detected_objs)[i])
}

# 1.3 low quality cell filtering ====
# apply the function to each object in the list
subset_objs <- mapply(qc_filtering, raw_doublet_detected_objs, names(raw_doublet_detected_objs), SIMPLIFY = FALSE)

# Number of cells in d1A4 : 4347 
# Number of cells in d1B3 : 3684 
# Number of cells in d1B5 : 4316 
# Number of cells in d4A3 : 7932 
# Number of cells in d4A4 : 4028 
# Number of cells in d4B3 : 4775 

#### finalize objects ====
# 2.1 merge into one object ====
# merge together
data_obj <- merge(subset_objs$d1A4, y = c(subset_objs$d1B3, subset_objs$d1B5, subset_objs$d4A3, subset_objs$d4A4, subset_objs$d4B3),
                  add.cell.ids = c("d1A4", "d1B3", "d1B5", "d4A3", "d4A4", "d4B3"),
                  project = "pdo")
data_obj #29,082 cells (already undergone filtering and doublet removal)

# 2.2 add group labels ====
data_obj@meta.data$barcode = rownames(data_obj@meta.data)

data_obj@meta.data$sample = rownames(data_obj@meta.data)
data_obj@meta.data$sample = gsub("_.*", "", data_obj@meta.data$sample)
data_obj@meta.data$sample = gsub("d", "", data_obj@meta.data$sample)

data_obj@meta.data$treatment[data_obj@meta.data$sample %in% c("1A4", "1B3", "1B5")] <- "PD1"
data_obj@meta.data$treatment[data_obj@meta.data$sample %in% c("4A3", "4A4", "4B3")] <- "PD1+Corin"

table(data_obj@meta.data$treatment)
# PD1 PD1+Corin 
# 12347     16735 

# 2.3 check QC now that everything is together ====
p1 <- scCustomize::QC_Plots_Genes(seurat_object = data_obj, low_cutoff = 100, high_cutoff = 8500)
p2 <- scCustomize::QC_Plots_UMIs(seurat_object = data_obj, low_cutoff = 500, high_cutoff = 50000)
p3 <- scCustomize::QC_Plots_Mito(seurat_object = data_obj, high_cutoff = 20)
p4 <- scCustomize::QC_Plots_Complexity(seurat_object = data_obj, high_cutoff = 0.75)

wrap_plots(p1, p2, p3, p4, ncol = 4) #looks great

data_obj@meta.data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    min_value = min(log10GenesPerUMI),
    max_value = max(log10GenesPerUMI),
    median_value = median(log10GenesPerUMI)
  )

# 3.5 save base object ====
saveRDS(data_obj, file=file.path(save_path, "corin_pd1_obj.rds"))
# ---------------------------
# DoubletFinder and merge SeuratObject
# ---------------------------
library(harmony)
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(stringi)
library(plyr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(viridis)
library(reshape2)
library(scales)
library(rlang)
library(dendextend)
library(DoubletFinder)
library(future)
library(ggsci)
library(ggpubr)
library(Matrix)
library(irlba)
library(Nebulosa)

##### --------- part1.Doubletfinder ---------#####
# Parameter Set
sample_list <- c("Stress_1","Stress_2","Stress_3","Control_1","Control_2","Control_3")
res.use <- 0.5
dim.use <- 30
seed.use <- 88
project.name <- "mdd.liver"

# path set
data_dir_base <- "/public/home/chidm/Workspace/06_mdd.liver/singlesample"
out_dir <- "/public/home/chidm/Workspace/06_mdd.liver/paper/doublet-new"
setwd(out_dir)
qc_plot_dir <- file.path(out_dir, "qc_plots")
dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)

# doublet rate set（reference: doubletfinder）
get_doublet_rate <- function(cell_number) {
  rate_table <- data.frame(
    rate = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )
  idx <- which.min(abs(rate_table$recovered - cell_number))
  return(rate_table$rate[idx])
}

sample_list_after_doublet <- list()
qc_summary <- data.frame()

for (i in sample_list) {
  count.data <- Read10X(data.dir = file.path(data_dir_base, i, "2.Basic_analysis/2.2.filtered_feature_bc_matrix"))
  obj <- CreateSeuratObject(counts = count.data,min.cells = 3,min.features = 200, project = project.name)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  # QC plot：Vlnplot
  png(filename = file.path(qc_plot_dir, paste0(i, "_before_filter_violin.png")), width = 10, height = 4, units = "in", res = 400)
  p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) + plot_layout(guides = "collect")
  print(p1)
  dev.off()
  
  # record  cell number before filtration 
  cells_before <- ncol(obj)
  
  # filter low-quality cells
  obj <- subset(obj, subset =
                  nFeature_RNA >= 200 & nFeature_RNA <= 10000 & #不过滤nCount了
                  percent.mt <= 25)
  
  cells_after_qc <- ncol(obj)
  
  # standard pipline
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:dim.use)
  obj <- FindClusters(obj, resolution = res.use)
  obj <- RunUMAP(obj, dims = 1:dim.use)
  
  # DoubletFinder
  sweep.res.list <- paramSweep(obj, PCs = 1:dim.use)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  pdf(paste0(i, "_find_pK_plot.pdf"))
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  
  pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"] %>% as.character() %>% as.numeric()
  cell_count <- ncol(obj)
  rate <- get_doublet_rate(cell_count)
  nExp_poi <- round(rate * cell_count)
  homotypic.prop <- modelHomotypic(obj@active.ident)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  obj <- doubletFinder(obj, PCs = 1:dim.use, pN = 0.25, pK = pK, nExp = nExp_poi)
  obj <- doubletFinder(obj, PCs = 1:dim.use, pN = 0.25, pK = pK, nExp = nExp_poi.adj,
                       reuse.pANN = grep("pANN", colnames(obj@meta.data), value = TRUE))
  
  hi_lo <- obj@meta.data[, grep("^DF\\.classifications", names(obj@meta.data))] == "Singlet"
  hi_lo <- as.data.frame(hi_lo) + 0
  obj$DF_hi.lo <- ifelse(rowSums(hi_lo) == 2, "Singlet",
                         ifelse(rowSums(hi_lo) == 1, "Doublet_lo", "Doublet_hi"))
  
  # filter Doublet_hi，keep Singlet and Doublet_lo
  obj <- subset(obj, subset = DF_hi.lo %in% c("Singlet", "Doublet_lo"))
  cells_after_doublet <- ncol(obj)
  
  # Renormalization & clustering
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, vars.to.regress = "nCount_RNA")
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:dim.use)
  obj <- FindClusters(obj, resolution = res.use)
  obj <- RunUMAP(obj, dims = 1:dim.use)
  
  # save UMAP plot
  png(paste0(out_dir, "/", i, "_UMAP_after_doublet_removal.png"), width = 6, height = 5, units = "in", res = 400)
  p <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.3)
  print(p)
  dev.off()
  
  # save object after doublet
  sample_list_after_doublet[[i]] <- obj
  
  # QC table
  qc_summary <- rbind(qc_summary, data.frame(
    sample = i,
    cells_before = cells_before,
    cells_after_qc = cells_after_qc,
    cells_after_doublet_removal = cells_after_doublet,
    estimated_doublet_rate = rate
  ))
}
write.csv(qc_summary, file = file.path(out_dir, "QC_summary_all_samples.csv"),  row.names = FALSE)

##### -------- merge Seurat object -------- ######
sample_order <- c("Stress_1", "Stress_2", "Stress_3", "Control_1", "Control_2", "Control_3")
sample_list_doubletRM <- list()
## ---- rename barcode and save as list ----
for (i in 1:length(sample_order)) {
  sample <- sample_order[i]
  suffix <- i  
  obj <- sample_list_after_doublet[[sample]]
  new_colnames <- paste0(gsub("-1$", "", colnames(obj)), "-", suffix)
  colnames(obj) <- new_colnames
  sample_list_doubletRM[[sample]] <- obj
}

### merge Seuratobject
merged_seurat <- merge(
  x = sample_list_doubletRM[[1]],
  y = sample_list_doubletRM[2:length(sample_list_doubletRM)],
  add.cell.ids = names(sample_list_doubletRM),
  project = "mdd.liver"
)
### add metadata of sample
merged_seurat$sample <-paste0(str_split(colnames(merged_seurat),"_",simplify=T)[,1],"_",
                              str_split(colnames(merged_seurat),"_",simplify=T)[,2])

# add group information
merged_seurat$group <- ifelse(grepl("^Stress", merged_seurat$sample), "CRS", "Control")

save(merged_seurat, file = file.path(out_dir, "mdd_liver_doublet_removed_merged.RData"))#24883 features across 52848 samples

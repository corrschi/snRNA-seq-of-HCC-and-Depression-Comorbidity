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
#参数设置
sample_list <- c("Stress_1","Stress_2","Stress_3","Control_1","Control_2","Control_3")
res.use <- 0.5
dim.use <- 30
seed.use <- 88
project.name <- "mdd.liver"

# 路径设置
data_dir_base <- "/public/home/chidm/Workspace/ouyang/06_mdd.liver/公司/singlesample"
out_dir <- "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/doublet-new"
setwd(out_dir)
qc_plot_dir <- file.path(out_dir, "qc_plots")
dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)

# doublet rate设置（参照doubletfinder官网表格）
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
  
  # QC图：小提琴图
  png(filename = file.path(qc_plot_dir, paste0(i, "_before_filter_violin.png")), width = 10, height = 4, units = "in", res = 400)
  p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1) + plot_layout(guides = "collect")
  print(p1)
  dev.off()
  
  # 记录过滤前细胞数
  cells_before <- ncol(obj)
  
  # 过滤低质量细胞
  obj <- subset(obj, subset =
                  nFeature_RNA >= 200 & nFeature_RNA <= 10000 & #不过滤nCount了
                  percent.mt <= 25)
  
  cells_after_qc <- ncol(obj)
  
  # 标准流程
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:dim.use)
  obj <- FindClusters(obj, resolution = res.use)
  obj <- RunUMAP(obj, dims = 1:dim.use)
  
  # DoubletFinder过程
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
  
  # 过滤掉 Doublet_hi，保留 Singlet 和 Doublet_lo
  obj <- subset(obj, subset = DF_hi.lo %in% c("Singlet", "Doublet_lo"))
  cells_after_doublet <- ncol(obj)
  
  # 重新标准化 & 聚类
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, vars.to.regress = "nCount_RNA")
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:dim.use)
  obj <- FindClusters(obj, resolution = res.use)
  obj <- RunUMAP(obj, dims = 1:dim.use)
  
  # 保存 UMAP 图
  png(paste0(out_dir, "/", i, "_UMAP_after_doublet_removal.png"), width = 6, height = 5, units = "in", res = 400)
  p <- DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 0.3)
  print(p)
  dev.off()
  
  # 保存处理后对象
  sample_list_after_doublet[[i]] <- obj
  
  # 生成QC总表
  qc_summary <- rbind(qc_summary, data.frame(
    sample = i,
    cells_before = cells_before,
    cells_after_qc = cells_after_qc,
    cells_after_doublet_removal = cells_after_doublet,
    estimated_doublet_rate = rate
  ))
}
# 保存QC总表
write.csv(qc_summary, file = file.path(out_dir, "QC_summary_all_samples.csv"),  row.names = FALSE)

##### -------- 合并 Seurat 对象 -------- ######
sample_order <- c("Stress_1", "Stress_2", "Stress_3", "Control_1", "Control_2", "Control_3")
sample_list_doubletRM <- list()
# ---- 重命名 barcode 并存入列表 ----
for (i in 1:length(sample_order)) {
  sample <- sample_order[i]
  suffix <- i  # 用1–6作为后缀
  obj <- sample_list_after_doublet[[sample]]
  # 改列名，加后缀
  new_colnames <- paste0(gsub("-1$", "", colnames(obj)), "-", suffix)
  colnames(obj) <- new_colnames
  sample_list_doubletRM[[sample]] <- obj
}

### 2.1 合并去重后的 Seurat 对象，并统一处理 metadata
sample_order <- c("Stress_1", "Stress_2", "Stress_3", "Control_1", "Control_2", "Control_3")
sample_list_doubletRM <- list()

## 重命名 barcode 加后缀
for (i in 1:length(sample_order)) {
  sample <- sample_order[i]
  suffix <- i  # 用1–6作为后缀
  obj <- sample_list_after_doublet[[sample]]
  
  # 修改列名，加后缀
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

# 添加 group 信息
merged_seurat$group <- ifelse(grepl("^Stress", merged_seurat$sample), "CRS", "Control")

save(merged_seurat, file = file.path(out_dir, "mdd_liver_doublet_removed_merged.RData"))#24883 features across 52848 samples

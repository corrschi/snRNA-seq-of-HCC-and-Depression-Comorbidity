#####-------ncount, nfeature, mt distribution----------#####
setwd("/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/doublet-new/qc-plots-chi")
library(Seurat)
library(dplyr)
library(ggplot2)

sample_list <- c( "Stress_1", "Stress_2", "Stress_3","Control_1", "Control_2", "Control_3")
base_path <- "/public/home/chidm/Workspace/ouyang/06_mdd.liver/公司/singlesample"
project.name <- "mdd.liver"


sample_objs <- list()
qc_meta_all <- data.frame()

# loop read data by sample and build Seurat object
for (sample in sample_list) {
  data_dir <- file.path(base_path, sample, "2.Basic_analysis/2.2.filtered_feature_bc_matrix")
  count.data <- Read10X(data.dir = data_dir)
  
  obj <- CreateSeuratObject(counts = count.data, min.cells = 3, min.features = 200, project = project.name)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  
  obj$sample <- sample
  sample_objs[[sample]] <- obj
  
  qc_meta_all <- rbind(qc_meta_all, obj@meta.data)
}

write.csv(qc_meta_all, "QC_all_samples_meta.csv",row.names = F)

# Summarize the QC statistics for each sample
library(dplyr)
qc_summary <- qc_meta_all %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    cell_number = n(),
    nFeature_min = min(nFeature_RNA),
    nFeature_median = median(nFeature_RNA),
    nFeature_max = max(nFeature_RNA),
    nCount_median = median(nCount_RNA),
    percent_mt_median = median(percent.mt)
  )

# save QC summary table
write.csv(qc_summary, "QC_summary_per_sample.csv",row.names = F)

# 可视化 violin plot
pdf("QC_violin_nFeature.pdf", width = 8, height = 5)
ggplot(qc_meta_all, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(trim = TRUE) + theme_minimal() + labs(title = "nFeature_RNA distribution")
dev.off()

pdf("QC_violin_nCount.pdf", width = 8, height = 5)
ggplot(qc_meta_all, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(trim = TRUE) + theme_minimal() + labs(title = "nCount_RNA distribution")
dev.off()

pdf("QC_violin_percentMT.pdf", width = 8, height = 5)
ggplot(qc_meta_all, aes(x = sample, y = percent.mt, fill = sample)) +
  geom_violin(trim = TRUE) + theme_minimal() + labs(title = "percent.mt distribution")
dev.off()


library(ggplot2)

# 图1：nFeature_RNA 分布
# pdf("QC_density_nFeature-basic.pdf", width = 5, height = 5)
# ggplot(qc_meta_all, aes(x = nFeature_RNA)) +
#   geom_density(fill = "lightgreen", color = "black") +
#   theme_minimal() +
#   labs(x = "nFeature RNA", y = NULL)
# dev.off()

windowsFonts(Times = windowsFont("Times New Roman"))  # 注册字体
pdf("sFig5a-QC_density_nFeature.pdf", width = 4, height = 5)
ggplot(qc_meta_all, aes(x = nFeature_RNA)) +
  geom_density(fill = "lightgreen", color = "black") +
  scale_x_continuous(limits = c(0, 10000)) +
  theme_minimal(base_family = "Times") +  # 使用注册名 Times
  labs(x = "nFeature RNA", y = NULL) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(family = "Times",size=16,face = "bold", color = "black"),
    axis.title = element_text(family = "Times",size=14,face = "bold", color = "black"),
    axis.text = element_text(family = "Times",size=14,face = "bold", color = "black")
  )
dev.off()

# 图2：percent.mt 分布
# pdf("QC_density_percentMT-baic.pdf", width = 5, height = 5)
# ggplot(qc_meta_all, aes(x = percent.mt)) +
#   geom_density(fill = "lightgreen", color = "black") +
#   theme_minimal() +
#   labs(x = "Mitoch genes (%)", y = NULL)
# dev.off()

windowsFonts(Times = windowsFont("Times New Roman"))  # 注册字体
pdf("sFig5b-QC_density_percentMT.pdf", width = 4, height = 5)
ggplot(qc_meta_all, aes(x = percent.mt)) +
  geom_density(fill = "lightgreen", color = "black") +
  scale_x_continuous(limits = c(0, 50)) +  # 限制横坐标范围
  theme_minimal(base_family = "Times") +
  labs(x = "Mitoch genes (%)", y = NULL) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # 黑色边框
    plot.title = element_text(family = "Times",size=16,face = "bold", color = "black"),
    axis.title = element_text(family = "Times", size = 14,face = "bold", color = "black"),  # 坐标轴标题字体大小
    axis.text = element_text(family = "Times", size = 14,face = "bold", color = "black")    # 坐标轴刻度字体大小
  )

dev.off()

#### final过滤标准：
obj <- subset(obj, subset = 
                nFeature_RNA >= 200 & nFeature_RNA <= 10000 &
                percent.mt <= 25)  # mt 通常 <25% 肝脏样本可以接受

# 备用查看常用
summary(macro.filter2$percent.mt)
summary(macro.filter2$nFeature_RNA)
summary(macro.filter2$nCount_RNA)
macro.filter2@meta.data[which(macro.filter2$percent.mt>=20 & macro.filter2$percent.mt <= 25), "Anno.macro"]%>%table

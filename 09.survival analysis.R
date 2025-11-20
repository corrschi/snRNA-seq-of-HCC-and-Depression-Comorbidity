#####--------part1 TCGA data process---------#####
library(IOBR)
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(stringr)
library(ggsignif)
library(ggsci)
library(survival)
library(survminer)
library(pROC)

####--------1.prepare TCGA_LIHC expr_matrix + clinical_data-------#####
### 1.1 LIHC expr prepare
load("~/Documents/tcga/TCGA-LIHC/TCGA-LIHC_mRNA.Rdata")
ls()
names(assays(data))
rowdata <- rowData(data) 

# Extract raw count expression data and process
mrna <- data[data@rowRanges$gene_type == "protein_coding",]
count <- assay(mrna,"unstranded")
# tpm <- assay(mrna, "tpm_unstrand") #if tpm
# fpkm <- assay(mrna,"fpkm_unstrand") #if fpkm
count = data.frame(count)
count[1:4,1:2]

colnames(count) <- gsub("\\.", "-", colnames(count))
count[1:4,1:2]

# colnames process
colnames(count) <- sapply(strsplit(colnames(count), "-"), 
                          function(x) paste(x[1:4], collapse = "-"))
count[1:4,1:2]

# Ensembl_ID process
eset <- count %>% rownames_to_column(var = "Ensembl_ID") 
eset$Ensembl_ID <-substring(eset$Ensembl_ID, 1, 15) 
head(eset)[1:4,1:4]

eset_mean <- eset %>%
  group_by(Ensembl_ID) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  as.data.frame() 
eset_mean[1:4,1:4]

# set Ensembl_ID as rownames
rownames(eset_mean) <- eset_mean$Ensembl_ID
eset_mean$Ensembl_ID <- NULL
head(eset_mean)[1:5]

# Integrating into the IOBR workflow: count2tpm()
eset_lihc = count2tpm(countMat = eset_mean,
                      org="hsa",
                      idType = "Ensembl",
                      gene_symbol = "symbol",
                      source = "local") 
eset_lihc[1:5,1:5]
# eset_lihc %>% as.matrix() %>% as.vector() %>% summary()

eset_lihc <- eset_lihc[, colnames(eset_lihc) != "TCGA-DD-AACA-02B"] 
##Note: TCGA-DD-AACA-02A and TCGA-DD-AACA-02B represent two batches (A and B). In this analysis, only batch A will be retained.

colnames(eset_lihc) <- substr(colnames(eset_lihc), 1, 15)
eset_lihc[1:5,1:5]
eset_lihc_clean <- eset_lihc

### 1.2 LIHC clinical_data prepare
pdata <- readRDS(file = "~/Documents/tcga/TCGA-LIHC/TCGA-LIHC.clinical_patient.rds")
colnames(pdata)
head(pdata[,1:6])

pdata <- pdata[,c("bcr_patient_barcode", "stage_event_pathologic_stage",
                  "gender","age_at_initial_pathologic_diagnosis")] #colnames(pdata)有各种临床信息
colnames(pdata) <- c("bcr_patient_barcode", "pathologic_stage","gender","age")
str(pdata)
table(pdata$pathologic_stage)
pdata$pathologic_stage <- as.character(pdata$pathologic_stage)
pdata <- pdata[pdata$pathologic_stage!="",]
pdata <- na.omit(pdata) 
head(pdata)

# stage I 、II、III、IV process
pdata$stage <- pdata$pathologic_stage
pdata$pathologic_stage%>%table %>% names %>% dput()
pdata$stage[grepl("Stage I$|Stage IA$|Stage IB$",pdata$pathologic_stage)] <- "Stage I"
pdata$stage[grepl("Stage II$|Stage IIA$|Stage IIB$",pdata$pathologic_stage)] <- "Stage II"
pdata$stage[grepl("Stage III$|Stage IIIA$|Stage IIIB$|Stage IIIC$",pdata$pathologic_stage)] <- "Stage III"
pdata$stage[grepl("Stage IV$|Stage IVA$|Stage IVB$|Stage IVC$",pdata$pathologic_stage)] <- "Stage IV"
table(pdata$stage)
table(pdata$pathologic_stage,pdata$stage)
pdata$stage <- factor(pdata$stage, levels = c("Stage I","Stage II","Stage III","Stage IV"))
head(pdata) 

# download from xena
survival_data <- read.table(file="~/Documents/tcga/TCGA-LIHC/survival_LIHC_survival_xena.txt", header = T,sep="\t")
rownames(survival_data) <- survival_data$sample 

survival_data$sample%>%length #438sample
survival_data$X_PATIENT%>%unique%>%length #377patient

# pdata and survival_data combine
overlap_survival_data_pdata <-intersect(str_sub(rownames(survival_data), 1,12), 
                                        pdata$bcr_patient_barcode)
overlap_survival_data_pdata%>%length #352patient

kv <- pdata %>% 
  dplyr::filter(bcr_patient_barcode %in% overlap_survival_data_pdata) %>%as.data.frame()
kv%>%head()
kv%>%dim
survival_data %>% head()
survival_data_overlap <- survival_data %>% filter(X_PATIENT %in% overlap_survival_data_pdata)
survival_data_overlap$X_PATIENT%>%unique%>%length ##352patient

survival_data_overlap$gender <- plyr::mapvalues(str_sub(rownames(survival_data_overlap), 1,12),
                                                from = kv$bcr_patient_barcode %>% as.character(),
                                                to = kv$gender %>% as.character())

survival_data_overlap$age <- plyr::mapvalues(str_sub(rownames(survival_data_overlap), 1,12),
                                             from = kv$bcr_patient_barcode %>% as.character(),
                                             to = kv$age %>% as.character())
survival_data_overlap$stage <- plyr::mapvalues(str_sub(rownames(survival_data_overlap), 1,12),
                                               from = kv$bcr_patient_barcode %>% as.character(),
                                               to = kv$stage %>% as.character())

survival_data_overlap %>% head #352patient，403sample
survival_data_clean <- survival_data_overlap

# save(eset_lihc_clean,survival_data_clean,file='~/Workspace/ouyang/06_mdd.liver/IOBR/LIHC/1_TCGA_epxr_and_survival_data_clean.RData')

####---------part2 EPIC part--------#####
library(Seurat) # v5.0.0
library(EPIC)
library(dplyr)
setwd("~/Workspace/ouyang/06_mdd.liver/paper/tcga")

#1. load TCGA data
load('./IOBR/LIHC/1_TCGA_epxr_and_survival_data_clean.RData')
survival_data_clean%>%head  
eset_lihc_clean[1:5,1:5]  

#2. load scRNA data
load("./epi.RData")
epi@active.ident%>%table

load("./macro.RData")
macro@active.ident%>%table

####------2.Estimating cell proportions using EPIC.-------######
###------2.1 extract signature matrix-------####
### 2.1.1 epi_c8
epi_c8_obj <- subset(epi.filter1, idents = "Epi.c8")
avg_expr <- AverageExpression(object = epi_c8_obj)$RNA
colnames(avg_expr)[1] <- "epi_c8" 

path_epi.filter1_DEGs <- "./epi.filter1_res0.5_dim30_culster_all_DEGs.csv"
marker_epi.filter1 <- read.csv(path_epi.filter1_DEGs,header = T)

marker_epi.c8 <- marker_epi.filter1 %>% filter(cluster=='8')  %>% head(n=10)
mouse_marker.genes <- c("Epcam","Krt18","Chrna9",marker_epi.c8$gene)

# Step 2: Mouse-to-human gene name mapping (using an offline method)
library(biomaRt)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
entrez_mouse <- mapIds(org.Mm.eg.db, keys = mouse_marker.genes, keytype = "SYMBOL", column = "ENTREZID")
human_symbols <- mapIds(org.Hs.eg.db, keys = entrez_mouse, keytype = "ENTREZID", column = "SYMBOL")
human_genes <- unique(na.omit(human_symbols))
epi.human_marker.genes <- human_genes
# load("epi.human_marker.genes.RData")

rownames(avg_expr) <- toupper(rownames(avg_expr))
signature_matrix <- avg_expr[rownames(avg_expr) %in% epi.human_marker.genes, , drop = FALSE]
write.csv(signature_matrix, "signature_epi_c8.csv", quote = FALSE)

### 2.1.2 macro
macro_filter2_obj<- subset(macro, idents = "Macro_Cxcl1")
avg_expr <- AverageExpression(object = macro_filter2_obj)$RNA
colnames(avg_expr)[1] <- "macro_cxcl1"
#signature matrix of macro
path_macro.cxcl1_DEGs <- "./macro_res0.5_dim30_culster_all_DEGs.csv"
marker_macro <- read.csv(path_macro.cxcl1_DEGs,header = T)
marker_macro <- marker_macro   %>% filter(cluster=='3')  %>% head(n=10) #c3 is macro_cxcl1
mouse_marker.genes  <- c("Cd68","Cd14","C1qa","C1qc",marker_macro$gene)

entrez_mouse <- mapIds(org.Mm.eg.db, keys = mouse_marker.genes, keytype = "SYMBOL", column = "ENTREZID")
human_symbols <- mapIds(org.Hs.eg.db, keys = entrez_mouse, keytype = "ENTREZID", column = "SYMBOL")
human_genes <- unique(na.omit(human_symbols))
macro.human_marker.genes <- human_genes
# load("macro.human_marker.genes.RData")

rownames(avg_expr) <- toupper(rownames(avg_expr))
signature_matrix <- avg_expr[rownames(avg_expr) %in% macro.human_marker.genes , drop = FALSE]
write.csv(signature_matrix, "signature_macro_cxcl1.csv", quote = FALSE)


###------2.2 环境B（deseq2-EPIC）推算细胞比例-------####
# 自定义函数：完成自定义signature后EPIC推算bulk-RNAseq细胞比例推算
run_EPIC_with_signature <- function(signature_file, bulk_matrix) {
  sig_mat <- read.csv(signature_file, row.names = 1)
  rownames(sig_mat) <- toupper(rownames(sig_mat))
  common_genes <- intersect(rownames(sig_mat), rownames(bulk_matrix))
  sig_mat <- sig_mat[common_genes, , drop = FALSE]
  bulk <- bulk_matrix[common_genes, ]
  
  reference <- list(
    refProfiles = sig_mat,
    sigGenes = rownames(sig_mat)
  )
  
  out <- EPIC(bulk = as.matrix(bulk), reference = reference)
  celltype_name <- gsub("signature_|\\.csv", "", basename(signature_file))
  res <- as.data.frame(out$cellFractions)
  
  if (celltype_name %in% colnames(res)) {
    res <- res[, celltype_name, drop = FALSE]
  } else {
    res <- res[, 1, drop = FALSE] 
    colnames(res) <- celltype_name
  }
  
  res$SampleID <- rownames(res)
  return(res)
}

signature_files <- c( "signature_epi_c8.csv", "signature_macro_cxcl1.csv")
results_list <- lapply(signature_files, 
                       run_EPIC_with_signature, 
                       bulk_matrix = eset_lihc_clean)
final_result <- Reduce(function(x, y) full_join(x, y, by = "SampleID"), 
                       results_list) %>% select(SampleID, epi_c8, macro_cxcl1)
head(final_result)
write.csv(final_result, "EPIC_custom_signatures_fraction-epi.c8_macro.cxcl1.csv", row.names = FALSE)

####------3.ligand and Receptor Score -------######
cell_fraction <- read.csv(paste0("./EPIC_custom_signatures_fraction-epi.c8_macro.cxcl1.csv"),header = T)
cell_fraction%>%head

ligand_genes <- c("COL4A2", "COL4A1", "COL1A2", "COL1A1")
receptor_genes <- c("SDC4", "SDC1", "CD44", "ITGA9", "ITGB1")
all_genes <- unique(c(ligand_genes, receptor_genes))

expr_subset <- eset_lihc_clean[intersect(all_genes, rownames(eset_lihc_clean)), ]
expr_df <- as.data.frame(t(expr_subset))  # 样本 × 基因格式

df <- merge(expr_df, cell_fraction, by.x = "row.names", by.y = "SampleID")
rownames(df) <- df$Row.names; df$Row.names <- NULL

## Step 2：calculate interaction score
pairs <- list(
  c("COL4A2", "SDC4"), c("COL4A2", "SDC1"), c("COL4A2", "CD44"), c("COL4A2", "ITGA9"), c("COL4A2", "ITGB1"),
  c("COL4A1", "SDC4"), c("COL4A1", "SDC1"), c("COL4A1", "CD44"), c("COL4A1", "ITGA9"), c("COL4A1", "ITGB1"),
  c("COL1A2", "SDC1"), c("COL1A2", "SDC4"), c("COL1A2", "CD44"), c("COL1A2", "ITGA9"), c("COL1A2", "ITGB1"),
  c("COL1A1", "SDC4")
)


calculate_score_top3_sum <- function(row) {
  scores <- sapply(pairs, function(pair) {
    ligand <- pair[1]
    receptor <- pair[2]
    
    ligand_expr <- if (!is.na(row[[ligand]])) row[[ligand]] * row[["epi_c8"]] else 0
    receptor_expr <- if (!is.na(row[[receptor]])) row[[receptor]] * row[["macro_cxcl1"]] else 0
    
    return(ligand_expr * receptor_expr)
  })
  
  top3 <- sort(scores, decreasing = TRUE)[1:3]
  return(sum(top3))
}


df$interaction_score <- apply(df, 1, calculate_score_top3_sum )
df$interaction_score <- log2(df$interaction_score + 1)


df$SampleID <- rownames(df)
df%>%head
## Step 3：survival data
colnames(survival_data_clean)[colnames(survival_data_clean) == "sample"] <- "SampleID"
# only keep tumor
survival_data_clean <- survival_data_clean[grepl("-01$", survival_data_clean$SampleID), ]
common_samples <- intersect(survival_data_clean$SampleID, df$SampleID)
df_use <- df[df$SampleID %in% common_samples, ]
surv_use <- survival_data_clean[survival_data_clean$SampleID %in% common_samples, ]
merged_data <- merge(surv_use, df_use, by = "SampleID")
merged_data %>% head()
####------4.survial -------######

library(survival)
library(survminer)
merged_data$OS_months <- merged_data$OS.time / 30 
res.cut <- surv_cutpoint(merged_data, time = "OS_months", event = "OS", variables = "interaction_score")
cutoff <- res.cut$cutpoint$cutpoint
merged_data$group <- ifelse(merged_data$interaction_score >= cutoff, "High", "Low") 

surv_obj <- Surv(time = merged_data$OS_months, event = merged_data$OS)

fit <- survfit(surv_obj ~ group, data = merged_data)
# plot
png(paste0("1.LIHC-Survival-Months-result-score.png"), height = 7, width = 5, units = "in", res = 400)
p <- ggsurvplot(
  fit,
  data = merged_data,
  pval = TRUE,
  risk.table = TRUE,
  xlab = "Time (months)",
  palette = c("#E41A1C", "#377EB8"),
  title = "Epi.c8-Macro.cxcl1 interaciton score",
  legend.title = "Score",
  legend.labs = c("High", "Low")
)
print(p)
dev.off()


# Forest plot
library(survival)
library(broom)
library(ggplot2)
merged_data$group <- relevel(factor(merged_data$group), ref = "Low") ## 注意设置比较顺序

merged_data$age_group <- ifelse(merged_data$age >= 65, "≥65", "<65")
merged_data$age_group <- factor(merged_data$age_group, levels = c("<65", "≥65"))

merged_data$stage_group <- ifelse(merged_data$stage %in% c("Stage III", "Stage IV"), "T3+T4", "T1+T2")
merged_data$stage_group <- factor(merged_data$stage_group, levels = c("T1+T2", "T3+T4"))

cox_model <- coxph(Surv(OS_months, OS) ~ group + age_group + gender + stage_group, data = merged_data)

df <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
df <- subset(df, term != "(Intercept)")

df$term <- gsub("grouphigh", "Group: score high vs low", df$term)
df$term <- gsub("genderMALE", "Gender: male", df$term)
df$term <- gsub("age_group≥65", "Age ≥65", df$term)
df$term <- gsub("stage_groupT3+T4", "Stage: T3+T4 vs T1+T2", df$term)

df$HR_label <- sprintf("%.1f", df$estimate)
df$P_label <- ifelse(df$p.value < 0.0001, "<0.0001", sprintf("%.4f", df$p.value))
df$label <- paste(df$P_label, df$HR_label)

df$term <- factor(df$term, levels = rev(df$term))

pdf("2.LIHC-Cox-ForestPlot-with-Stage.pdf", width = 5.5, height = 4)
ggplot(df, aes(x = estimate, y = term)) +
  geom_point(shape = 15, size = 3.5, color = "#1b7837") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "solid", color = "skyblue", linewidth = 0.8) +
  geom_text(aes(x = 3.5, label = label), hjust = 0, size = 4.0) +  # 关键：x=2.5
  scale_x_continuous(trans = "log10", limits = c(0.5, 10), breaks = c(1, 2, 4, 8)) +
  labs(x = "Hazard ratio", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )
dev.off()

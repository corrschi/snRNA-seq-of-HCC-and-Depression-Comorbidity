library(org.Mm.eg.db) 
library(AnnotationDbi)
library(AnnoProbe)
library(infercnv)
library(Seurat)
library(gplots)
library(ggplot2)
library(dplyr)
library(harmony)
library(ggsci)


load("./mdd.filter1.RData")
p <- DimPlot(mdd.filter1,pt.size = 0.05,label=T,label.size=4,cols = color_all,raster=F)
ggsave("mdd.filter1-UMAP_checkdata.png", p, width = 9, height = 5)

epi <- subset(mdd.filter1,idents="Epithelial cell")
## subset tcell
tcell <- subset(mdd.filter1,idents=c("T cell"))

# t cell subcluster
dim.use <- 30
res.use <- 0.5
seed.use <- 88
object.name <- "tcell" #修改

MERGE <- tcell
# MERGE <- NormalizeData(object = MERGE, normalization.method = "LogNormalize", scale.factor = 1e4)
MERGE <- FindVariableFeatures(object = MERGE, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
MERGE <- ScaleData(object = MERGE, features = rownames(x = MERGE), 
                   # vars.to.regress = c("nCount_RNA"), 
                   verbose = T)
MERGE <- RunPCA(object = MERGE, features = VariableFeatures(object = MERGE), verbose = FALSE)

png(paste0(object.name,"_RunHarmony_idx.png"), width = 10, height = 10, res = 400, units = "in")
MERGE <- RunHarmony(MERGE,c("sample"), plot_convergence = TRUE,verbose = FALSE) #最后按sample进行harmony
dev.off()

MERGE <- FindNeighbors(MERGE, reduction = "harmony", verbose = FALSE, dims = 1:dim.use)
MERGE <- FindClusters(MERGE, resolution = res.use, verbose = FALSE, random.seed = seed.use) 

MERGE <- RunUMAP(MERGE, dims = 1:dim.use, reduction = "harmony",umap.method = "uwot",
                 n.neighbors = 20L, min.dist = 0.5)

tcell <- MERGE

p <- DimPlot(tcell,pt.size = 0.05,label=T,label.size=4,cols = color_used,raster=F)
ggsave("tcell-UMAP_reharmony_num.png", p, width = 9, height = 5)


p <- FeaturePlot(tcell,features = c("Cd3e","Mki67" ), 
                 pt.size = 0.1,ncol = 2, cols = c("lightgrey", "#FF0000"),raster = FALSE)
p <- p & theme(
  panel.border = element_rect(color = "black", size = 1, fill = NA),
  panel.background = element_blank(),
  strip.background = element_blank(),
  panel.grid = element_blank(),
  strip.text = element_text(family = "Times", face = "bold", color = "black", size = 12))
ggsave("tcell-Featureplot-Mki67.png", p, width = 8, height =3)

tcell <- subset(tcell,idents=c("0","2","4","5","7","12")) # remove proliferative t cells

sce <- merge(epi,tcell)

# 2.read data
table(Idents(sce))
table(sce@meta.data$cluster)
sce$orig.ident <- sce@meta.data %>% rownames()
sce@meta.data$patient <- sce@meta.data$sample
table(sce@meta.data$patient)

# expression matrix
dat <- GetAssayData(sce,layer = 'counts',assay = 'RNA')
dat[1:4,1:4]     
## cell type
groupinfo <- data.frame( v1 = colnames(dat), v2 = sce@meta.data$cluster,
                         check.names = FALSE) 
groupinfo$v2%>%table #check data 
## genomic coordinates
library(org.Mm.eg.db) 
library(AnnotationDbi)
library(AnnoProbe)
geneInfor <- annoGene(rownames(dat),"SYMBOL","mouse") 
# Retained mouse autosomes（1–19）
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
colnames(geneInfor)
geneInfor <- geneInfor %>% arrange(chr_num, start) %>% select(SYMBOL, chr, start, end)
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]

length(unique(geneInfor[,1]))
head(geneInfor)

# 3.ranking
dat%>%dim
dat <- dat[rownames(dat) %in% geneInfor[,1],]
dat <- dat[match(geneInfor[,1],rownames(dat)),]
# check data
dim(dat)
head(groupinfo)
dat[1:4,1:4]
table(groupinfo$v2)
dim(groupinfo)

# 4.save files
expFile <- 'expFile.txt' 
colnames(dat) <- gsub("-", "_", colnames(dat))
write.table(dat, file = expFile, sep = '\t', quote = F)# for a long time

groupFiles <- 'groupFiles.txt'
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
write.table(groupinfo,file = groupFiles, sep = '\t',quote = F, col.names = F, row.names = F)

head(geneInfor)
geneFile <- 'geneFile.txt'
write.table(geneInfor, file = geneFile, sep = '\t',quote = F, col.names = F, row.names = F)

# 5. inferCNV process
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = "expFile.txt", 
                                     annotations_file  = groupFiles,
                                     delim             = "\t",
                                     gene_order_file   = "geneFile.txt",
                                     ref_group_names   = c("T cell"))  

infercnv_obj2 <- infercnv::run(infercnv_obj,
                               cutoff            = 0.1,
                               out_dir           = "infercnv_output",
                               cluster_by_groups = TRUE,       
                               hclust_method     = "ward.D2",
                               analysis_mode     = "samples",   
                               denoise           = TRUE,
                               HMM               = FALSE,       
                               plot_steps        = FALSE,
                               num_threads       = 4)
# read in inferCNV object
infer_CNV_obj<-readRDS("./infercnv_output/run.final.infercnv_obj")
expr <- infer_CNV_obj@expr.data
expr[1:4,1:4]

data_cnv <- as.data.frame(expr)
dim(expr)
colnames(data_cnv)
rownames(data_cnv)
colnames(sce) <- gsub("-","_",colnames(sce))
meta <- sce@meta.data

##################3.distinguish malignant and non-malignant cells############
library(dplyr)
library(mclust)   
cnv_mat <- as.matrix(infer_CNV_obj@expr.data)   # genes x cells
stopifnot(all(gsub("-", "_", colnames(sce)) %in% colnames(cnv_mat) |
                colnames(cnv_mat) %in% gsub("-", "_", colnames(sce))))

colnames(sce) <- gsub("-", "_", colnames(sce))

cnv_score <- colMeans(abs(cnv_mat - 1))
cnv_df <- data.frame(cell = names(cnv_score), cnv_score, row.names = NULL)

fit <- Mclust(log1p(cnv_df$cnv_score), G = 2)
hi_class <- which.max(fit$parameters$mean)      
cnv_df$CNV_call <- ifelse(fit$classification == hi_class, 
                          "CandidateMalig", "CandidateNonMalig")

sce$CNV_score <- cnv_df$cnv_score[match(colnames(sce), cnv_df$cell)]
sce$CNV_call  <- cnv_df$CNV_call[ match(colnames(sce), cnv_df$cell)]

epi <- subset(sce, idents = "Epithelial cell")

# 6)  tSNE/UMAP
p <- DimPlot(epi, reduction = "umap", group.by = "CNV_call",pt.size=0.5,
             cols = c("CandidateMalig"    = "#49C39E",
                      "CandidateNonMalig" = "#F07F5A")) +
  labs(color = "Result Type")
ggsave("epi_CNV_malignancy_umap.png", p, width = 6, height = 4, dpi = 300)

epi.metadata <- epi@meta.data[,c("orig.ident", "group", "cluster","subcluster", "sample",  "patient", "CNV_score", "CNV_call")]
write.csv(epi.metadata,file="epi.CNV.csv",row.names = F)
epi$CNV_call %>% table () 

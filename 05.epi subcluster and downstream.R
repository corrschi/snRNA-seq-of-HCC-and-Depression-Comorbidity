####------------epithelial subcluster --------######
setwd("~/Workspace/ouyang/06_mdd.liver/paper/clustering-epi")
epi<- subset(mdd.liver,idents="Epithelial cell")

dim.use <- 30
res.use <- 0.5
seed.use <- 88
object.name <- "epi" #修改
MERGE <- epi
# MERGE <- NormalizeData(object = MERGE, normalization.method = "LogNormalize", scale.factor = 1e4)
MERGE <- FindVariableFeatures(object = MERGE, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
MERGE <- ScaleData(object = MERGE, features = rownames(x = MERGE), 
                   verbose = T)
MERGE <- RunPCA(object = MERGE, features = VariableFeatures(object = MERGE), verbose = FALSE)
png(paste0(object.name,"_RunHarmony_idx.png"), width = 10, height = 10, res = 400, units = "in")
MERGE <- RunHarmony(MERGE,c("sample"), plot_convergence = TRUE,verbose = FALSE)
dev.off()

MERGE <- FindNeighbors(MERGE, reduction = "harmony", verbose = FALSE, dims = 1:dim.use)
MERGE <- FindClusters(MERGE, resolution = res.use, verbose = FALSE, random.seed = seed.use) 

MERGE <- RunUMAP(MERGE, dims = 1:dim.use, reduction = "harmony",umap.method = "uwot",
                 n.neighbors = 20L, min.dist = 0.5)


epi <- MERGE

## Anno.epi修改
epi$Anno.epi <- paste0("Epi.c",epi$seurat_clusters)
epi$Anno.epi <- factor(epi$Anno.epi,levels = paste0("Epi.c",c(0:12)))
epi$cluster <- epi$Anno.epi 
Idents(epi) <- "Anno.epi"

## group修改
epi.filter1$group <- as.character(epi.filter1$group)  # 确保是字符向量
epi.filter1$group[epi.filter1$group == "Stress"] <- "CRS"
epi.filter1$group <- factor(epi.filter1$group, levels = c("Control", "CRS")) 

#epi.filter1 meta.data重新存储
# save
subset_cells <- epi.filter1
subset_cells@assays$RNA@data <- as.matrix(0)
subset_cells@assays$RNA@scale.data <- as.matrix(0)
save(subset_cells, file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/epi.filter1.RData")

meta.data <- epi.filter1@meta.data
save(meta.data,file="/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/epi.filter1_meta.data.RData")

#Plots
# Fig5e
png(paste0("Fig5e-",object.name,"_res", res.use ,"_dim", dim.use ,"_all_by_harmony_dimplot-use.png"), width = 7, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.2,label=T,label.size=5,cols = color_epi,raster=F)
dev.off()

# sFig5a
p <- FeaturePlot(MERGE, 
                 features = c("Krt8", "Alb", "Epcam"),
                 pt.size = 0.1,
                 ncol = 3,
                 cols = c("lightgrey", "#FF0000"),
                 raster = FALSE)

p <- p + theme(
  plot.title = element_text(family = "Times New Roman", face = "bold", size = 14),
  axis.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  axis.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  strip.text = element_text(family = "Times New Roman", face = "bold", size = 12)) 
ggsave("sFig5a.FeaturePlote-epi.Marker.png", plot = p, width = 12, height = 3, units = "in", dpi = 400)


# Fig6f
library(Seurat)
library(Nebulosa)
png("Fig6f.密度图-Chrna9.png",  width = 5, height = 4, res = 400, units = "in")
plot_density(MERGE,"Chrna9")
dev.off()

####-------irGSEA----------######
setwd("~Workspace/ouyang/06_mdd.liver/paper/clustering-epi/irGSEA")

library(Seurat)
library(irGSEA)
library(tidyverse)

pbmc3k.final<- epi

pbmc3k.final$Anno.epi <- paste0("Epi.c", pbmc3k.final$seurat_clusters)
pbmc3k.final$Anno.epi <- factor(pbmc3k.final$Anno.epi,levels = paste0("Epi.c",c(0:12)))

Idents(pbmc3k.final) <- 'Anno.epi'
# remove Epi.c12 because cell number too small
pbmc3k.final <- subset(pbmc3k.final,idents="Epi.c12",invert=T)
pbmc3k.final@active.ident%>%table

pbmc3k.final$seurat_annotations <- pbmc3k.final$Anno.epi


#2 Calculate enrichment scores 
pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", 
                             slot = "data", 
                             seeds = 123, 
                             ncores = 4,
                             min.cells = 3,
                             min.feature = 0,
                             custom = F, geneset = NULL, 
                             msigdb = T, 
                             species = "Mus musculus",  
                             category = "C2",             
                             subcategory = "CP:KEGG",     
                             geneid = "symbol",      
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea", "JASMINE", "viper"), 
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')     

Seurat::Assays(pbmc3k.final)

#3 Integrate differential gene set
result.dge <- irGSEA.integrate(object = pbmc3k.final,
                               group.by = "seurat_annotations",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

class(result.dge)

save(pbmc3k.final,file="./pbmc3k.final-KEGG.RData")
save(result.dge,file = "./result.dge-KEGG.RData")

#4 可视化
# load("./pbmc3k.final-KEGG.RData")
# load("./result.dge-KEGG.RData")

# sFig5e:ridgeplot
ridgeplot <- irGSEA.ridgeplot(
  object = pbmc3k.final,
  method = "singscore",
  show.geneset = "KEGG-NEUROACTIVE-LIGAND-RECEPTOR-INTERACTION"
)

ridgeplot <- ridgeplot + theme(
  plot.title = element_text(family = "Times New Roman", face = "bold", size = 14),
  axis.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  axis.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  strip.text = element_text(family = "Times New Roman", face = "bold", size = 12))
png("sFig5e.irGSEA.ridgeplot-KEGG-NEUROACTIVE-LIGAND-RECEPTOR-INTERACTION.png", 
    width = 7, height = 7, units = "in", res = 400)
print(ridgeplot)
dev.off()

# Fig6g
library(Seurat)
library(ggplot2)
library(patchwork)

# 设置 theme
vln_theme <- theme(
  plot.title = element_text(hjust = 0, vjust = 1, size = 12, family = "Times New Roman", face = "bold"),
  axis.title.x = element_text(size = 10, family = "Times New Roman"),
  axis.title.y = element_text(size = 10, family = "Times New Roman"),
  axis.text.x = element_text(size = 12, family = "Times New Roman"), #x轴坐标分组
  axis.text.y = element_text(size = 10, family = "Times New Roman"),
  legend.text = element_text(size = 10, family = "Times New Roman"),
  legend.title = element_text(size = 12, family = "Times New Roman")
) &
  plot_layout(guides = "collect") &
  plot_annotation(title = NULL)


dat_plot <- subset(epi.filter1, idents = "Epi.c8")
dat_plot$group <- factor(dat_plot$group, levels = c("Control", "CRS"))
Idents(dat_plot) <- 'group'

png("Fig6g-Epi.c8-Vln-Chrna9.png", width = 3, height = 3, res = 400, units = "in")
p <- VlnPlot(dat_plot,
             features = c("Chrna9"),
             pt.size = 0.001,
             ncol = 1)
p <- p & vln_theme  
print(p)
dev.off()

####----------sFig6i ClusterGVis for epi.filter1 Heatmap---------#######
# load data
setwd("~/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/ClusterGVis")
library(ClusterGVis)
library(Seurat)
pbmc <-epi
## data check
new.cluster.ids <- pbmc@active.ident%>%table%>%names
new.cluster.ids

p <- DimPlot(pbmc, reduction = "umap",label = T,label.size=6) + NoLegend()
ggsave("1.seurat-umap-check.png", plot = p, width = 6, height = 6, units = "in", dpi = 400)

# find markers for every cluster compared to all remaining cells
# report only the positive ones
pbmc.markers.all <- Seurat::FindAllMarkers(pbmc,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

# # get top 10 genes
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = pbmc,
                                diffData = pbmc.markers,
                                showAverage = TRUE)
str(st.data)

# add gene name
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = F)]

# line plot
p<- visCluster(object = st.data,
               plot.type = "line")
ggsave(paste0( "sFig5d.Line plot.png"),p,width = 12,height = 9)


# heatmap plot
pdf('Heatmap.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:13))#调整聚类顺序。
dev.off()




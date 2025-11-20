##----------------------------------------------
## cellchat of epi.c8 with immune cells related sFig6b
##----------------------------------------------
# Seurat4 env of conda
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
ptm = Sys.time()
library(ComplexHeatmap)
# rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("~/Workspace/ouyang/06_mdd.liver/paper/clustering-total/cellchat")
setwd("./epi.c8-immune")

# load mdd.filter1
load("~/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver.filter1/mdd.filter1.RData")
load("~/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver.filter1/mdd.filter1_meta.data.RData")

DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))

subset_cells@meta.data <- meta.data
scRNA <- subset_cells #24286 features across 50336 samples
scRNA$celltype <- scRNA$Anno.chi 
table(scRNA$celltype)

p <- DimPlot(scRNA,pt.size = 0.05,label=T,label.size=4,cols = color_all,raster=F)
ggsave("mdd.filter1-UMAP_checkdata.png", p, width = 9, height = 5)

# load epi.filter1
load("~/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/epi.filter1.RData")
load("~/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/epi.filter1_meta.data.RData")

DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))

subset_cells@meta.data <- meta.data
epi.filter1 <- subset_cells #24286 features across 5558 samples
epi.filter1 @active.ident%>%table
p <- DimPlot(epi.filter1,pt.size = 0.05,label=T,label.size=4,cols = color_epi,raster=F)
ggsave("epi.filter1-UMAP_checkdata.png", p, width = 6, height = 5)


epi.c8_barcode <- subset(epi.filter1,idents=("epi.c8"))%>%colnames #208 cells
all.immune_barcode <- subset(scRNA,idents=c("Epithelial cell","Fibroblast",  "Endothelial cell"),
                                           invert = T)%>%colnames
barcodes <- c(epi.c8_barcode,all.immune_barcode) #41872 cells

scRNA <- subset(scRNA,cells = barcodes) 
scRNA@active.ident%>%table

scRNA$celltype <-  as.character(scRNA$celltype)
scRNA$celltype%>%str

scRNA@meta.data$celltype[scRNA@meta.data$celltype == "Epithelial cell"] <- "Epi.c8"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "T cell"] <- "T"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "NK cell"] <- "NK"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "Mono/Mac cell"] <- "Mono_Mac"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "Neutrophil cell"] <- "Neutro"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "B cell"] <- "B_plasma"
scRNA@meta.data$celltype[scRNA@meta.data$celltype == "Plasma cell"] <- "B_plasma"

Idents(scRNA) <- 'celltype'
scRNA@active.ident%>%table
table(scRNA$celltype,scRNA$group)
table(scRNA$group)

seurat.Control <- subset(scRNA, subset=group=="Control") #NL
seurat.CRS <- subset(scRNA, subset=group=="CRS") #LS

### 创建cellchat对象
cellchat.Control <- createCellChat(seurat.Control@assays$RNA@data, meta = seurat.Control@meta.data, group.by = "celltype")
cellchat.CRS <- createCellChat(seurat.CRS@assays$RNA@data, meta = seurat.CRS@meta.data, group.by = "celltype")

# 2.2 Control的细胞通讯网络
cellchat <- cellchat.Control
cellchat@DB <- CellChatDB.mouse #小鼠
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_Control.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp<- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_Control.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.Control <- cellchat
save(cellchat.Control, file="cellchat.Control.RData")


# 2.3 CRS的细胞通讯网络
cellchat <- cellchat.CRS
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_CRS.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_CRS.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.CRS <- cellchat
save(cellchat.CRS, file="cellchat.CRS.RData")

# 2.4 合并cellchat对象
load("cellchat.Control.RData")
load("cellchat.CRS.RData")
object.list <- list(Control = cellchat.Control, CRS = cellchat.CRS)
# object.list <- list(NL = cellchat.NL, LS = cellchat.LS)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

######---------- 3. 可视化--------#############
pdf("sFig6b.netVisual_number_strength.pdf",width = 8,height =5)
par(mfrow = c(1,2))
p1 <- netVisual_diffInteraction(cellchat, weight.scale = T)
p2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
print(p)
dev.off()


##----------------------------------------------------------
## macro.filter2 plots-Cellchat related Fig6e and sFig6d,e
##----------------------------------------------------------
setwd("~Workspace/ouyang/06_mdd.liver/paper/clustering-macro/figures/cellchat")
setwd("./epi.c8-macro.allsub")

## 创建cellchat对象
# 加载数据
total_path <- "~Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver.filter1/"
epi_path <- "~Workspace/ouyang/06_mdd.liver/paper/clustering-epi/"
macro_path <- "~Workspace/ouyang/06_mdd.liver/paper/clustering-macro/"
# load mdd.filter1
load(paste0(total_path,"mdd.filter1.RData"))
load(paste0(total_path,"mdd.filter1_meta.data.RData"))#meta.data
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))
subset_cells@meta.data <- meta.data

scRNA <- subset_cells #24286 features across 50336 samples
scRNA$celltype <- scRNA$Anno.chi 
table(scRNA$celltype)

# load epi
load(paste0(epi_path,"epi.filter1.RData"))
load(paste0(epi_path,"epi.filter1_meta.data.RData"))#meta.data
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))
subset_cells@meta.data <- meta.data

epi.filter1 <- subset_cells #24286 features across 5558 samples
epi.filter1$celltype <- epi.filter1$Anno.epi
table(epi.filter1$celltype)

# load macro
load(paste0(macro_path,"macro.filter2.RData"))
load(paste0(macro_path,"macro.filter2_meta.data.RData"))#meta.data
DefaultAssay(subset_cells) <- "RNA"
subset_cells@assays$RNA@data <- subset_cells@assays$RNA@counts
subset_cells <- NormalizeData(subset_cells)
subset_cells <- ScaleData(subset_cells, features=VariableFeatures(subset_cells))
subset_cells@meta.data <- meta.data

macro.filter2 <- subset_cells #24286 features across 52020 samples
macro.filter2$celltype <- macro.filter2$Anno.macro
table(macro.filter2$celltype)

epi.c8_barcode <- subset(epi.filter1,idents=("epi.c8"))%>%colnames
barcodes <- c(epi.c8_barcode,colnames(macro.filter2)) 
scRNA <- subset(scRNA,cells = barcodes) #24286 features across 12285 samples

## kv for mapvalue
scRNA$celltype %>% table
macro.filter2$celltype <- macro.filter2$Anno.macro %>% as.character()

epi.c8 <- subset(epi.filter1,idents=("Epi.c8"))
epi.c8$celltype <-"Epi.c8"

kv <- rbind(macro.filter2@meta.data%>%select(celltype),epi.c8@meta.data %>% select(celltype)) #12285 samples
kv$barcode <- rownames(kv)

scRNA$celltype <- plyr::mapvalues(scRNA %>% colnames(),
                                  from = kv$barcode,
                                  to = kv$celltype)
scRNA$celltype %>% table
Idents(scRNA) <- 'celltype'
table(scRNA$celltype,scRNA$group)
table(scRNA$group)

## 2.1 cellchat pipline
seurat.Control <- subset(scRNA, subset=group=="Control") #NL
seurat.CRS <- subset(scRNA, subset=group=="CRS") #LS

cellchat.Control <- createCellChat(seurat.Control@assays$RNA@data, meta = seurat.Control@meta.data, group.by = "celltype")
cellchat.CRS <- createCellChat(seurat.CRS@assays$RNA@data, meta = seurat.CRS@meta.data, group.by = "celltype")

# 2.2 Control的细胞通讯网络
cellchat <- cellchat.Control
cellchat@DB <- CellChatDB.mouse #小鼠
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_Control.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp<- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_Control.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.Control <- cellchat
save(cellchat.Control, file="cellchat.Control.RData")


# 2.3 CRS的细胞通讯网络
cellchat <- cellchat.CRS
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr_CRS.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway_CRS.csv")

cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.CRS <- cellchat
save(cellchat.CRS, file="cellchat.CRS.RData")

# 2.4 合并cellchat对象
load("cellchat.Control.RData")
load("cellchat.CRS.RData")
object.list <- list(Control = cellchat.Control, CRS = cellchat.CRS)
# object.list <- list(NL = cellchat.NL, LS = cellchat.LS)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# sFig6d
levels(cellchat@idents$joint)
pdf("sFig6d.number_strength-netVisual-epi.c8 and macro.pdf", width = 8, height = 5)
par(mfrow = c(1, 2))

p1 <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "count", sources.use = "Epi.c8")
p2 <- netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", sources.use = "Epi.c8")
dev.off()

# Fig6e
df.net <- subsetCommunication(cellchat, thresh = 0.05)

df.net$CRS$pathway_name%>%unique
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("COLLAGEN"))
pairLR.use %>% head()

p <- netVisual_bubble(cellchat, 
                      sources.use = 1, 
                      targets.use = c(4,5), 
                      comparison = c(1, 2), 
                      max.dataset = 2,  
                      title.name = "Increased signaling in CRS",
                      angle.x = 90, 
                      pairLR.use = pairLR.use,
                      remove.isolate = F)
ggsave("sFig6e-COLLAGEN-LR_bubble-Increased.pdf", p, width = 4, height =5)

## sFig6e
pdf("sFig6e.COLLAGEN_chord.pdf", width = 12, height = 6)
pathways.show <- c("COLLAGEN")
par(mfrow = c(1,2), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "chord", 
                      pt.title = 3, 
                      title.space = 0.05,
                      edge.width.max = 12,
                      vertex.label.cex = 0.9, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


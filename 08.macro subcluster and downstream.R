# macro subcluster and annotation
setwd("~/Workspace/ouyang/06_mdd.liver/paper/clustering-macro")
## load mdd.filter1
load("./mdd.filter1.RData")
mdd.filter1 #24286 features across 50336 samples

macro <- subset(mdd.filter1,idents="Mono_Macro cell")

dim.use <- 30
res.use <- 0.5
seed.use <- 88
object.name <- "macro" #修改

MERGE <- macro
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


macro <- MERGE

#Plots
png(paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_all_by_harmony_dimplot.png"), width = 6, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.2,label=T,label.size=4,cols = color_macro,raster=F)
dev.off() 

#FindMarkers
object.name <- "mdd.liver"
# MERGE.sampling <- subset(MERGE, cells=WhichCells(MERGE, downsample=500, seed = seed.use))#抽样
MERGE.sampling <- MERGE
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGs.csv"),sep = ",",row.names = T,quote = F)
save(MERGE.markers,file = "~/Workspace/ouyang/06_mdd.liver/paper/clustering-total/total_clustering.markers.RData")

#######--mdd.filter2注释--------########
Idents(MERGE) <- "seurat_clusters"
MERGE <- RenameIdents(object = MERGE, 
                      `0` = "Macro_Ccl8", 
                      `1` = "Macro_Spp1",
                      `2` = "Mono_Ly6c2",`5` = "Mono_Ly6c2",
                      `3` = "Macro_Cxcl1",
                      `4` = "Macro_Cxcl9",
                      `6` = "Kuppfle_Clec4f", `11` = "Kuppfle_Clec4f",
                      `7` = "Macro_Stmn1",
                      `8` = "Macro_Pdpk1",
                      `9` = "Mono_S100a9",
                      `10` = "Macro_Cxcl13")
MERGE[["Anno.macro"]] <- Idents(object = MERGE)

color_macro <-color_used
png(paste0("Fig6a.macro_all_by_harmony_dimplot-Anno-Nolabel.png"), width = 7, height = 5, res = 400, units = "in")
DimPlot(MERGE,pt.size = 0.2,label=F,label.size=4,cols = color_macro,raster=F)
dev.off() 

macro <- MERGE
save(macro, file = "./macro.RData")

#####------ Fig6b-----####
windowsFonts(Times = windowsFont("Times New Roman"))
p <- DotPlot(MERGE, 
             features = c("Ccl8","Ccl7",
                          "Spp1","Mmp12",
                          "Ly6c2","Plac8",
                          "Cxcl1","Cxcl2",
                          "Cxcl9","Ccl24",
                          "Clec4f","Timd4",
                          "Stmn1","Hmgb2",
                          "Pdpk1","Pecam1",
                          "S100a9","S100a8",
                          "Cxcl13","Serpinb2"       
             )) +
  theme_bw(base_family = "Times") +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(family = "Times", face = "bold", color = "black", 
                               angle = 0, hjust = 1, vjust = 0.95, size = 8),
    axis.text.x = element_text(family = "Times", face = "bold", color = "black",
                               angle = 45, hjust = 1, vjust = 0.95, size = 8),
    axis.title = element_text(family = "Times", face = "bold", color = "black", size = 10),
    panel.border = element_rect(color = "black", size = 1.5)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c('#FFFFF0', '#FFCC66', '#FF9966', '#FF0000')) +
  scale_size(range = c(0.1, 6)) +
  scale_y_discrete(limits = rev(levels(Idents(MERGE))))
ggsave("Fig7b.Dotplot-Anno-macro分群主Marker.png", p, width = 8, height = 4)

#####------ Fig6c-----####
library(pheatmap)
library(viridis)

dat_plot <- MERGE
cluster.average <- AverageExpression(dat_plot)[[1]] #The return value of AverageExpression is a list!

genes.select <- c("Cxcl2", "Cxcl1","Ccl3","Cxcl3","Ccl4", "Ccl9", "Ccl2", "Ccl6","Ccl5",
                  "Il1a","Il1b","Il6","Tnf")

cluster.average.select <- cluster.average[genes.select,]

cluster.average.select.transposed <- t(cluster.average.select)  
cell_size <- 10 
pdf_width <- nrow(cluster.average.select.transposed) * cell_size / 50 + 5 
pdf_height <- ncol(cluster.average.select.transposed) * cell_size / 100 + 0 

pdf("Fig6c.Heatmap of Chemokines and Inflammatory Cytokines.pdf", width = 9, height = 6)
pheatmap(cluster.average.select.transposed,
         cluster_row = F, cluster_col = FALSE,  
         scale = "column",
         # color = colorRampPalette(c("#6d8bc3", "white", "#e3716e"))(100),
         color = colorRampPalette(c('#000c71','white','#c5121f'))(100), 
         border_color = "white",
         cellwidth = cell_size,  
         cellheight = cell_size,  
         gaps_col = c(9), 
         treeheight_col = 0, treeheight_row = 0,  
         angle_col = 90,  
         fontsize_row = 10,  
         fontsize_col = 10)
dev.off()


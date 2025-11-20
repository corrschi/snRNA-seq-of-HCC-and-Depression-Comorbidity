# ---------------------------------------------
# Seurat harmony and clustering and Annotation
# --------------------------------------------
setwd("~/Workspace/06_mdd.liver/paper/clustering-total")
load(file.path(out_dir, "mdd_liver_doublet_removed_merged.RData"))

dim.use <- 30
res.use <- 0.6
seed.use <- 88
object.name <- "mdd.liver" #修改

MERGE <- merged_seurat

MERGE <- NormalizeData(object = MERGE, normalization.method = "LogNormalize", scale.factor = 1e4)
MERGE <- FindVariableFeatures(object = MERGE, selection.method = 'vst', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
MERGE <- ScaleData(object = MERGE, features = rownames(x = MERGE), 
                   # vars.to.regress = c("nCount_RNA"), 
                   verbose = T)
MERGE <- RunPCA(object = MERGE, features = VariableFeatures(object = MERGE), verbose = FALSE)

png(paste0(object.name,"_RunHarmony_idx.png"), width = 10, height = 10, res = 400, units = "in")
MERGE <- RunHarmony(MERGE,c("sample"), plot_convergence = TRUE,verbose = FALSE) #注意匹配
dev.off()

MERGE <- FindNeighbors(MERGE, reduction = "harmony", verbose = FALSE, dims = 1:dim.use)
MERGE <- FindClusters(MERGE, resolution = res.use, verbose = FALSE, random.seed = seed.use) 

MERGE <- RunUMAP(MERGE, dims = 1:dim.use, umap.method = "uwot",reduction = "harmony",
                 n.neighbors = 20L, min.dist = 0.5)

merged_seurat <- MERGE #修改

#FindMarkers
object.name <- "mdd.liver"
# MERGE.sampling <- MERGE
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGs.csv"),sep = ",",row.names = T,quote = F)
escc.makers <- MERGE.markers
save(escc.makers,file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/total_clustering.markers.RData")


####----------分群常用Marker总结（小鼠）------###########
png(paste0("Feature-main-3-T-sub.png"),width = 16, height = 9, units = "in", res = 400)
FeaturePlot(MERGE, 
            features = c("Cd3d", "Cd8a", "Cd4", "Foxp3", "Il2ra", "Ccr7", "Tcf7", "Lef1", 
                         "Tbx21", "Klrg1", "Itgae", "Nr4a1", "Mki67",  "Pdcd1", "Slc4a10"
            ), 
            pt.size = 0.1,
            ncol=4)
dev.off()

png(paste0("Feature-main-4-kuppfle-all.png"),width = 16, height = 9, units = "in", res = 400)
FeaturePlot(MERGE, 
            features = c("Cd68","C1qa","C1qb","C1qc",
                         "Adgre1","Clec4f","Vsig4","Timd4","Marco","Fcgr1","Spic"
            ), 
            pt.size = 0.1,
            ncol=4)
dev.off()

##########---------Annotation----------###########
Idents(MERGE) <- "seurat_clusters"
MERGE <- RenameIdents(object = MERGE, 
                      `15` = "Fibroblast", 
                      `2` = "Epithelial cell",`13` = "Epithelial cell",`21` = "Epithelial cell",
                      `11` = "Endothelial cell",
                      `9` = "T cell",`1` = "T cell",`4` = "T cell",`8` = "T cell",`17` = "T cell",
                      `23` = "Mix T and Neutrophil cell",
                      `14` = "NK cell",
                      `3` = "Mono/Mac cell",`5` = "Mono/Mac cell",`7` = "Mono/Mac cell",
                      `16` = "Mono/Mac cell",`20` = "Mono/Mac cell",
                      `6` = "Mono/Mac cell",
                      `12` = "cDC",
                      `18` = "pDC",
                      `0` = "Neutrophil cell",`22` = "Neutrophil cell",
                      `10` = "B cell",
                      `19` = "Plasma cell"
)
MERGE[["Anno.chi"]] <- Idents(object = MERGE)

# save
subset_cells <- MERGE
subset_cells@assays$RNA@data <- as.matrix(0)
subset_cells@assays$RNA@scale.data <- as.matrix(0)
save(subset_cells, file = "~/Workspace/06_mdd.liver/paper/clustering-total/mdd.liver-v0/mdd.liver.RData")

meta.data <- MERGE@meta.data
save(meta.data, file = "~/Workspace/06_mdd.liver/paper/clustering-total/mdd.liver-v0/mdd.liver_meta.data.RData")


##########---------Plots----------###########
# Fig5a
p <- DimPlot(MERGE,pt.size = 0.05,label=T,label.size=4,cols = color_all,raster=F)
ggsave("Fig5a.MDD.liver-UMAP-Anno.png", p, width = 9, height = 5)
# p <- DimPlot(MERGE,pt.size = 0.05,label=F,label.size=4,cols = color_all,raster=F)
# ggsave("Fig5a.MDD.liver-UMAP-Anno-noLabel.png", p, width = 8, height = 5)

# sFig5c,d
Idents(MERGE) <- 'group'
p <- DimPlot(MERGE,pt.size = 0.05,label=F,label.size=4,cols = color_list,raster=F)
ggsave("sFig5c.MDD.liver-UMAP-byGroup.png", p, width = 5, height = 4)
Idents(MERGE) <- 'sample'
p <- DimPlot(MERGE,pt.size = 0.05,label=F,label.size=4,cols = color_sample[2:10],raster=F)
ggsave("sFig5d.MDD.liver-UMAP-bySample.png", p, width = 5, height = 4)

# Fig5b
windowsFonts(Times = windowsFont("Times New Roman"))
p <- DotPlot(MERGE, 
             features = c("Col3a1","Col1a1",       # fibroblast
                          "Krt18", "Epcam",        # epithelial
                          "Pecam1", "Eng",         # endothelial
                          "Ptprc",                 # immune
                          "Cd3e", "Cd3d",          # T cells
                          "Klrb1c", "Nkg7",        # NK
                          "Cd68", "Cd14", "C1qa", "C1qc",  # mono/macro
                          "H2-Ab1", "Cst3",        # cDC
                          "Siglech", "Bst2",       # pDC
                          "Csf3r", "Clec4d",       # neutrophil
                          "Cd79a", "Ms4a1",        # B cells
                          "Jchain", "Mzb1"         # plasma
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
ggsave("Fig5b.mdd.filter2-Dotplot-Anno-总体分群主Marker.png", p, width = 10, height = 4)
ggsave("Fig5b.mdd.filter2-Dotplot-Anno-总体分群主Marker.pdf", p, width = 10, height = 4)

# sFig5e-FeaturePlot
library(Seurat)
library(ggplot2)
windowsFonts(Times = windowsFont("Times New Roman"))
p <- FeaturePlot(MERGE, 
                 features = c("Col3a1","Col1a1",      # fibroblast
                              "Krt18", "Epcam", "Alb",# epithelial
                              "Pecam1", "Eng",        # endothelial
                              "Ptprc",                # immune
                              "Cd3e", "Cd3d",         # T
                              "Klrb1c","Nkg7",        # NK
                              "Cd68","Cd14","C1qc",   # Mono/macro
                              "H2-Ab1","Cst3",        # cDC
                              "Siglech","Bst2",       # pDC
                              "Csf3r","Clec4d",       # Neutro
                              "Cd79a","Ms4a1",        # B
                              "Jchain","Mzb1"         # Plasma
                 ),
                 pt.size = 0.1,
                 ncol = 5,
                 cols = c("lightgrey", "#FF0000"),
                 raster = FALSE)

# 统一主题设置（标题字体、边框、背景）
p <- p & theme(
  panel.border = element_rect(color = "black", size = 1, fill = NA),
  panel.background = element_blank(),
  strip.background = element_blank(),
  panel.grid = element_blank(),
  strip.text = element_text(family = "Times", face = "bold", color = "black", size = 12)
)

ggsave("sFig5e.mdd.filter2-Featureplot-Anno-总体分群主Marker-附图-有边框.png", p, width = 20, height = 15)


#FindMarkers
object.name <- "mdd.liver"
MERGE.sampling <- subset(MERGE, cells=WhichCells(MERGE, downsample=500, seed = seed.use))#抽样
# MERGE.sampling <- MERGE
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs-sampling500-ANNO.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGssampling500-ANNO.csv"),sep = ",",row.names = T,quote = F)
save(MERGE.markers,file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/total_clustering.markers-ANNO.RData")


##########---------save mdd.liver ----------###########
subset_cells <- MERGE
subset_cells@assays$RNA@data <- as.matrix(0)
subset_cells@assays$RNA@scale.data <- as.matrix(0)
save(subset_cells, file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver_use.RData")

meta.data <- MERGE@meta.data
save(meta.data, file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver_meta.data_use.RData")

######----------plots Fig5e.Pieplot-----------#####
# Fig5e
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

MERGE$cluster <- MERGE$Anno.chi
bb <- table(MERGE$group,MERGE$cluster)%>%as.matrix
bb_rowSum <- rowSums(bb)
join_col <- rep(bb_rowSum,length(table(MERGE$cluster))) #rep 分群数量，需修改
bb <- as.data.frame(bb)
bb <- cbind(bb,join_col)
colnames(bb) <- c("group","variable","freq","sum")
percent <- round((bb$freq / bb$sum) * 100, 2)
bb <- cbind(bb,percent)

dat_plot <- bb

# 颜色定义
celltype_colors <- c("Fibroblast"="#ea5c6f", "Epithelial cell"="#f7905a", 
                     "Endothelial cell"="#e187cb", "T cell"="#fb948d", 
                     "NK cell"="#ebed6f", "Mono_Macro cell"="#b2db87", "cDC"="#7ee7bb",
                     "pDC"="#64cccf",  "Neutrophil cell"= "#a9dce6", 
                     "B cell"="#a48cbe", "Plasma cell"="#e4b7d6")
# 绘制饼图
p <- dat_plot %>%
  ggplot(aes(x = "", y = percent, fill = variable)) +
  geom_bar(stat = "identity", width = 1, color = "black") + 
  coord_polar(theta = "y") +                               
  facet_wrap(~ group, ncol = 2) +                          
  scale_fill_manual(values = celltype_colors) +            
  theme_minimal() +                                         
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA), 
    axis.text = element_blank(),                          
    axis.ticks = element_blank(),                          
    panel.grid = element_blank(),                          # 移除网格
    strip.text = element_text(size = 12, face = "bold")    
  ) +
  labs(fill = NULL, y = NULL, x = NULL,                    
       title = "")  # Percentage of different cell types in different group

ggsave("Fig5e.tissue_celltype_piechart.png", p, width = 9, height = 6)
write.csv(dat_plot,file='Fig5e.Source_data_tissue_celltype_piechart.csv')

######----------sFig5g.细胞比例 Facet-----------#####
library(ggpubr)
library(viridis)
dat <- MERGE
dat$cluster <- dat$Anno.chi %>% as.character()
Idents(dat) <- "cluster"
dat$sample <- dat$sample %>% as.character()

Box_dat_produce <- function(x){
  bb <- table(x$sample,x$cluster)%>%as.matrix
  bb_rowSum <- rowSums(bb)
  join_col <- rep(bb_rowSum,length(table(x$cluster))) 
  bb <- as.data.frame(bb)
  bb <- cbind(bb,join_col)
  colnames(bb) <- c("sample","variable","freq","sum")
  percent <- bb$freq/bb$sum
  bb <- cbind(bb,percent)
  group <- rep(x$group%>%as.character%>%unique,length(rownames(bb))) 
  bb <- cbind(bb,group)
}


data1 <- subset(dat,subset=group=="Control")
data2 <- subset(dat,subset=group=="CRS")


dat <- rbind(Box_dat_produce(data1),
             Box_dat_produce(data2))


dat$group <- factor(dat$group,levels=c("Control","CRS"))

my_comparisons <- list( c("Control","CRS"))

dat_plot <- dat


# Plots
library(ggpubr)
library(viridis)

kv <- data.frame(k=c("B cell", "cDC", "Endothelial cell", "Epithelial cell", "Fibroblast", 
                     "Mono_Macro cell", "Neutrophil cell", 
                     "NK cell", "pDC", "Plasma cell", "T cell"),
                 v=c("B","cDC","Endo","Epi","Fibro","Mono/Mac","Neu","NK","pDC","Plasma","T"))

dat_plot$variable1 <- plyr::mapvalues(dat_plot$variable,from = kv$k,to=kv$v)
dat_plot$cluster <- dat_plot$variable
dat_plot$variable <- dat_plot$variable1 

png("sFig5g.Propertion of cell clusters_wilcox.test-facet.png",width =10,height = 6,units = "in", res = 400)
ggboxplot(dat_plot,
          x="variable",y="percent",
          color = "group", palette = "jama",
          add = "point", outlier.colour=NULL) +
  facet_wrap(~ variable, scales = "free_x",ncol=6) +  # 按cluster分面
  theme(axis.text.x=element_text(angle =0,hjust = 0.4,vjust = 1))+
  # stat_compare_means(method = "wilcox.test")
  stat_compare_means(aes(group = group ), 
                     # label = "p",
                     label = "p.signif",
                     method = "wilcox.test",
                     label.y = c(0.4)
  )+
  scale_colour_manual(values = c("#00A087FF" , "#3C5488FF"))+
  labs(x="",y = "percent")
dev.off()


####---------------Fig5c.d.ACH related genes expression------------######
# # Chrn开头
# chrn_genes <- grep("^Chrn", gene_names, value = TRUE)
# # chrna,chrnb,chrnd,chrng开头
# chrna_genes <- grep("^Chrna", gene_names, value = TRUE)
# chrnb_genes <- grep("^Chrnb", gene_names, value = TRUE)

#####--- Chrn开头(常用)
chrn_genes <-  c("Chrna1","Chrna2","Chrna4","Chrna5", "Chrna7","Chrna9", "Chrna10",
                 "Chrnb1","Chrnb2", "Chrnb3", "Chrnb4",
                 "Chrnd","Chrne","Chrng")
chrm_genes <- c("Chrm1", "Chrm2", "Chrm3", "Chrm4")          

all_genes <- c(chrn_genes, chrm_genes)
# 每页绘制6个基因
n_per_page <- 6
n_pages <- ceiling(length(all_genes) / n_per_page)
# sFig5f
for (i in 1:n_pages) {
  gene_subset <- all_genes[((i - 1) * n_per_page + 1):min(i * n_per_page, length(all_genes))]
  png(paste0("sFig5f.Vln_Chr_Page", i, ".png"), width = 10, height = 8, res = 300, units = "in")
  print(VlnPlot(MERGE, features = gene_subset, pt.size = 0.001, ncol = 2))
  dev.off()
}

png(paste0("sFig5f.Vln_Chr_Page1 in all.png"), width = 24, height = 12, res = 400, units = "in")
VlnPlot(MERGE,                          
        pt.size = 0.0001,ncol=6,                   
        features = all_genes)
dev.off()


#Fig5c
# color_tmp <- c("Fibroblast" = "#F8766D","Epithelial cell" = "#D89000","Endothelial cell" = "#A3A500","T cell" = "#39B600",
#                "NK cell" = "#00BF7D","Mono_Macro cell" = "#00BFC4", "cDC" = "#00B0F6", "pDC" = "#9590FF","Neutrophil cell" = "#E76BF3",
#                "B cell" = "#FF62BC", "Plasma cell" = "#B79F00")

png(paste0("Fig5c.Vln_Chrna9.png"), width = 6, height = 4, res = 300, units = "in")
print(VlnPlot(MERGE, features = "Chrna9", pt.size = 0.001, ncol = 1))
dev.off()
png(paste0("Fig5d.Vln-Chrna9-byGroup.png"), width = 9, height = 4, res = 400, units = "in")
VlnPlot(MERGE,                          
        pt.size = 0.0001,                         
        group.by = "cluster", 
        split.by = "group",
        features =  c("Chrna9"))
dev.off() 

####----------Fig5d,6g,6l,sFig5f DEG calculation----#####
# library(Seurat)
# library(dplyr)
clusters <- levels(MERGE)

all_degs <- lapply(clusters, function(cluster) {
  subset_data <- subset(MERGE, idents = cluster)
  Idents(subset_data) <- "group"
  
  degs <- FindMarkers(
    subset_data,
    ident.1 = "CRS",  
    ident.2 = "Control",  
    min.pct = 0,
    logfc.threshold = 0,
    test.use = "wilcox" 
  )
  
  degs$p_adj_BH <- p.adjust(degs$p_val, method = "BH")
  
  degs$cluster <- cluster
  degs$gene <- rownames(degs)
  
  return(degs)  
})

final_degs <- bind_rows(all_degs)

head(final_degs)

final_degs <- final_degs[,c("cluster","gene", "p_val","avg_log2FC", "pct.1", "pct.2", "p_val_adj",
                            "p_adj_BH")]
colnames(final_degs) <- c("cluster","gene", "p_val","avg_log2FC", "pct.1-Suicide", "pct.2-Control","p_val_adj", 
                          "p_adj_BH")

write.csv(final_degs, "DEGs_results_every_clusters-Stress vs Control.csv", row.names = F)


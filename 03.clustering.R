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
MERGE <- subset(MERGE, subset = nFeature_RNA >= 200 & nFeature_RNA <= 10000 & percent.mt <= 25) 

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

merged_seurat <- MERGE 

#FindMarkers
object.name <- "mdd.liver"
MERGE.sampling <- MERGE
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGs.csv"),sep = ",",row.names = T,quote = F)
escc.makers <- MERGE.markers
save(escc.makers,file = "~/Workspace/06_mdd.liver/paper/clustering-total/total_clustering.markers.RData")


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


Idents(MERGE) <- 'group'
p <- DimPlot(MERGE,pt.size = 0.05,label=F,label.size=4,cols = color_list,raster=F)
ggsave("MDD.liver-UMAP-byGroup.png", p, width = 5, height = 4)
Idents(MERGE) <- 'sample'
p <- DimPlot(MERGE,pt.size = 0.05,label=F,label.size=4,cols = color_sample[2:10],raster=F)
ggsave("MDD.liver-UMAP-bySample.png", p, width = 5, height = 4)

# Extended Fig5a
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
ggsave("mdd.filter2-Dotplot-Anno-global cluster Marker.png", p, width = 10, height = 4)
ggsave("mdd.filter2-Dotplot-Anno-global cluster Marker.pdf", p, width = 10, height = 4)


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

# set theme
p <- p & theme(
  panel.border = element_rect(color = "black", size = 1, fill = NA),
  panel.background = element_blank(),
  strip.background = element_blank(),
  panel.grid = element_blank(),
  strip.text = element_text(family = "Times", face = "bold", color = "black", size = 12)
)

ggsave("mdd.filter2-Featureplot-Anno-global cluster Marker.png", p, width = 20, height = 15)


#FindMarkers
object.name <- "mdd.liver"
MERGE.sampling <- subset(MERGE, cells=WhichCells(MERGE, downsample=500, seed = seed.use))
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGssampling500-ANNO.csv"),sep = ",",row.names = T,quote = F)
save(MERGE.markers,file = "~/06_mdd.liver/paper/clustering-total/total_clustering.markers-ANNO.RData")


##########---------save mdd.liver ----------###########
subset_cells <- MERGE
save(subset_cells, file = "~/Workspace/06_mdd.liver/paper/clustering-total/mdd.liver_use.RData")


## Supplementary Figure 11b
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

# color set
celltype_colors <- c("Fibroblast"="#ea5c6f", "Epithelial cell"="#f7905a", 
                     "Endothelial cell"="#e187cb", "T cell"="#fb948d", 
                     "NK cell"="#ebed6f", "Mono_Macro cell"="#b2db87", "cDC"="#7ee7bb",
                     "pDC"="#64cccf",  "Neutrophil cell"= "#a9dce6", 
                     "B cell"="#a48cbe", "Plasma cell"="#e4b7d6")
# pieplot
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

ggsave("Supple Fig11b.tissue_celltype_piechart.png", p, width = 9, height = 6)
write.csv(dat_plot,file='upple Fig11b.Source_data_tissue_celltype_piechart.csv')

######----------cell percent Facet-----------#####
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


data1 <- subset(dat,subset=group=="HCC")
data2 <- subset(dat,subset=group=="HCC+CRS")


dat <- rbind(Box_dat_produce(data1),
             Box_dat_produce(data2))


dat$group <- factor(dat$group,levels=c("HCC","HCC+CRS"))

my_comparisons <- list( c("HCC","HCC+CRS"))

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

png("Supplementary Figure 12a-Propertion of cell clusters_wilcox.test-facet.png",width =10,height = 6,units = "in", res = 400)
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


####--------------- ACH related genes expression------------######
# # grep Chrn,Chrm
chrn_genes <- grep("^Chrn", gene_names, value = TRUE)
chrn_genes <- grep("^Chrm", gene_names, value = TRUE)
all_genes <- c(chrn_genes, chrm_genes)

# 6 genes for every page
n_per_page <- 6
n_pages <- ceiling(length(all_genes) / n_per_page)

for (i in 1:n_pages) {
  gene_subset <- all_genes[((i - 1) * n_per_page + 1):min(i * n_per_page, length(all_genes))]
  png(paste0("Vln_Chr_Page", i, ".png"), width = 10, height = 8, res = 300, units = "in")
  print(VlnPlot(MERGE, features = gene_subset, pt.size = 0.001, ncol = 2))
  dev.off()
}

png(paste0("Vln_Chr_Page1 in all.png"), width = 24, height = 12, res = 400, units = "in")
VlnPlot(MERGE,                          
        pt.size = 0.0001,ncol=6,                   
        features = all_genes)
dev.off()


#Fig5c
png(paste0("Fig5b.Vln-Chrna9-byGroup.png"), width = 9, height = 4, res = 400, units = "in")
VlnPlot(MERGE,                          
        pt.size = 0.0001,                         
        group.by = "cluster", 
        split.by = "group",
        features =  c("Chrna9"))
dev.off() 


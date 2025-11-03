# ---------------------------------------------
# Seurat harmony and clustering and Annotation
# --------------------------------------------
setwd("/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total")
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
MERGE.sampling <- subset(MERGE, cells=WhichCells(escc, downsample=500, seed = seed.use))#抽样
# MERGE.sampling <- MERGE
MERGE.markers <- FindMarkers_parallel(MERGE.sampling, mc.cores = 5)
MERGE.markers %>% TOP_N(50, pct.1 = 0.2) -> top50
MERGE.markers <- MERGE.markers %>% TOP_N(5000)
write.table(top50,file = paste0(object.name,"_res", res.use ,"_dim", dim.use ,"_culster_top50_DEGs-sampling500.csv"),sep = ",", row.names = T, quote = F)
write.table(MERGE.markers,file = paste0(object.name,"_res",res.use ,"_dim", dim.use ,"_culster_all_DEGssampling500.csv"),sep = ",",row.names = T,quote = F)
escc.makers <- MERGE.markers
save(escc.makers,file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/total_clustering.markers-sampling500.RData")


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
save(subset_cells, file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver-v0/mdd.liver.RData")

meta.data <- MERGE@meta.data
save(meta.data, file = "/public/home/chidm/Workspace/ouyang/06_mdd.liver/paper/clustering-total/mdd.liver-v0/mdd.liver_meta.data.RData")

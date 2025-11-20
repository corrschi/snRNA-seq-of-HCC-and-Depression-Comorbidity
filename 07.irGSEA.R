setwd("~/Workspace/ouyang/06_mdd.liver/paper/clustering-epi/irGSEA")

library(Seurat)
library(irGSEA)
library(tidyverse)
# load data
load("./epi.RData")
epi@active.ident%>%table

pbmc3k.final <- epi
pbmc3k.final$Anno.epi <- paste0("Epi.c", pbmc3k.final$seurat_clusters)
pbmc3k.final$Anno.epi%>%table

Idents(pbmc3k.final) <- 'Anno.epi'
pbmc3k.final @active.ident%>%table

p <- DimPlot(pbmc3k.final, reduction = "umap",label = T,label.size=6) + NoLegend()
ggsave("1.seurat-umap-check.png", plot = p, width = 6, height = 6, units = "in", dpi = 400)
pbmc3k.final$seurat_annotations <- pbmc3k.final$Anno.epi #对应软件里的分群命名


### irGSEA 
# Calculate enrichment scores 
pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", 
                             slot = "data", 
                             seeds = 123, 
                             ncores = 4,
                             min.cells = 3,
                             min.feature = 0,
                             custom = F, geneset = NULL, 
                             msigdb = T, 
                             species = "Mus musculus", 
                             category = "H",            
                             subcategory = NULL,      
                             geneid = "symbol",     
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea", "JASMINE", "viper"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')    

Seurat::Assays(pbmc3k.final)
save(pbmc3k.final,file="pbmc3k.final-Hallmark.RData")
# Integrate differential gene set
result.dge <- irGSEA.integrate(object = pbmc3k.final,
                               group.by = "seurat_annotations",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

class(result.dge)
#> [1] "list"

save(result.dge,file = "result.dge-HallMark.RData")

#4 可视化
load("./result.dge-KEGG.RData")
# load("./result.dge-HallMark.RData")
# result.dge %>% str

pdf(paste0("sFig6a.epi.filter1-irGSEA.heatmap.plot-KEGG.pdf"),width = 9, height = 6)
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
print(irGSEA.heatmap.plot)
dev.off()


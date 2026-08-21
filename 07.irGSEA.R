setwd("~/Workspace/06_mdd.liver/paper/clustering-epi/irGSEA")

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
pbmc3k.final$seurat_annotations <- pbmc3k.final$Anno.epi 


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

# ridgeplot
ridgeplot <- irGSEA.ridgeplot(
  object = pbmc3k.final,
  method = "singscore",
  show.geneset = "KEGG-NEUROACTIVE-LIGAND-RECEPTOR-INTERACTION"
)

# 添加 Times New Roman 字体样式
ridgeplot <- ridgeplot + theme(
  plot.title = element_text(family = "Times New Roman", face = "bold", size = 14),
  axis.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  axis.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.text = element_text(family = "Times New Roman", face = "bold", size = 10),
  legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
  strip.text = element_text(family = "Times New Roman", face = "bold", size = 12))
png("Fig5e.irGSEA.ridgeplot.png", 
    width = 7, height = 7, units = "in", res = 400)
print(ridgeplot)
dev.off()


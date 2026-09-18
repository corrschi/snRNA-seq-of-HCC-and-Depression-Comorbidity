###----------basic loading----#####
library(harmony)
library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(parallel)
library(stringi)
library(plyr)
library(ggthemes)
library(cowplot)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(viridis)
library(reshape2)
library(scales)
library(rlang)
library(dendextend)
library(DoubletFinder)
library(future)
library(ggsci)
library(ggpubr)
library(Matrix)
library(irlba)
library(Nebulosa)
library(grDevices)

options(future.globals.maxSize = 10*1000 * 1024^2)
plan("multicore", workers = 8) 
options(future.globals.maxSize = 10000000000000)

load("/public/home/chidm/Documents/Function_inhouse/Seurat_fun_inhouse.RData")
load("/public/home/chidm/Documents/Function_inhouse/TOP_N_new.RData")
load("/public/home/chidm/Documents/Function_inhouse/FindMarkers_parallel.RData")

###----------color set----#####
library(ggsci)
color_all <- c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f",
               "#b2db87","#7ee7bb","#64cccf","#a9dce6","#a48cbe","#e4b7d6")


color_epi <- c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", 
               "#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", 
               "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", 
               "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35", 
               "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD", 
               "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175")

color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
color_macro <- color_used
color_sample <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
                  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
                  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", 
                  "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", 
                  "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                  "#E31A1C")
sessionInfo()


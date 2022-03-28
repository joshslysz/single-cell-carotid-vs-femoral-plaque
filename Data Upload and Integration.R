## Load libraries ##
library(Seurat)
library(DietSeurat)
library(tidyverse)
library(ggpubr)
library(limma)
library(SingleR)
library(enrichR)
library(cowplot)
library(patchwork)
library(VennDiagram)
library(pheatmap)
library(RColorBrewer)
library(monocle3)
library(msigdbr)
library(scmap)
library(scater)
library(msigdb)
library(SCINA)
library(devtools)
library(cerebroApp)
library(shiny)
library(SeuratWrappers)
library(UpSetR)
library(harmony)
library(muscat)
library(glmmTMB)


## Get the data ##
## Femoral data ##
#FEMPIC01
femPIC01.data <- Read10X(data.dir = "PIC01_Carotid_OUTS")
srtfemPIC01<- CreateSeuratObject(counts = femPIC01.data, project = "femoral")
#FEMPIC01
femPIC42.data <- Read10X(data.dir = "PIC42_Carotid_OUTS")
srtfemPIC42<- CreateSeuratObject(counts = femPIC42.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_srt <- merge(x = srtfemPIC01, y = srtfemPIC42)
#FEMPIC45
femPIC45.data <- Read10X(data.dir = "PIC45_Carotid_OUTS")
srtfemPIC45<- CreateSeuratObject(counts = femPIC45.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_srtfemPIC45_srt <- merge(x = femPIC01_femPIC42_srt, y = srtfemPIC45)
#FEMPIC46
femPIC46.data <- Read10X(data.dir = "PIC46_Carotid_OUTS")
srtfemPIC46<- CreateSeuratObject(counts = femPIC46.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_srtfemPIC45_srtfemPIC46_srt <- merge(x = femPIC01_femPIC42_srtfemPIC45_srt, y = srtfemPIC46)
#FEMPIC48
femPIC48.data <- Read10X(data.dir = "PIC48_Carotid_OUTS")
srtfemPIC48<- CreateSeuratObject(counts = femPIC48.data, project = "femoral")
#COMBINE
full_femoral_srt <- merge(x = femPIC01_femPIC42_srtfemPIC45_srtfemPIC46_srt, y = srtfemPIC48)

## carotid data ##
car1.data <- Read10X(data.dir = "PIC01_Carotid_OUTS")
srtcarPIC01<- CreateSeuratObject(counts = carPIC01.data, project = "carotid")


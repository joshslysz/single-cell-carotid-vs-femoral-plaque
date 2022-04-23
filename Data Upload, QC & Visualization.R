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

################# DATA UPLOAD AND QC ON SEPARATE SAMPLES ######################
#### FEMORAL SAMPLES ####
### Upload Femoral data ###
#FEMPIC01
femPIC01.data <- Read10X(data.dir = "PIC01_Femoral_OUTS")
srtfemPIC01<- CreateSeuratObject(counts = femPIC01.data, project = "femoral")
#FEMPIC01
femPIC42.data <- Read10X(data.dir = "PIC42_Femoral_OUTS")
srtfemPIC42<- CreateSeuratObject(counts = femPIC42.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_srt <- merge(x = srtfemPIC01, y = srtfemPIC42)
#FEMPIC45
femPIC45.data <- Read10X(data.dir = "PIC45_Femoral_OUTS")
srtfemPIC45<- CreateSeuratObject(counts = femPIC45.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_femPIC45_srt <- merge(x = femPIC01_femPIC42_srt, y = srtfemPIC45)
#FEMPIC46
femPIC46.data <- Read10X(data.dir = "PIC46_Femoral_OUTS")
srtfemPIC46<- CreateSeuratObject(counts = femPIC46.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_femPIC45_femPIC46_srt <- merge(x = femPIC01_femPIC42_femPIC45_srt, y = srtfemPIC46)
#FEMPIC48
femPIC48.data <- Read10X(data.dir = "PIC48_Femoral_OUTS")
srtfemPIC48<- CreateSeuratObject(counts = femPIC48.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_femPIC45_femPIC46_femPIC48_srt <- merge(x = femPIC01_femPIC42_femPIC45_femPIC46_srt, y = srtfemPIC48)
#FEMPIC50
femPIC50.data <- Read10X(data.dir = "PIC50_Femoral_OUTS")
srtfemPIC50<- CreateSeuratObject(counts = femPIC50.data, project = "femoral")
#COMBINE
femPIC01_femPIC42_femPIC45_femPIC46_femPIC48_femPIC50_srt <- merge(x = femPIC01_femPIC42_femPIC45_femPIC46_femPIC48_srt, y = srtfemPIC50)
#FEMPIC53
femPIC53.data <- Read10X(data.dir = "PIC53_Femoral_OUTS")
srtfemPIC53<- CreateSeuratObject(counts = femPIC53.data, project = "femoral")
#COMBINE
full_femoral_srt <- merge(x = femPIC01_femPIC42_femPIC45_femPIC46_femPIC48_femPIC50_srt, y = srtfemPIC53)

### Quality Control Femoral data ###
# Add in the Mitochondrial PCT% information
full_femoral_srt$percent.mt <- PercentageFeatureSet(full_femoral_srt, pattern = "^MT-")

# nCount_RNA is the number of UMI counts in a cell
hist(full_femoral_srt$nCount_RNA)

# nFeature_RNA is the number of different genes that had any reads
hist(full_femoral_srt$nFeature_RNA)

# percent.mt is the percent mitochondrial reads
hist(full_femoral_srt$percent.mt)

# Make a violin plot of the QC columns
VlnQCplt <- VlnPlot(full_femoral_srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3)
VlnQCplt

ggsave(filename = "QC_Vln_Femoral.png", plot = plt, width = 7, height = 3.5)
ggsave(filename = "QC_Vln_Femoral_aftersubet.png", plot = VlnQCplt, width = 7, height = 3.5)

# Filter out unwanted cells from full femoral srt
full_femoral_srt <- subset(full_femoral_srt, subset = 
                             nFeature_RNA > 200 & nFeature_RNA < 2500 
                           & nCount_RNA > 200 & nCount_RNA < 10000 
                           & percent.mt < 10) 

#How many cells were lost? #should now be 36063 cells
full_femoral_srt

#### CAROTID SAMPLES ####
### Upload Carotid data ###
#CARPIC01
carPIC01.data <- Read10X(data.dir = "PIC01_Carotid_OUTS")
srtcarPIC01<- CreateSeuratObject(counts = carPIC01.data, project = "carotid")
#CARPIC49
carPIC49.data <- Read10X(data.dir = "PIC49_Carotid_OUTS")
srtcarPIC49<- CreateSeuratObject(counts = carPIC49.data, project = "carotid")
#COMBINE
carPIC01_carPIC49_srt <- merge(x = srtcarPIC01, y = srtcarPIC49)
#CARPIC51
carPIC51.data <- Read10X(data.dir = "PIC51_Carotid_OUTS")
srtcarPIC51<- CreateSeuratObject(counts = carPIC51.data, project = "carotid")
#COMBINE
full_carotid_srt <- merge(x = carPIC01_carPIC49_srt, y = srtcarPIC51)

###Quality Control Carotid data ###
# Add in the Mitochondrial PCT% information
full_carotid_srt$percent.mt <- PercentageFeatureSet(full_carotid_srt, pattern = "^MT-")

# nCount_RNA is the number of UMI counts in a cell
hist(full_carotid_srt$nCount_RNA)

# nFeature_RNA is the number of different genes that had any reads
hist(full_carotid_srt$nFeature_RNA)

# percent.mt is the percent mitochondrial reads
hist(full_carotid_srt$percent.mt)

# Make a violin plot of the QC columns
VLnplt_Carotid <- VlnPlot(full_carotid_srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3)
VLnplt_Carotid 
ggsave(filename = "QC_Full_Carotid_after.png", plot = VLnplt_Carotid, width = 7, height = 3.5)

# Filter out unwanted cells from full carotid srt
full_carotid_srt <- subset(full_carotid_srt, subset = 
                             nFeature_RNA > 200 & nFeature_RNA < 2500 
                           & nCount_RNA > 200 & nCount_RNA < 10000 
                           & percent.mt < 10) 
#How many cells were lost? #should now be 36063 cells
full_carotid_srt

####################### END DATA UPLOAD AND QUALITY CONTROL ON SEPARATE SAMPLES #######################

####################### COMBINE SAMPLES, NORMALIZE, SCALE ######################
### combine full carotid srt and full femoral srt ###
Carotid_Femoral_srt <- merge(x= full_femoral_srt, y= full_carotid_srt)
### Normalize and Scale the data ###
# Log-transform the counts
Carotid_Femoral_srt <- NormalizeData(Carotid_Femoral_srt)
# Find Variable Features
Carotid_Femoral_srt <- FindVariableFeatures(Carotid_Femoral_srt)
# Scale the data
Carotid_Femoral_srt <- ScaleData(Carotid_Femoral_srt)

####################### Visualize Carotid vs femoral with PCA BEFORE Integration and correcting for batch ######################
# Run PCA
Carotid_Femoral_srt <- RunPCA(Carotid_Femoral_srt)
# Plot the PCA
DimPlot(Carotid_Femoral_srt, reduction = "pca", group.by = "orig.ident")
# Plot the loadings
VizDimLoadings(Carotid_Femoral_srt, dims = 1:2, reduction = "pca")

####################### Visualize Carotid vs femoral with UMAP BEFORE Integration and correcting for batch ######################
# Choose the number of principle components to keep
ElbowPlot(Carotid_Femoral_srt,ndims = 50)
# Find nearest neighbors and construct the graph
Carotid_Femoral_srt <- FindNeighbors(Carotid_Femoral_srt, k.param = 50, dims = 1:40)
# Find the clusters
Carotid_Femoral_srt <- FindClusters(Carotid_Femoral_srt, resolution = 0.35)
# Get the UMAP embedding
Carotid_Femoral_srt <- RunUMAP(Carotid_Femoral_srt, dims = 1:40)
# Plot the UMAP with clustering
DimPlot(Carotid_Femoral_srt, reduction = "umap", label = TRUE)
# Dim Plots
DimPlot(Carotid_Femoral_srt, reduction = "umap", group.by = "orig.ident") 
DimPlot(Carotid_Femoral_srt, group.by = "orig.ident") + DimPlot(Carotid_Femoral_srt, group.by = "seurat_clusters")




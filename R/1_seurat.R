library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(clustree)

# Set here your working directory
# setwd("~/Desktop/eae_microglia/")
###Read data and generate seurat objects#####
#load data
load(file.path("data", "SC_NT2.rda"))

### First check
SC_NT2[["percent.mt"]] <- PercentageFeatureSet(SC_NT2, pattern = "^mt-")
SC_NT2[["percent.Rps"]] <- PercentageFeatureSet(SC_NT2, pattern = "^Rps")
SC_NT2[["percent.Rpl"]] <- PercentageFeatureSet(SC_NT2, pattern = "^Rpl")
PlotQC1_Raw<-VlnPlot(SC_NT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rps", "percent.Rpl"), ncol = 5)
ggsave( file.path("plots", "QC", "Qc1.pdf"),PlotQC1_Raw, device="pdf", width = 15)

### Filter according to RNA, mt, Rps, Rpl 
SC_NT2 <- subset(SC_NT2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
PlotQC2 <-VlnPlot(SC_NT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rps", "percent.Rpl"), ncol = 5)
ggsave( file.path("plots", "QC", "Qc2_Filtered.pdf"),PlotQC2, device="pdf", width = 15)

### create UMAP
SC_NT2 <- NormalizeData(SC_NT2, normalization.method = "LogNormalize")
SC_NT2 <- FindVariableFeatures(SC_NT2, selection.method = "vst", nfeatures = 10000)
SC_NT2 <- ScaleData(SC_NT2)
SC_NT2 <- RunPCA(SC_NT2)
ElbowPlot(SC_NT2)
SC_NT2 <- RunUMAP(SC_NT2, dims = 1:15)
DimPlot(SC_NT2, reduction = "umap")
ggsave( file.path("plots", "QC", "UMAP_WithDoublets.pdf"), width = 8, height=8)

## Remove doublets
sce <- SingleCellExperiment(assays=list(counts=SC_NT2@assays$RNA@counts))
sce <- scDblFinder(sce)
SC_NT2s<-SC_NT2[,sce@colData$scDblFinder.class == "singlet"]
SC_NT2s[["percent.mt"]] <- PercentageFeatureSet(SC_NT2s, pattern = "^mt-")
SC_NT2s[["percent.Rps"]] <- PercentageFeatureSet(SC_NT2s, pattern = "^Rps")
SC_NT2s[["percent.Rpl"]] <- PercentageFeatureSet(SC_NT2s, pattern = "^Rpl")
PlotQC3_NoDbl <-VlnPlot(SC_NT2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rps", "percent.Rpl"), ncol = 5)
ggsave( file.path("plots", "QC", "Qc3_NoDoublets.pdf"),PlotQC3_NoDbl, device="pdf", width = 15)
Umap_NoDbl <- DimPlot(SC_NT2s, reduction = "umap")
ggsave( file.path("plots", "QC", "UMAP_NoDoublets.pdf"),Umap_NoDbl, device="pdf", width = 8, height=8)




# clustering 
SC_NT2s <- FindNeighbors(SC_NT2s, dims = 1:15, verbose = FALSE)
SC_NT2s <- FindClusters(SC_NT2s, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(SC_NT2s)

SC_NT2s <- FindClusters(SC_NT2s, resolution = 1.6)

DimPlot(SC_NT2s, label = TRUE) + NoLegend()
ggsave( file.path("plots", "umap", "UMAP_Clusters.pdf"),Umap_NoDbl, device="pdf", width = 6, height=6)


#save data
save(SC_NT2s, file = file.path("data", "seurat_NOdoublets_clusters.RData"))


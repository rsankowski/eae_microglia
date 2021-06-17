library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(ggplot2)
library(clustree)

# Set here your working directory
# setwd("~/Desktop/eae_microglia/")
###Read data and generate seurat objects#####
#load data
load(file.path("data", "Microglia_Fillatreau", "SC_NT2.rda"))

SC_NT2[["percent.mt"]] <- PercentageFeatureSet(SC_NT2, pattern = "^mt-")
PlotBefore<-VlnPlot(SC_NT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave( file.path("plots", "QC", "Qc_with_doublets.pdf"),PlotBefore, device="pdf", width = 12)
SC_NT2 <- NormalizeData(SC_NT2, normalization.method = "LogNormalize")
SC_NT2 <- FindVariableFeatures(SC_NT2, selection.method = "vst", nfeatures = 2000)
SC_NT2 <- ScaleData(SC_NT2)
SC_NT2 <- RunPCA(SC_NT2)
ElbowPlot(SC_NT2)
SC_NT2 <- RunUMAP(SC_NT2, dims = 1:15)
UmapBefore <- DimPlot(SC_NT2, reduction = "umap")
ggsave( file.path("plots", "QC", "UMAP_with_doublets.pdf"),UmapBefore, device="pdf", width = 6, height=6)

## Remove doublets
sce <- SingleCellExperiment(assays=list(counts=SC_NT2@assays$RNA@counts))
sce <- scDblFinder(sce)
SC_NT2s<-SC_NT2[,sce@colData$scDblFinder.class == "singlet"]
SC_NT2s[["percent.mt"]] <- PercentageFeatureSet(SC_NT2s, pattern = "^mt-")
PlotAfter<-VlnPlot(SC_NT2s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave( file.path("plots", "QC", "Qc_without_doublets.pdf"),PlotAfter, device="pdf", width = 12)
UmapAfter <- DimPlot(SC_NT2s, reduction = "umap")
ggsave( file.path("plots", "QC", "UMAP_without_doublets.pdf"),UmapAfter, device="pdf", width = 6, height=6)




# clustering 
SC_NT2s <- FindNeighbors(SC_NT2s, dims = 1:15, verbose = FALSE)
SC_NT2s <- FindClusters(SC_NT2s, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(SC_NT2s)

SC_NT2s <- FindClusters(SC_NT2s, resolution = 0.6)

Umap <- DimPlot(SC_NT2s, label = TRUE) + NoLegend()
ggsave( file.path("plots", "umap", "UMAP_Clusters.pdf"),UmapAfter, device="pdf", width = 6, height=6)


#save data
save(SC_NT2s, file = file.path("data", "seurat_NOdoublets_clusters.RData"))

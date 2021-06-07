library(Seurat)
library(cowplot)
library(patchwork)
library(harmony)
library(tidyverse)
library(viridis)
library(Matrix)
library(clustree)
library(assertthat)
library(ggtree)
library(DoubletFinder)

source(file.path("R", "functions.R"))

date = Sys.Date()

###Read data and generate seurat objects#####
#load data
load(file.path("data", "Microglia_Fillatreau", "SC_NT2.rda"))

#regress out mitochondrial
SC_NT2 <- PercentageFeatureSet(SC_NT2, pattern = "^mt-", col.name = "percent.mt")

# run sctransform
SC_NT2 <- SCTransform(SC_NT2, vars.to.regress = "percent.mt", verbose = FALSE)

#pca etc
SC_NT2 <- RunPCA(SC_NT2, verbose = FALSE)

#run elbow plot to find most relevant PCs
ElbowPlot(SC_NT2)

#run UMAP and clustering on chosen PCs
SC_NT2 <- RunUMAP(SC_NT2, dims = 1:30, verbose = FALSE)
SC_NT2 <- FindNeighbors(SC_NT2, dims = 1:30, verbose = FALSE)
SC_NT2 <- FindClusters(SC_NT2, resolution = seq(from=0.2,to=2,by=0.2))

#run clustree to find a cluster resolution
#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(SC_NT2)

SC_NT2 <- FindClusters(SC_NT2, resolution = 1.4)

DimPlot(SC_NT2, label = TRUE) + NoLegend()

#save data
save(SC_NT2, file = file.path("data", "seurat_with_doublets.RData"))

library(Seurat)
library(SingleR)
library(RaceID)
library(scuttle)
library(tidyverse)
library(SingleCellExperiment)

# Set working directory
# setwd("~/Desktop/eae_microglia/")

### load Data to classify
load(file.path("data", "seurat_NOdoublets_clusters.RData"))
### load Reference dataset
load(file.path("data", "RaceID_data_pv.RData"))
load(file.path("data", "metadata_pv.RData"))

#extract variable features 
varfeat <- VariableFeatures(SC_NT2s)
varfeat <- varfeat[varfeat %in% rownames(sc@expdata)]

PeakData <- SingleCellExperiment(assays=list(counts=SC_NT2s@assays$RNA@counts[varfeat,]))
PeakDataNorm <- logNormCounts(PeakData)  

#reference dataset
Ref <- SingleCellExperiment(assays=list(counts=sc@expdata[varfeat, metadata$ID] ))
RefNorm <- logNormCounts(Ref)     


PredMicroglia <- SingleR(test=PeakDataNorm, ref=RefNorm, labels=metadata$Subpopulation,de.method="wilcox")

SC_NT2s@meta.data$LabelsSingleR <- PredMicroglia$labels
DimPlot(SC_NT2s, reduction = "umap",group.by = "LabelsSingleR",label = T)
ggsave( file.path("plots", "umap", "UMAP_LabelsSingleR.pdf"),device="pdf", width = 10, height=10)

#save data
save(SC_NT2s, file = file.path("data", "seurat_NOdoublets_clusters.RData"))

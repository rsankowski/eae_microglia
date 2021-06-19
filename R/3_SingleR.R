library(Seurat)
library(SingleR)
library(RaceID)
library(scuttle)
library(ggplot2)

# Set working directory
# setwd("~/Desktop/eae_microglia/")

### Data to classify
load(file.path("data", "seurat_NOdoublets_clusters.RData"))
PeakData <- SingleCellExperiment(assays=list(counts=SC_NT2s@assays$RNA@counts))
PeakDataNorm <- logNormCounts(PeakData)  



### Reference dataset
load(file.path("data", "Microglia_Roman", "RaceID_data_pv.RData"))
load(file.path("data", "Microglia_Roman", "metadata_pv.RData"))
Ref <- SingleCellExperiment(assays=list(counts=sc@expdata[, metadata$ID] ))
RefNorm <- logNormCounts(Ref)     


PredMicroglia <- SingleR(test=PeakDataNorm, ref=RefNorm, labels=metadata$Subpopulation,de.method="wilcox")

SC_NT2s@meta.data$LabelsSingleR <- PredMicroglia$labels
Umap_Labels <- DimPlot(SC_NT2s, reduction = "umap",group.by = "LabelsSingleR",label = T)
ggsave( file.path("plots", "umap", "UMAP_LabelsSingleR.pdf"),Umap_Labels, device="pdf", width = 10, height=10)

#save data
save(SC_NT2s, file = file.path("data", "seurat_NOdoublets_clusters.RData"))

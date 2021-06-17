library(Seurat)
library(tidyverse)
library(readxl)

# Set working directory
# setwd("~/Desktop/eae_microglia/")
load(file.path("data", "seurat_NOdoublets_clusters.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= SC_NT2s@meta.data[,"seurat_clusters"], row.names = rownames(SC_NT2s@meta.data)) %>%
  bind_cols(as.data.frame(t(SC_NT2s[["RNA"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()
  
rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(SC_NT2s) <- order_clusters

#extract metadata
metadata <- SC_NT2s@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(SC_NT2s))

# cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), object = SC_NT2s, point_size = .5, .retain_cl = levels(SC_NT2s)) + labs(subtitle=paste0(as.character(i), ' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap",  paste0(i,"_signature_all.pdf")), useDingbats=F)
} 

library(Seurat)
library(tidyverse)
library(readxl)

# Set working directory
# setwd("~/Desktop/eae_microglia/")
load(file.path("data", "seurat_NOdoublets_clusters.RData"))
source(file.path("R", "functions.R"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= SC_NT2s@meta.data[,"seurat_clusters"], row.names = rownames(SC_NT2s@meta.data)) %>%
  bind_cols(as.data.frame(t(SC_NT2s[["RNA"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()
  
rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(SC_NT2s) <- rev(order_clusters)

#reorder levels of the SingleR labels
SC_NT2s@meta.data$LabelsSingleR <- factor(SC_NT2s$LabelsSingleR, levels = c("hpvMF", "dapvMF1", "pMC1", "pMC2", "pMC3", "pMC4", "pMC5", "hMG1", "daMG1", "daMG2", "daMG3", "daMG4", "Lympho", "Granulocytes"))

#extract metadata
metadata <- SC_NT2s@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(SC_NT2s))
metadata$LabelsSingleR <- factor(metadata$LabelsSingleR, levels = c("hpvMF", "dapvMF1", "pMC1", "pMC2", "pMC3", "pMC4", "pMC5", "hMG1", "daMG1", "daMG2", "daMG3", "daMG4", "Lympho", "Granulocytes"))

# cell signature umaps
signature_genes <-  read_excel(file.path("data","cell_signatures.xlsx"), "Core signature", skip = 2)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), object = SC_NT2s, point_size = .5, .retain_cl = levels(SC_NT2s)) + labs(subtitle=paste0(as.character(i), ' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap",  paste0(i,"_signature_all.pdf")), useDingbats=F)
} 

#plot singleR classifications
#umap
DimPlot(SC_NT2s, group.by = "LabelsSingleR", label = T) +
  scale_color_manual(values=colors_many)
ggsave(file.path("plots", "umap",  "singleR_labels.pdf"), useDingbats=F)

#marimekko
mosaicGG2(metadata, "seurat_clusters", "LabelsSingleR")
ggsave(file.path("plots", "others",  "singleR_labels_marimekko.pdf"), useDingbats=F)

#differential gene expression in subsets
SC_NT2s_2 <- SC_NT2s
Idents(SC_NT2s_2) <- SC_NT2s$LabelsSingleR

#find_markers for the cell Subsets
if (!file.exists(file.path("data", "markers_cell_subsets.RData"))) {
  all.markers <- FindAllMarkers(SC_NT2s_2) 

  save(all.markers, file = file.path("data", "markers_cell_subsets.RData"))
  write.csv(all.markers, file.path("data", "markers_cell_subsets.csv"))
  
} else {
  load(file.path("data", "markers_cell_subsets.RData"))
}

top10 <- all.markers %>% 
  #remove non annotated genes
  filter(!gene %in% grep("(^RP|mt-|Gm)", .$gene, value = T)) %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(SC_NT2s_2,features = top10$gene, group.colors = colors_many) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(file.path("plots", "heatmaps", "top10_subset_markers_heatmap.pdf"), useDingbats=F)

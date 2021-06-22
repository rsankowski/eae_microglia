
# Set here your working directory
# setwd("~/Desktop/eae_microglia/")

#setup the folders
dir.create("data")
dir.create(file.path("data", "Microglia_Fillatreau"))
dir.create("R")
dir.create("plots")
dir.create(file.path("plots", "QC"))
dir.create(file.path("plots", "umap"))
# dir.create(file.path("plots", "heatmaps"))
# dir.create(file.path("plots", "others"))

# the data object SC_NT2.rda is in the Microglia_Fillatreau folder

#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#the following code was adjusted from url https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c(
  "SingleCellExperiment",
  "SingleR",
  "Seurat",
  "RColorBrewer",
  "tidyverse",
  "assertthat",
  "viridis",
  "readxl",
  "scDblFinder",
  "clustree",
  "scuttle",
  "scran"
  
)


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)


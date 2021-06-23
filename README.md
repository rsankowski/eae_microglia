# immune cells from Spinal Cord
Analysis of 10x immune cells from Spinal Cord


# Instruction
This repository contains 4 numbered scripts that should be run sequentially to perform the analysis. The purpose of each script is briefly described here:
0_setup.R  --  it creates the folders which contains data and results and it make sure you have all the required R packages installed. After running this script, you should put the single cell objects in the folder "data"
1_seurat.R  --  quality control, doublet removal and clustering is performed on the Seurat object. A new Seruat object is created and saved
2_SingleR.R  --  classification of the cells is perfomed with the SingleR package
3_plotting.R  --  UMAP plots that help the identification of the various cell types using known markers and the classification done with SingleR

In addition to these scripts, the "functions.R" script is included in the repository. This file contains supplementary functions.



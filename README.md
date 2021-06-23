# immune cells from Spinal Cord
Analysis of 10x single cell data of CD45+ from Spinal Cord


# Instruction
This repository contains 4 numbered scripts that should be run sequentially to perform the analysis. The purpose of each script is briefly described here:

* 0_setup.R  --  it creates the folders which will be populated with data and results and it make sure you have all the required R packages installed. After running this script, you should put the single cell objects (SC_NT2.rda, RaceID_data_pv.RData, metadata_pv.RData) in the folder "data".
* 1_seurat.R  --  quality control, doublet removal and clustering is performed on the Seurat object. The resolution for clustering is selected using the clustree package. A new Seurat object is created and saved.
* 2_SingleR.R  --  classification of the cells is performed with the SingleR package.
* 3_plotting.R  --  exploratory plots that help the identification of the various cell types using known markers and the classification done with SingleR. Also cluster and cell subset markers are being identified and plotted.

In addition to these scripts, the "functions.R" script is included in the repository. This file contains supplementary functions.



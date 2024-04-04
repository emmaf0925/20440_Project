# 20440_Project
### Evaluating the Correlation Between Stress and Age in the Mouse Brain at the Single Cell Transcriptomic Level

Emma Finburgh and Nickeisha Cuthbert

## Overview
This repository currently contains the code required to perform the processing, integration, dimensional reduction, and UMAP visualization of two single cell RNA sequencing datasets (old/young and stressed/baseline mice) of the mouse brain. The goal of the analysis is to determine what transcriptomic signatures are shared between aged and stressed mice, as this correlation has been alluded to in literature but not quanlified at the single cell transcriptomic level. This initial analysis involves filtering, processing, and normalizing the raw hypothalamus single cell RNA sequencing data for the stressed/baseline mice (Brivio et al., 2023) the same way as described by Ximerakis et al. for their old/young condition mice. Then, anchor integration is performed to first combine old/young conditions on their own and stress/baseline conditions on their own, and finally, to combine all age and stress conditions. After integration, principal component analysis (PCA) is performed using 20 PCs, followed by nearest neighbors and clustering analyses, allowing for visualization of all combined data on a UMAP. Suerat v5 (Hao et al., 2023) is utilized in RStudio to perform the analysis.

## Data
The first dataset (Ximerakis et al., 2019) contains processed single cell RNA sequencing data from the whole brains of young and old mice. 

The second dataset (Brivio et al., 2023) contains raw 10x files ("barcodes", "features", and "matrix") from the hypothalamuses of either baseline or chronic mild stress condition mice either untreated or subjected to acute restraint stress.


## Folder structure
#### Code:
This folder contains "20440_Project_RScript.R", the R script which performs the single cell RNA sequencing processing, integration, and dimensional reduction, and then finally creates the UMAP figure. 
  
#### Data:
Since the data is too big to be uploaded, there is a Data.md file specifying where to download the necessary data. The following describes how the Data folder should be created once the data is downloaded in order to be able to run with the above code.
  ##### data_old:
  This folder contains the pre-processed 10X .txt.gz files from Ximerakis et al. of single cell RNAseq reads from aged mice.
  
  ##### data_young:
  This folder contains the pre-processed 10X .txt.gz files from Ximerakie et al. of single cell RNAseq reads from young mice.
  
  ##### data_stress:
  ###### M_baseline_ARS:
  This folder contains the barcodes, features, and matrix files from the raw 10X single cell RNAseq hypothalamus data from      the Brivio et al. for the baseline mice subjected to acute restraint stress.
  ###### M_baseline_ctrl:
  This folder contains the barcodes, features, and matrix files from the raw 10X single cell RNAseq hypothalamus data from      the Brivio et al. for the baseline control condition mice.
  ###### M_CMS_ARS:
  This folder contains the barcodes, features, and matrix files from the raw 10X single cell RNAseq hypothalamus data from      the Brivio et al. for the chronic mild stress mice subjected to acute restraint stress.
  ###### M_CMS_ctrl:
  This folder contains the barcodes, features, and matrix files from the raw 10X single cell RNAseq hypothalamus data from      the Brivio et al. for the chronic mild stress control condition mice.
  
#### Figures:
This folder contains the Stress_age_combined_umap.png, the output of the code stored in the "Code" folder. This figure is a UMAP produced after principal component analysis (PCA) is performed on all conditions of the integrated data: old (whole brain), young (whole brain), baseline control (hypothalamus), chronic mild stress control (hypothalamus), baseline subjected to acute restraint stress (hypothalamus), and chronic mild stress subjected to acute restraint stress(hypothalamus). The figure plots UMAP dimension 1 vs. UMAP dimension 2, and results in an outputted 42 numbered clusters, which are visualized by color.

## Installation
To run this code, first set up a file structure comparable to that of the file structure on this repository. Then, download the "20440_Project_RScript.R" file and all of the data contained in the "Data" folder. Open the "20440_Project_RScript.R" file in RStudio (or other preferred method of running R scripts) and update the "setwd()" line to specify the appropriate file path to your the R script. From here, running the R script should output the UMAP of the data integrated from all mouse conditions. 

## Citations
Brivio, E., Kos, A., Ulivi, A. F., Karamihalev, S., Ressle, A., Stoffel, R., Hirsch, D., Stelzer, G., Schmidt, M. V., Lopez, J. P., & Chen, A. (2023). Sex shapes cell-type-specific transcriptional signatures of stress exposure in the mouse hypothalamus. Cell Reports, 42(8), 112874.

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2023). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology. https://doi.org/10.1038/s41587-023-01767-y

Ximerakis, M., Lipnick, S. L., Innes, B. T., Simmons, S. K., Adiconis, X., Dionne, D., Mayweather, B. A., Nguyen, L., Niziolek, Z., Ozek, C., Butty, V. L., Isserlin, R., Buchanan, S. M., Levine, S. S., Regev, A., Bader, G. D., Levin, J. Z., & Rubin, L. L. (2019). Single-cell transcriptomic profiling of the aging mouse brain. Nature Neuroscience, 22(10), 1696–1708.





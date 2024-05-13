# Install Seurat and load the packages
install.packages('Seurat')
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("mojaveazure/seurat-disk")
install.packages('Matrix')
install.packages('dplyr')
install.packages('magrittr')
install.packages('ggplot2')
install.packages('patchwork')
install.packages('scCustomize')
BiocManager::install("MAST")
BiocManager::install('EnhancedVolcano')
install.packages('devtools')
devtools::install_github('kevinblighe/EnhancedVolcano')
install.packages('scales')
install.packages('viridis')
install.packages('RColorBrewer')
install.packages('colorspace')
install.packages('tidyverse')
install.packages('reshape2')
install.packages('ggplot2')
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Matrix)
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(MAST)
library(EnhancedVolcano)
library(scales)
library(viridis)
library(RColorBrewer)
library(colorspace)
library(tidyverse)
library(reshape2)
library(ggplot2)


#### Set working directory (*different for each user*)
setwd("/Users/emmafinburgh/20440_Project/Code")


#### Combine the age data

# Retrieve list of files from the directory for young condition (already processed)
file_list_young <- list.files(path = '../Data/data_young', pattern = "*.txt.gz", full.names = TRUE)
seurat_young.list <- list()

# Iterate through the young mouse data using seq_along
for (i in seq_along(file_list_young)){
  # Assign file paths
  file_path <- file_list_young[i]
  # Read the sample data from .txt.gz files, convert to a matrix
  sample_data <- read.table(file = file_path, header = TRUE, sep = '\t')
  sample_data <- as.matrix(sample_data)
  sample_data <- Matrix(sample_data, sparse = TRUE)
  # Create a Seurat object for each sample
  sample_seurat <- CreateSeuratObject(counts = sample_data, project = paste0(basename(file_path)))
  # Assign the gene expression data to the "RNA" assay of the seurat object
  sample_seurat[["RNA"]]$data <- sample_data
  seurat_young.list[[i]] <- sample_seurat
}

# Retrieve list of files from the directory for old condition (already processed)
file_list_old <- list.files(path = '../Data/data_old', pattern = "*.txt.gz", full.names = TRUE)
seurat_old.list <- list()

# Iterate through the old mouse data using seq_along
for (i in seq_along(file_list_old)){
  # Assign file paths
  file_path <- file_list_old[i]
  # Read the sample data from .txt.gz files, convert to a matrix
  sample_data <- read.table(file = file_path, header = TRUE, sep = '\t')
  sample_data <- as.matrix(sample_data)
  sample_data <- Matrix(sample_data, sparse = TRUE)
  # Create a Seurat object for each sample
  sample_seurat <- CreateSeuratObject(counts = sample_data, project = paste0(basename(file_path)))
  # Assign the gene expression data to the "RNA" assay of the seurat object
  sample_seurat[["RNA"]]$data <- sample_data
  seurat_old.list[[i]] <- sample_seurat
}
  
# Find and integrate anchors for young data
age_young.anchors <- FindIntegrationAnchors(object.list = seurat_young.list, dims = 1:30)
age_young.integrated <- IntegrateData(anchorset = age_young.anchors, dims = 1:30)
saveRDS(age_young.integrated, file.path("../RDS_files", "age_young.rds"))

# Find and integrate anchors for old data
age_old.anchors <- FindIntegrationAnchors(object.list = seurat_old.list, dims = 1:30)
age_old.integrated <- IntegrateData(anchorset = age_old.anchors, dims = 1:30)
saveRDS(age_old.integrated, file.path("../RDS_files", "age_old.rds"))

# Combine old and young data by finding and integrating anchored data
age.list<-c(age_old.integrated, age_young.integrated)
age_combined.anchors <- FindIntegrationAnchors(object.list = age.list, dims = 1:30)
age_combined.integrated <- IntegrateData(anchorset = age_combined.anchors, dims = 1:30)
saveRDS(age_combined.integrated, file.path("../RDS_files", "age_combined.rds"))



#### Process raw hypothalamus stress data

# Read in the 4 conditions of hypothalamus/stress data
CMS_ctrl <- Read10X('../Data/data_stress/M_CMS_ctrl')
baseline_ctrl <- Read10X('../Data/data_stress/M_baseline_ctrl')
CMS_ARS <- Read10X('../Data/data_stress/M_CMS_ARS')
baseline_ARS <- Read10X('../Data/data_stress/M_baseline_ARS')

# Create Seurat objects for each hypothalamus/stress conditions
# Exclude reads with <3 cells and cells with <200 genes
CMS_ctrl <- CreateSeuratObject(counts = CMS_ctrl, project = '../Data/data_stress/M_CMS_ctrl', min.cells = 3, min.features = 200)
baseline_ctrl <- CreateSeuratObject(counts = baseline_ctrl, project = '../Data/data_stress/M_baseline_ctrl', min.cells = 3, min.features = 200)
CMS_ARS <- CreateSeuratObject(counts = CMS_ARS, project = '../Data/data_stress/M_CMS_ARS_ctrl', min.cells = 3, min.features = 200)
baseline_ARS <- CreateSeuratObject(counts = baseline_ARS, project = '../Data/data_stress/M_baseline_ARS_ctrl', min.cells = 3, min.features = 200)

# Extract % mitochondrial DNA for hypothalamus/stress conditions
CMS_ctrl[["percent.mt"]] <- PercentageFeatureSet(CMS_ctrl, pattern = "^mt-")
baseline_ctrl[["percent.mt"]] <- PercentageFeatureSet(baseline_ctrl, pattern = "^mt-")
CMS_ARS[["percent.mt"]] <- PercentageFeatureSet(CMS_ARS, pattern = "^mt-")
baseline_ARS[["percent.mt"]] <- PercentageFeatureSet(baseline_ARS, pattern = "^mt-")

# Filter the hypothalamus/stress data as performed on the age data
# Set minimum number of genes = 250, maximum number of genes = 6000, and maximum % mitochondrial DNA = 30%
CMS_ctrl <- subset(CMS_ctrl, subset = nFeature_RNA > 250 & nFeature_RNA < 6000  & percent.mt < 30)
baseline_ctrl <- subset(baseline_ctrl, subset = nFeature_RNA > 250 & nFeature_RNA < 6000  & percent.mt < 30)
CMS_ARS <- subset(CMS_ARS, subset = nFeature_RNA > 250 & nFeature_RNA < 6000  & percent.mt < 30)
baseline_ARS <- subset(baseline_ARS, subset = nFeature_RNA > 250 & nFeature_RNA < 6000  & percent.mt < 30)

# log-normalize the hypothalamus/stress data
hypo.list<-c(CMS_ctrl, baseline_ctrl, CMS_ARS, baseline_ARS)
for (i in 1:4) {
  hypo.list[[i]] <- NormalizeData(hypo.list[[i]], verbose = FALSE)
  hypo.list[[i]] <- FindVariableFeatures(hypo.list[[i]], x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)
}

# Combine hypothalamus/stress data by finding and integrating anchored data
hypo.anchors <- FindIntegrationAnchors(object.list = hypo.list, dims = 1:30)
hypo.integrated <- IntegrateData(anchorset = hypo.anchors, dims = 1:30)
saveRDS(hypo.integrated, file.path("../RDS_files", "stress_hypo.rds"))



#### Combine hypothalamus/stress and age data

## First combine stress and young conditions
stress_young.list <- c(hypo.integrated, age_young.integrated)
stress_young.anchors <- FindIntegrationAnchors(object.list = stress_young.list, dims = 1:30)
stress_young.integrated <- IntegrateData(anchorset = stress_young.anchors, dims = 1:30)
# Scale the data
DefaultAssay(stress_young.integrated) <- "RNA"
stress_young.integrated <- ScaleData(stress_young.integrated, verbose = FALSE)
DefaultAssay(stress_young.integrated) <- "integrated"
stress_young.integrated <- ScaleData(stress_young.integrated, verbose = FALSE)
# Perform a PCA with 20 PCs
stress_young.integrated <- RunPCA(stress_young.integrated, npcs = 20, verbose = FALSE) # PCs specified from Ximerakis et al
# Visualize the proportion of variance explained by each PC
ElbowPlot(stress_young.integrated)
# Conduct nearest neighbors analysis and find clusters
stress_young.integrated <- FindNeighbors(stress_young.integrated, dims = 1:20)
stress_young.integrated <- FindClusters(stress_young.integrated, resolution = 1.6) # resolution specified from Ximerakis et al.
# Create a UMAP from the PCA reduction
stress_young.integrated <- RunUMAP(stress_young.integrated, reduction = "pca", dims = 1:20)
saveRDS(stress_young.integrated, file.path("../RDS_files", "stress_young.rds"))
# Save the UMAP
stress_young_umap <- DimPlot(stress_young.integrated)
ggsave("../Figures/Stress_young_umap.png", stress_young_umap, dpi = 300)

## Then combine stress and old conditions
stress_old.list <- c(hypo.integrated, age_old.integrated)
stress_old.anchors <- FindIntegrationAnchors(object.list = stress_old.list, dims = 1:30)
stress_old.integrated <- IntegrateData(anchorset = stress_old.anchors, dims = 1:30)
# Scale the data
DefaultAssay(stress_old.integrated) <- "RNA"
stress_old.integrated <- ScaleData(stress_old.integrated, verbose = FALSE)
DefaultAssay(stress_old.integrated) <- "integrated"
stress_old.integrated <- ScaleData(stress_old.integrated, verbose = FALSE)
# Perform a PCA with 20 PCs
stress_old.integrated <- RunPCA(stress_old.integrated, npcs = 20, verbose = FALSE) # PCs specified from Ximerakis et al
# Visualize the proportion of variance explained by each PC
ElbowPlot(stress_old.integrated)
# Conduct nearest neighbors analysis and find clusters
stress_old.integrated <- FindNeighbors(stress_old.integrated, dims = 1:20)
stress_old.integrated <- FindClusters(stress_old.integrated, resolution = 1.6) # resolution specified from Ximerakis et al.
# Create a UMAP from the PCA reduction
stress_old.integrated <- RunUMAP(stress_old.integrated, reduction = "pca", dims = 1:20)
saveRDS(stress_old.integrated, file.path("../RDS_files", "stress_old.rds"))
# Save the UMAP
stress_old_umap <- DimPlot(stress_old.integrated)
ggsave("../Figures/Stress_old_umap.png", stress_old_umap, dpi = 300)

## Finally, combine all data
stress_age.list <-c(hypo.integrated, age_combined.integrated)
stress_age.anchors <- FindIntegrationAnchors(object.list = stress_age.list, dims = 1:30)
stress_age.integrated <- IntegrateData(anchorset = stress_age.anchors, dims = 1:30)
# Scale the data
DefaultAssay(stress_age.integrated) <- "RNA"
stress_age.integrated <- ScaleData(stress_age.integrated, verbose = FALSE)
DefaultAssay(stress_age.integrated) <- "integrated"
stress_age.integrated <- ScaleData(stress_age.integrated, verbose = FALSE)
# Perform a PCA with 20 PCs
stress_age.integrated <- RunPCA(stress_age.integrated, npcs = 20, verbose = FALSE) # PCs specified from Ximerakis et al
# Visualize the proportion of variance explained by each PC
ElbowPlot(stress_age.integrated)
# Conduct nearest neighbors analysis and find clusters
stress_age.integrated <- FindNeighbors(stress_age.integrated, dims = 1:20)
stress_age.integrated <- FindClusters(stress_age.integrated, resolution = 1.6) # resolution specified from Ximerakis et al.
# Create a UMAP from the PCA reduction
stress_age.integrated <- RunUMAP(stress_age.integrated, reduction = "pca", dims = 1:20)
saveRDS(stress_age.integrated, file.path("../RDS_files", "stress_age.rds"))
# Save the UMAP
stress_age_umap <- DimPlot(stress_age.integrated)
ggsave("../Figures/Stress_age_combined_umap.png", stress_age_umap, width = 8, height = 6, dpi = 300)



#### Cell type cluster labeling

# Find top markers for each cluster
stress_age <- readRDS("../RDS_files/stress_age.rds")
stress_age.markers <- FindAllMarkers(stress_age, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Find the top 10 upregulated genes per cluster
top10 <- stress_age.markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)
# Narrow in to the top 3 upregulated genes for cell type identification purposes
top3 <- stress_age.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# Rename each cluster with top3 markers
top3 <- top3 %>%  group_by(cluster) %>%  
  dplyr::mutate(markers = paste0(gene, collapse = "/")) %>% dplyr::slice(1)  
marker.names <- top3$markers
current.cluster.ids <- as.character(0:(length(unique(Idents(stress_age)))-1))
new.cluster.ids <- marker.names
Idents(stress_age) <- plyr::mapvalues(Idents(stress_age),
                                    from = current.cluster.ids, to = new.cluster.ids)

# Save the top ten upregulated genes per cluster as a .csv file
write.csv(top10, file.path("../CSV_files", "Stress_age_cluster_markers.csv"))
stress_age$orig.cluster<-Idents(stress_age)

# Create a cluster dataframe and label the clusters
clustavg<-AverageExpression(stress_age)
clustdf<-as.data.frame(clustavg$RNA)
clustdf<-t(clustdf)
clustdf<-as_tibble(clustdf, rownames = "cluster")
head(clustdf)
dim(clustdf)

# List the cell identification markers (genes) from Ximerakis et al.
marks<-c("Pdgfra","Cldn11","Npy","Thbs4","Cd44",
         "Gja1",
         "Cdk1",
         "Sox11",
         "Syt1",
         "Baiap3",
         "Ccdc153",
         "Sspo",
         "Rax",
         "Ttr",
         "Cldn5",
         "Kcnj8",
         "Acta2",
         "Alas2",
         "Slc6a13",
         "Tmem119",
         "Plac8",
         "Pf4")

# Update the dataframe
df2 <- clustdf %>% dplyr::select("cluster",all_of(marks))
df2$orig_cluster<-df2$cluster
# Save the updated dataframe as a .csv file
write.csv(df2, file.path("../CSV_files", "Clusters_by_cell_type_df.csv"))

# Define the cell specific gene markers and assign clusters to cell types (Based on Ximerakis et al.)
df2<- df2 %>% mutate(ID = ifelse(Cldn11 > 2, "OLG", "Unlabeled"))
df2<- df2 %>% mutate(ID = ifelse(Pdgfra > 2, "OPC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Npy > 3, "OEG", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Thbs4 > 1, "NSC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Cd44 > 1, "ARP", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Gja1 > 1, "ASC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Gja1 > 1 & Cldn11 > 1, "ABC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Sox11 > 2, "ImmN", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Sox11 > 2 & Cdk1 > 1, "NRP", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Syt1 > 0.9, "mNEUR", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Syt1 > 0.9 & Baiap3 > 0.5, "PC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Ccdc153 > 3, "EPC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Ccdc153 > 1 & Gja1 > 1, "EPC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Sspo > 3, "VLMC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Rax > 3, "VLMC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Ttr > 3, "CPC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Cldn5 > 3, "EC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Cldn5 > 3 & Alas2 > 1, "Hb_VC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Cldn5 >1 & Gja1 > 1, "OEG", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Kcnj8 > 2, "PC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Kcnj8 > 1 & Cldn5 > 3,  "PC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Acta2 > 2, "VSMC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Acta2 > 1 & Cldn5 > 3, "VSMC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Slc6a13 > 2, "VLMC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Tmem119 > 2, "MG", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Plac8 > 3, "MNC", df2$ID))
df2<- df2 %>% mutate(ID = ifelse(Pf4 > 2, "MAC", df2$ID))

y<-dplyr::select(df2,cluster)
Idents(stress_age)<-plyr::mapvalues(Idents(stress_age), from = df2$cluster, to = df2$ID)
stress_age$group<-Idents(stress_age)

# Reorder groupings
levels(stress_age) <- c("OLG",
                      "OPC",
                      "OEG", 
                      "NSC",
                      "ARP",
                      "ASC",
                      "ABC",
                      "ImmN",
                      "NRP",
                      "mNEUR",
                      "PC",
                      "EPC",
                      "VLMC",
                      "CPC",
                      "EC",
                      "Hb_VC",
                      "VSMC",
                      "MG",
                      "MNC",
                      "MAC",
                      "Unlabeled")

# Replot UMAP with cell type annotations
stress_age_labeled_umap <- DimPlot(stress_age)
ggsave("../Figures/Stress_age_labeled_umap.png", stress_age_labeled_umap, dpi = 300)

# Label the umap with defined colors
cluster_colors <- c("OLG" = "#999800", "OPC" = "#0000FF", "OEG" = "#705600", "ASC" = "#700180",
                    "ABC" = "#909099", "ImmN" = "#00FFFF", "NRP" = "#9000F8", "mNEUR" = "#0FFF00",
                    "PC" = "#008080", "EPC" = "#FF00FF", "VLMC" = "#FF0000", "CPC" = "#FFF000",
                    "EC" = "#A32A2A", "Hb_VC" = "#000194", "VSMC" = "#FFC0CB", "MG" = "#099400", "MAC" = "#FF8C00")
                    
stress_age_labeled_colors_umap <- DimPlot(stress_age, cols = cluster_colors)
# Save the labeled UMAP
ggsave("../Figures/Stress_age_labeled_colors_umap.png", stress_age_labeled_colors_umap, dpi = 300)
saveRDS(stress_age, file.path("../RDS_files", "stress_age_cell_types.rds"))



#### Label the umap by mouse condition
## Baseline_ARS umap
find_baseline_ARS_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("data_stress/M_baseline_ARS_ctrl"))
baseline_ARS_cells <- list("baseline_ARS"= find_baseline_ARS_cells)
baseline_ARS_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = baseline_ARS_cells,
                                         highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/baseline_ARS_umap.png", baseline_ARS_umap, dpi = 300)

## Baseline_ctrl umap
find_baseline_ctrl_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("data_stress/M_baseline_ctrl"))
baseline_ctrl_cells <- list("baseline_ctrl"= find_baseline_ctrl_cells)
baseline_ctrl_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = baseline_ctrl_cells,
                                         highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/baseline_ctrl_umap.png", baseline_ctrl_umap, dpi = 300)

## CMS_ARS umap
find_CMS_ARS_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("data_stress/M_CMS_ARS_ctrl"))
CMS_ARS_cells <- list("CMS_ARS"= find_CMS_ARS_cells)
CMS_ARS_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = CMS_ARS_cells,
                                         highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/CMS_ARS_umap.png", CMS_ARS_umap, dpi = 300)

## CMS_ctrl umap
find_CMS_ctrl_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("data_stress/M_CMS_ctrl"))
CMS_ctrl_cells <- list("CMS_ctrl"= find_CMS_ctrl_cells)
CMS_ctrl_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = CMS_ctrl_cells,
                                         highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/CMS_ctrl_umap.png", CMS_ctrl_umap, dpi = 300)

## Old umap
find_old_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("X37", "X38", "X33", "X34", "X39", "X40", "X43", "X44"))
old_cells <- list("Old"= find_old_cells)
old_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = old_cells,
                                         highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/old_umap.png", old_umap, dpi = 300)

## Young umap
find_young_cells <- WhichCells(object = stress_age, expression = orig.ident %in% c("X6", "X7", "X19", "X20", "X21", "X22", "X27", "X28"))
young_cells <- list("Young"= find_young_cells)
young_umap <- Cell_Highlight_Plot(seurat_object=stress_age, cells_highlight = young_cells,
                                highlight_color = "Red")
# Save the labeled UMAP
ggsave("../Figures/young_umap.png", young_umap, dpi = 300)



#### Differential Gene Analysis
# Load the stress_age cell type annotated data
Age_data_annotated <- readRDS("../RDS_files/stress_age_cell_types.rds")

# Create column of names for each mouse condition with cell type
Age_data_annotated@meta.data$cell_condition <- rownames(Age_data_annotated@meta.data)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("X7", "X6", "X19", "X20", "X21", "X22", "X27", "X28"), paste("Young", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("X37", "X38", "X33", "X34", "X39", "X40", "X43", "X44"), paste("Old", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_CMS_ctrl"), paste("CMS_ctrl", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_CMS_ARS_ctrl"), paste("CMS_ARS", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_baseline_ctrl"), paste("baseline_ctrl", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
Age_data_annotated@meta.data$cell_condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_baseline_ARS_ctrl"), paste("baseline_ARS", Age_data_annotated@meta.data$group, sep = "_"), Age_data_annotated@meta.data$cell_condition)
# Create column of names for each mouse condition without cell type
Age_data_annotated@meta.data$condition <- rownames(Age_data_annotated@meta.data)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("X7", "X6", "X19", "X20", "X21", "X22", "X27", "X28"), "Young", Age_data_annotated@meta.data$condition)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("X37", "X38", "X33", "X34", "X39", "X40", "X43", "X44"),"Old", Age_data_annotated@meta.data$condition)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_CMS_ctrl"), "CMS_ctrl", Age_data_annotated@meta.data$condition)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_CMS_ARS_ctrl"), "CMS_ARS", Age_data_annotated@meta.data$condition)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_baseline_ctrl"), "baseline_ctrl", Age_data_annotated@meta.data$condition)
Age_data_annotated@meta.data$condition <- ifelse(Age_data_annotated@meta.data$orig.ident %in% c("data_stress/M_baseline_ARS_ctrl"),"baseline_ARS", Age_data_annotated@meta.data$condition)
# Save RDS file with new annotations
saveRDS(Age_data_annotated, file.path("../RDS_files", "stress_age_cell_types_annotated.rds")

        
## Function for performing MAST DGE (adapted from Hajdarovic and Yu et al.)
do_mast <- function(de_obj, ident1, ident2, condition="cell_condition") {
  DefaultAssay(de_obj) <-"integrated"
  Idents(de_obj) <- condition
  run_MAST <- FindMarkers(de_obj, ident.1 = ident1 , ident.2 = ident2, verbose = FALSE, test.use = "MAST", only.pos = FALSE)
  return(run_MAST)
}


## Differential genes expressed per cell type
# Re-run for each identity1 and identity2 of interest
identity1= paste("CMS_ctrl")
identity2 = paste("baseline_ctrl")
conditiontype = "cell_condition"

cell_DE_values <- list()
cell_names <- levels(Age_data_annotated@meta.data$group)

# Loop through all cell types
for (i in seq_along(cell_names)){
  cell_type <- cell_names[[i]]
  identity1_w_cell = paste(identity1, cell_type, sep = "_")
  identity2_w_cell = paste(identity2, cell_type, sep = "_")
  # Perform mast of old vs. young on all cell types
  mast_run <- do_mast(de_obj=Age_data_annotated, ident1=identity1_w_cell,
                      ident2=identity2_w_cell, condition=conditiontype)
  mast_run$cell_type <- cell_type
  # Save the gene names in a table
  mast_run <- as_tibble(mast_run, rownames = "gene")
  cell_DE_values[[i]] <- mast_run
}
# Save DGE as an RDS file
filename1 <- paste0("output_", identity1,"_vs_", identity2)
saveRDS(cell_DE_values, file.path("../RDS_files", paste0(filename1, ".rds")))

# Bind data together
cell_DE_values_merge <- bind_rows(cell_DE_values)
cell_DE_values_merge <- group_by(cell_DE_values_merge, cell_type)
# Save merged DGE values as a csv file
write.csv(cell_DE_values_merge, file.path("../CSV_files/", paste0(filename1, ".csv")))
saveRDS(cell_DE_values_merge, file.path("../RDS_files", paste0(filename1, "_merge.rds")))

# Define a color palette
require(scales)
cell_type <- levels(Age_data_annotated@meta.data$group)
cols <- hue_pal()(length(cell_type))
col_data <- as_tibble(cbind(cols,cell_type))
cols

# Merge the color data with the DGE values
cell_DE_values_merge <-full_join(cell_DE_values_merge, col_data, by = "cell_type")

# Color only the differentially expressed genes (|log2FC| > 1 and p < 0.05)
cell_DE_values_merge <-mutate(cell_DE_values_merge, FC = ifelse(avg_log2FC <=-1 | avg_log2FC >=1, TRUE, FALSE))
cell_DE_values_merge <-mutate(cell_DE_values_merge, pval_yes = ifelse(p_val_adj <0.05, TRUE, FALSE))
cell_DE_values_merge <-mutate(cell_DE_values_merge, fillyn = ifelse(FC== TRUE & pval_yes==TRUE, TRUE, FALSE))

cell_DE_values_merge <- mutate(cell_DE_values_merge, col_use = ifelse(fillyn==TRUE, cols, "dark gray"))
cell_DE_values_merge <-mutate(cell_DE_values_merge, updown = ifelse(p_val_adj <0.05 & avg_log2FC >=1, "Upregulated", ifelse(p_val_adj <0.05 & avg_log2FC <=-1, "Downregulated", "NS")))
cell_DE_values_merge
# Put the colored cell types in a DE table
DE_table <- table(cell_DE_values_merge$updown, cell_DE_values_merge$cell_type)
DE_table <- as.data.frame.matrix(DE_table)
DE_table
# Save DE by cell type as a .csv file
write.csv(DE_table, file.path("../CSV_files/", paste0("DE_table_", identity1,"_vs_", identity2, ".csv")))

plotting_cols <- as.character(cell_DE_values_merge$col_use)

# Assign cell type as x-values
x_val <- unique(cell_DE_values_merge$cell_type)
cell_DE_values_merge$cell_type <- factor(cell_DE_values_merge$cell_type, levels=x_val)
table(cell_DE_values_merge$fillyn, cell_DE_values_merge$cell_type)

# Plotting DE as a Strip plot
strip_plot <- ggplot(cell_DE_values_merge, aes(cell_type, y = avg_log2FC)) +
  geom_jitter(position = position_jitter(0.3), color = plotting_cols, alpha = 0.5) +
  geom_text(aes(label = ifelse(fillyn & (avg_log2FC >3 | avg_log2FC < -3), gene, "")), position = position_jitter(0.3), vjust = -0.5, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
strip_plot
filename2 <- paste0("../Figures/", identity1,"_vs_", identity2, "_strip_plot.png")
ggsave(filename2, strip_plot, width = 8, height = 6, dpi = 500)



#### Volcano plots

## Stress conditions vs. age
stress_conditions <- c("CMS_ctrl", "CMS_ARS", "baseline_ctrl", "baseline_ARS")
ages <- c("Young", "Old")

for (stress_condition in stress_conditions) {
  for (age in ages){
    # Prepare the identifiers
    ident1 <- stress_condition
    ident2 <- age
    
    # Run MAST
    mast_run <- do_mast(de_obj = Age_data_annotated, ident1 = ident1, ident2 = ident2, condition = "condition")
    mast_run_df <- as.data.frame(mast_run)
    
    # Create volcano plot
    volcano_plot <- EnhancedVolcano(mast_run_df, lab = rownames(mast_run_df), x = 'avg_log2FC', y = 'p_val_adj')
    
    # Save volcano plot
    ggsave(paste("../Figures/", stress_condition, "_", age, "_volcano_plot.png", sep = ""), volcano_plot, dpi = 500)
  }
}


## Cell type volcano plots
var_conditions_1 <- c("Old") 
var_conditions_2 <- c("Young")
volcano_cells <- c("ASC", "OLG", "MG","mNEUR", "OPC", "PC", "EC", "EPC")

for (var_condition_1 in var_conditions_1) {
  for (var_condition_2 in var_conditions_2){
    for (vol_cell in volcano_cells){
      # Prepare the identifiers
      ident1 <- paste(var_condition_1, vol_cell, sep = "_")
      ident2 <- paste(var_condition_2, vol_cell, sep = "_")
      
      # Run MAST
      mast_run <- do_mast(de_obj = Age_data_annotated, ident1 = ident1, ident2 = ident2, condition = "cell_condition")
      mast_run_df <- as.data.frame(mast_run)
      
      # Create volcano plot
      volcano_plot <- EnhancedVolcano(mast_run_df, lab = rownames(mast_run_df), x = 'avg_log2FC', y = 'p_val_adj', legendPosition = 'right')
      
      # Save volcano plot
      ggsave(paste("../Figures/", ident1, "_", ident2, "_volcano_plot.png", sep = ""), volcano_plot, width = 8, height = 6, dpi = 500)
    }
  }
}


## Labeled cell type Volcano plots
var_conditions_1 <- c("CMS_ctrl") # Switch to "Old"
var_conditions_2 <- c("baseline_ctrl") # Switch to "Young"
volcano_cells <- c("EPC") # Switch to cell type of interest

# Overlapping genes per cell type (choose based on cell type above)
ASC_lab = c("Gpr83", "Itih2", "Krt15", "Cd34", "D930028M14Rik",
  "Cxcl2", "Higd1b", "Avp", "Cd37", "Tek", "Bcl2a1b",
  "Syndig1l", "Spp1", "Gkn3", "Il1a", "BC028528",
  "Hbb-bt", "Hba-a2", "Hbb-bs")

OLG_lab = c("Otp", "Pdyn", "Crispld2", "Pmch", "Npas4",
            "Grp", "Foxq1")

MG_lab = c("Sntn", "Hdc", "Lgals3", "Myo16", "Syndig1l",
           "Cp", "Myh11", "Slc6a13", "Avp", "Fibcd1",
           "Pdyn", "Cenpf", "Pam", "Ccl7")

OPC_lab = c("Pdyn", "2410004P03Rik", "Calml4", "Sfrp2", "Casq2",
            "Htr3a", "Syt10", "Ucma", "Fam183b", "Hpgds",
            "Nos1", "Ccl4", "Cldn5")

PC_lab = c("Tagap", "Cxcl10", "Thbs1", "Nkx2-1", "Dmkn", "A2m", "Matn1")

EC_lab = c("Oprk1", "Cadps", "Fibcd1", "Prlr", "Avp",
           "Gabra1", "Stx1b", "Myt1l", "Nxph1")

EPC_lab = c("Gpr88", "Tac1", "Cxcl10", "Bin2", "Adgre1",
            "Tek", "Tmem119", "Ptprb", "Avp", "Sst",
            "Opcml", "Nhlh2", "Anln")

# Loop through conditions to compare
for (var_condition_1 in var_conditions_1) {
  for (var_condition_2 in var_conditions_2){
    for (vol_cell in volcano_cells){
      # Prepare the identifiers
      ident1 <- paste(var_condition_1, vol_cell, sep = "_")
      ident2 <- paste(var_condition_2, vol_cell, sep = "_")
      
      # Run MAST
      mast_run <- do_mast(de_obj = Age_data_annotated, ident1 = ident1, ident2 = ident2, condition = "cell_condition")
      mast_run_df <- as.data.frame(mast_run)
      
      # Create volcano plot
      volcano_plot <- EnhancedVolcano(mast_run_df, lab = rownames(mast_run_df), 
                                      x = 'avg_log2FC', y = 'p_val_adj', legendPosition = 'right',
                                      selectLab = EPC_lab, # specify overlapping genes per cell type
                                      drawConnectors = TRUE)
      
      # Save volcano plot
      ggsave(paste("../Figures/", ident1, "_", ident2, "_volcano_plot_annotated.png", sep = ""), volcano_plot, width = 8, height = 6, dpi = 500)
    }
  }
}



#### FeaturePlots
# Define features (genes) of interest
features <- c("Cd34", "Rgs1", "Adgrf5", "Cdh5", "Cxcl10", "Il1a", "Pdyn",
              "Spp1", "Higd1b", "Avp", "Tek", "Tagap", "Thbs1", "Krt15", 
              "Gpr88", "Nkx2-1", "Bin2", "Cldn5", "Gpr88", "Tac1", "Tmem119",
              "Dmkn", "Syt10", "Casq2")

# Define subsets to plot on a feature plot
CMS_subset <- subset(x=Age_data_annotated, subset = condition %in% c("CMS_ctrl"))
baseline_subset <- subset(x=Age_data_annotated, subset = condition %in% c("baseline_ctrl"))
Old_subset <- subset(x=Age_data_annotated, subset = condition %in% c("Old"))
Young_subset <- subset(x=Age_data_annotated, subset = condition %in% c("Young"))
stress_subset <- subset(x=Age_data_annotated, subset = condition %in% c("CMS_ctrl", "baseline_ctrl"))
age_subset <- subset(x=Age_data_annotated, subset = condition %in% c("Old", "Young"))

subsets <- list(CMS = CMS_subset, baseline = baseline_subset, Old = Old_subset, Young = Young_subset)
age_stress_subsets <- list(age_sub = age_subset, stress_sub = stress_subset)

# Loop through subsets and plot defined features on the UMAP
for (subset_name in names(subsets)) {
  subset <- subsets[[subset_name]]
  for (feature in features) {
    gene_feature_plot <- FeaturePlot(subset, features = feature, min.cutoff = -1, max.cutoff = 4)
    ggsave(paste0("../Figures/FeaturePlots_", subset_name, "_", feature, ".png"), gene_feature_plot, width = 6, height = 6, dpi = 300)
  }
}




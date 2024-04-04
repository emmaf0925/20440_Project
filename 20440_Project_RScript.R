# Install Seurat and load the packages
install.packages('Seurat')
install.packages("remotes")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("mojaveazure/seurat-disk")
install.packages('Matrix')
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Matrix)

# Set working directory (*different for each user*)
setwd("/Users/emmafinburgh/20440_Project/Code")

# Combine the age data

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

# Process raw hypothalamus stress data

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

# log-normalize the 
hypo.list<-c(CMS_ctrl, baseline_ctrl, CMS_ARS, baseline_ARS)
for (i in 1:4) {
  hypo.list[[i]] <- NormalizeData(hypo.list[[i]], verbose = FALSE)
  hypo.list[[i]] <- FindVariableFeatures(hypo.list[[i]], x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5)
}

# Combine hypothalamus/stress data by finding and integrating anchored data
hypo.anchors <- FindIntegrationAnchors(object.list = hypo.list, dims = 1:30)
hypo.integrated <- IntegrateData(anchorset = hypo.anchors, dims = 1:30)
saveRDS(hypo.integrated, file.path("../RDS_files", "stress_hypo.rds"))

# Combine hypothalamus/stress and age data

# First combine stress and young conditions
stress_young.list <- c(hypo.integrated, age_young.integrated)
stress_young.anchors <- FindIntegrationAnchors(object.list = stress_young.list, dims = 1:30)
stress_young.integrated <- IntegrateData(anchorset = stress_young.anchors, dims = 1:30)
# Scale the data
DefaultAssay(stress_young.integrated) <- "RNA"
stress_young.integrated <- ScaleData(stress_young.integrated, verbose = FALSE)
DefaultAssay(stress_young.integrated) <- "integrated"
stress_young.integrated <- ScaleData(stress_young.integrated, verbose = FALSE)


# Then combine stress and old conditions
stress_old.list <- c(hypo.integrated, age_old.integrated)
stress_old.anchors <- FindIntegrationAnchors(object.list = stress_old.list, dims = 1:30)
stress_old.integrated <- IntegrateData(anchorset = stress_old.anchors, dims = 1:30)
# Scale the data
DefaultAssay(stress_old.integrated) <- "RNA"
stress_old.integrated <- ScaleData(stress_old.integrated, verbose = FALSE)
DefaultAssay(stress_old.integrated) <- "integrated"
stress_old.integrated <- ScaleData(stress_old.integrated, verbose = FALSE)


# Finally, combine all data
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
png("../Figures/Stress_age_combined_umap.png")
stress_age_umap
dev.off()


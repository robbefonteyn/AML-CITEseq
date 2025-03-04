#install.packages(c("sp","cli","spam", "leidenbase"), repos='http://cran.us.r-project.org')
# Packages installation
#install.packages("Seurat", repos='http://cran.us.r-project.org')


library(Seurat)
library(dplyr)
library(ggplot2)

message(paste("Loading merged object..."))
path <- "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_Full_with_meta.rds"
merged_SAB <- readRDS(path)

message(paste("Finding variable features..."))
merged_SAB <- FindVariableFeatures(merged_SAB, selection.method = "vst", nfeatures = 2000)
gc()

message(paste("Scaling data..."))
merged_SAB <- ScaleData(merged_SAB, features = rownames(merged_SAB))
gc()

message(paste("Running PCA..."))
merged_SAB <- RunPCA(merged_SAB, npcs = 30)
gc()

message(paste("Running UMAP..."))
merged_SAB <- RunUMAP(merged_SAB, dims = 1:30)


message(paste("Running tSNE..."))
merged_SAB <- RunTSNE(merged_SAB, dims = 1:30)

#message(paste("Saving the new object"))
#saveRDS(merged_SAB, file="/kyukon/scratch/gent/458/vsc45888/SAB_Full_with_meta_clustering.rds")

message(paste("Generating UMAP plot..."))
DimPlot(merged_SAB, reduction = "umap", group.by = "Patient", label = TRUE) + ggtitle("UMAP Clustering")
ggsave("UMAP_clustering_merged_data_patient.jpeg", width=1600, height=900, units="px", scale=2.5)

DimPlot(merged_SAB, reduction = "umap",group.by = "exp", label = TRUE) + ggtitle("UMAP Clustering")
ggsave("UMAP_clustering_merged_data_exp.jpeg", width=1600, height=900, units="px", scale=2.5)

DimPlot(merged_SAB, reduction = "umap",group.by = "day", label = TRUE) + ggtitle("UMAP Clustering")
ggsave("UMAP_clustering_merged_data_day.jpeg", width=1600, height=900, units="px", scale=2.5)

DimPlot(merged_SAB, reduction = "umap",group.by = "predicted.celltype.l2", label = TRUE) + ggtitle("UMAP Clustering")
ggsave("UMAP_clustering_merged_data_celltype.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Generating tSNE plot..."))
DimPlot(merged_SAB, reduction = "tsne", label = TRUE) + ggtitle("tSNE Clustering")
ggsave("tSNE_clustering_merged_data.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Generating feature plot for CD34..."))
FeaturePlot(merged_SAB, features = "CD34", reduction = "umap")
ggsave("CD34_feature_plot.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Generating violin plot for CD34 expression..."))
VlnPlot(merged_SAB, features = "CD34", group.by = "seurat_clusters")
ggsave("CD34_violin_plot.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Generating ridge plot for CD34 expression..."))
RidgePlot(merged_SAB, features = "CD34", group.by = "seurat_clusters")
ggsave("CD34_ridge_plot.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Switching to ADT assay..."))
DefaultAssay(merged_SAB) <- "ADT"
message(paste("Finding variable features in ADT assay..."))
merged_SAB <- FindVariableFeatures(merged_SAB, selection.method = "vst", nfeatures = 50, assay= "ADT")
message(paste("Top variable features:"))
print(head(VariableFeatures(merged_SAB)))  # Show the top variable features
VariableFeaturePlot(merged_SAB)
ggsave("ADT_variable_features.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Switching to RNA assay..."))
DefaultAssay(merged_SAB) <- "RNA"
message(paste("Finding variable features in RNA assay..."))
merged_SAB <- FindVariableFeatures(merged_SAB)
VariableFeaturePlot(merged_SAB)
ggsave("RNA_variable_features.jpeg", width=1600, height=900, units="px", scale=2.5)

message(paste("Script completed successfully!"))


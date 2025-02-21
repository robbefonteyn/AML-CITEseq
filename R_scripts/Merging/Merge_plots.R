library(Seurat)
library(dplyr)
library(ggplot2)

# Plot of the merged object
path <- "C:/Users/Robbe Fonteyn/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/SAB_Full.rds"
merged_SAB <- readRDS(path)

merged_SAB <- FindVariableFeatures(merged_SAB, selection.method = "vst", nfeatures = 2000)

merged_SAB <- ScaleData(merged_SAB, features = rownames(merged_SAB))

# Plot of the merged object
merged_SAB <- RunPCA(merged_SAB, npcs = 30)

# UMAP plot
merged_SAB <- RunUMAP(merged_SAB, dims = 1:30)

# tSNE plot
merged_SAB <- RunTSNE(merged_SAB, dims = 1:30)


DimPlot(merged_SAB, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering")
ggsave("UMAP_clustering_merged_data", width=1600, height=900, units="px", scale=2.5)

# tSNE plot
DimPlot(merged_SAB, reduction = "tsne", label = TRUE) + ggtitle("tSNE Clustering")
ggsave("tSNE_clustering_merged_data", width=1600, height=900, units="px", scale=2.5)

# Feature plot for a gene (e.g., "CD34")
FeaturePlot(merged_SAB, features = "CD34", reduction = "umap")
ggsave("CD34_feature_plot", width=1600, height=900, units="px", scale=2.5)

# Violin plot for gene expression across clusters
VlnPlot(merged_SAB, features = "CD34", group.by = "seurat_clusters")
ggsave("CD34_violin_plot", width=1600, height=900, units="px", scale=2.5)

# Ridge plot for gene expression
RidgePlot(merged_SAB, features = "CD34", group.by = "seurat_clusters")
ggsave("CD34_ridge_plot", width=1600, height=900, units="px", scale=2.5)

# Heatmap of top differentially expressed genes
top_markers <- FindAllMarkers(merged_SAB)
top10_markers <- top_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(merged_SAB, features = top10_markers$gene)

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(2410)
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

#message(paste("Running tSNE..."))
#merged_SAB <- RunTSNE(merged_SAB, dims = 1:30)
message(paste("Skipping tSNE..."))

# Set raster = FALSE for high-quality images
theme_set(theme_minimal())

# Define function to save plots with high resolution
save_high_res_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300, units = "in")
}

message("Generating UMAP plots...")
p1 <- DimPlot(merged_SAB, reduction = "umap", group.by = "Patient", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Patient")
save_high_res_plot(p1, "UMAP_clustering_merged_data_patient.jpeg")

p2 <- DimPlot(merged_SAB, reduction = "umap", group.by = "exp", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Experiment")
save_high_res_plot(p2, "UMAP_clustering_merged_data_exp.jpeg")

p3 <- DimPlot(merged_SAB, reduction = "umap", group.by = "day", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Day")
save_high_res_plot(p3, "UMAP_clustering_merged_data_day.jpeg")

p4 <- DimPlot(merged_SAB, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Predicted Cell Type") + NoLegend()
save_high_res_plot(p4, "UMAP_clustering_merged_data_celltype.jpeg")

#message("Generating tSNE plot...")
#p5 <- DimPlot(merged_SAB, reduction = "tsne", label = TRUE, raster = FALSE) + 
#  ggtitle("tSNE Clustering")
#save_high_res_plot(p5, "tSNE_clustering_merged_data.jpeg")

message("Switching to ADT assay and analyzing variable features...")
DefaultAssay(merged_SAB) <- "ADT"
merged_SAB <- FindVariableFeatures(merged_SAB, selection.method = "vst", nfeatures = 50, assay= "ADT")

# Extract top 10 variable features for ADT
top10_ADT <- head(VariableFeatures(merged_SAB), 10)
message("Top 10 ADT variable features: ", paste(top10_ADT, collapse = ", "))

p11 <- VariableFeaturePlot(merged_SAB) + ggtitle("Variable Features in ADT Assay")
p11 <- LabelPoints(plot = p11, points = top10_ADT, repel = TRUE)
save_high_res_plot(p11, "ADT_variable_features.jpeg")

message("Switching to RNA assay and analyzing variable features...")
DefaultAssay(merged_SAB) <- "RNA"
merged_SAB <- FindVariableFeatures(merged_SAB)

# Extract top 10 variable features for RNA
top10_RNA <- head(VariableFeatures(merged_SAB), 10)
message("Top 10 RNA variable features: ", paste(top10_RNA, collapse = ", "))

p12 <- VariableFeaturePlot(merged_SAB) + ggtitle("Variable Features in RNA Assay")
p12 <- LabelPoints(plot = p12, points = top10_RNA, repel = TRUE)
save_high_res_plot(p12, "RNA_variable_features.jpeg")

message("Generating stacked feature plots for selected markers...")
p13 <- FeaturePlot(merged_SAB, features = c("CD34", "CD3E", "MS4A1"), reduction = "umap", raster = FALSE) &
  theme(legend.position = "right")
save_high_res_plot(p13, "Stacked_feature_plot.jpeg")

message("Script completed successfully!")
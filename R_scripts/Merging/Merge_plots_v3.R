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
p1 <- DimPlot(merged_SAB, reduction = "umap", group.by = "Sorting", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Sort")
save_high_res_plot(p1, "UMAP_clustering_merged_data_Sort.jpeg")

message("Script completed successfully!")

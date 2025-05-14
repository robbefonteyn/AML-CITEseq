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

message("Generating UMAP plots...")
p1 <- DimPlot(merged_SAB, reduction = "umap", group.by = "Sorting", label = TRUE, raster = FALSE) + 
  ggtitle("UMAP Clustering by Sort")
save_high_res_plot(p1, "UMAP_clustering_merged_data_Sort.jpeg")

# Get unique values from the "exp" column
unique_exps <- unique(merged_SAB@meta.data$exp)

# Create and save a UMAP plot for each "exp" value
for (exp_value in unique_exps) {
  message(paste("Generating UMAP for experiment:", exp_value))
  
  # Identify the cells belonging to the current "exp" group
  cells_to_highlight <- rownames(merged_SAB@meta.data[merged_SAB@meta.data$exp == exp_value, ])
  
  # Generate UMAP plot with highlighted group and others in grey
  p_umap <- DimPlot(merged_SAB, reduction = "umap", group.by = "exp",
                    cells.highlight = cells_to_highlight, 
                    cols.highlight = "orange",  # Highlighted points in red
                    cols = "lightgray",      # Background points in grey
                    raster = FALSE) +
    ggtitle(paste("UMAP Clustering -", exp_value)) + NoLegend()
  
  # Save plot with high resolution
  filename <- paste0("UMAP_", exp_value, ".jpeg")
  ggsave(filename, plot = p_umap, width = 8, height = 6, dpi = 300, units = "in")
}

# Define sorting categories
sorting_categories <- unique(merged_SAB@meta.data$Sorting)

# Define colors
highlight_color <- "orange"
background_color <- "lightgray"

# Function to generate and save UMAP plots for a given metadata column
generate_umap_plots <- function(column_name, merged_SAB) {
  unique_values <- unique(merged_SAB@meta.data[[column_name]])
  
  for (value in unique_values) {
    message(paste("Generating UMAP for", column_name, ":", value))
    
    # Identify cells belonging to the current category
    cells_to_highlight <- rownames(merged_SAB@meta.data[merged_SAB@meta.data[[column_name]] == value, ])
    
    # Generate UMAP plot
    p_umap <- DimPlot(merged_SAB, reduction = "umap",
                      cells.highlight = cells_to_highlight, 
                      cols.highlight = highlight_color,  # Highlighted in light blue
                      cols = background_color,           # Background in grey
                      raster = FALSE) +
      ggtitle(paste("UMAP -", value, "Highlighted")) + 
      NoLegend()  # Removes legend and axes labels
    
    # Save plot with high resolution
    filename <- paste0("UMAP_", gsub(" ", "_", column_name), "_", gsub(" ", "_", value), ".jpeg")
    ggsave(filename, plot = p_umap, width = 8, height = 6, dpi = 300, units = "in")
  }
}

# Generate UMAP plots for Sorting
generate_umap_plots("Sorting", merged_SAB)

# Generate UMAP plots for Sample Type
generate_umap_plots("sample type", merged_SAB)

message(paste("Script completed successfully!"))


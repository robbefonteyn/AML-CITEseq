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

message("Script completed successfully!")
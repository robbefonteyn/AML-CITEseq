# Load libraries
library(Seurat)
library(scater)
library(ggvenn)
library(ggpubr)
library(dplyr)

# Create directory for plots
output_dir <- "output_figures_dietSeurat"
dir.create(output_dir, showWarnings = FALSE)

# Load Seurat object
message("Loading Seurat object...")
path <- "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_Full_with_Triana.rds"
SAB <- readRDS(path)

message("Creating slimmed Seurat object...")
SAB_slim <- DietSeurat(SAB, layers = c("counts", "data"), assays = c("RNA", "HTO", "ADT"))

DefaultAssay(SAB_slim) <- "RNA"
SAB_slim <- NormalizeData(SAB_slim)

# Histogram of raw library sizes
message("Plotting histogram of raw library sizes...")
p1 <- ggplot(SAB_slim@meta.data, aes(nCount_RNA)) + geom_histogram(bins=60)
ggsave(file.path(output_dir, "hist_raw_counts.png"), p1)

# Histogram of normalized library sizes
message("Plotting histogram of normalized library sizes...")
counts.norm <- as.matrix(SAB_slim@assays$RNA@data)
nCount_RNALN <- data.frame(nCount_RNA=colSums(counts.norm))
p2 <- ggplot(nCount_RNALN, aes(nCount_RNA)) + geom_histogram(bins=60)
ggsave(file.path(output_dir, "hist_normalized_counts.png"), p2)

# Variable feature selection and plotting
message("Finding variable features...")
SAB_slim <- FindVariableFeatures(SAB_slim)
p3 <- VariableFeaturePlot(SAB_slim)
ggsave(file.path(output_dir, "variable_features.png"), p3)

# Scaling data and PCA
message("Scaling data and running PCA...")
SAB_slim <- ScaleData(SAB_slim)
SAB_slim <- RunPCA(SAB_slim, nfeatures.print = 10)

# PCA plots
message("Plotting PCA...")
p4 <- DimPlot(SAB_slim, reduction="pca", group.by="Patient")
ggsave(file.path(output_dir, "pca_plot.png"), p4)

# PCA heatmaps
message("Generating PCA heatmaps...")
# Use Cairo for headless PNG generation
if (!requireNamespace("Cairo", quietly = TRUE)) {
  install.packages("Cairo", repos = "http://cran.us.r-project.org")
}
library(Cairo)

message("Generating PCA heatmaps using Cairo...")

CairoPNG(filename = file.path(output_dir, "pca_heatmap_1_12.png"), width = 1000, height = 1000)
DimHeatmap(SAB_slim, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

CairoPNG(filename = file.path(output_dir, "pca_heatmap_13_24.png"), width = 1000, height = 1000)
DimHeatmap(SAB_slim, dims = 13:24, cells = 500, balanced = TRUE)
dev.off()

CairoPNG(filename = file.path(output_dir, "pca_heatmap_25_36.png"), width = 1000, height = 1000)
DimHeatmap(SAB_slim, dims = 25:36, cells = 500, balanced = TRUE)
dev.off()


# Elbow plot
message("Creating ElbowPlot...")
p5 <- ElbowPlot(SAB_slim, ndims=50)
ggsave(file.path(output_dir, "elbow_plot.png"), p5)

# Clustering and dimensionality reduction
message("Running clustering and dimensionality reduction...")
SAB_slim <- FindNeighbors(SAB_slim, dims=1:30)
SAB_slim <- RunTSNE(SAB_slim, dims=1:30, check_duplicates=FALSE)
SAB_slim <- RunUMAP(SAB_slim, dims=1:30, n.neighbors=20)

# t-SNE and UMAP plots
message("Plotting t-SNE and UMAP...")
p6 <- DimPlot(SAB_slim, reduction="tsne", label=TRUE, label.size=8, pt.size=2) + NoLegend()
p7 <- DimPlot(SAB_slim, reduction="tsne", pt.size=2)
p8 <- DimPlot(SAB_slim)
p9 <- DimPlot(SAB_slim, label=TRUE, label.size=8) + NoLegend() + ggtitle("UMAP_on_PCA")

ggsave(file.path(output_dir, "tsne_labeled.png"), p6)
ggsave(file.path(output_dir, "tsne_plain.png"), p7)
ggsave(file.path(output_dir, "umap_basic.png"), p8)
ggsave(file.path(output_dir, "umap_labeled.png"), p9)

# Clustering ID plots
message("Plotting cluster identity comparisons...")
Idents(SAB_slim) <- SAB_slim@meta.data$RNA_snn_res.0.8
p10 <- DimPlot(SAB_slim, label=TRUE) + NoLegend() + NoAxes()
p11 <- DimPlot(SAB_slim, group.by = "Patient") + NoAxes()
p12 <- ggarrange(p10, p11)

ggsave(file.path(output_dir, "cluster_comparison.png"), p12)

# Save final slimmed Seurat object
message("Saving slimmed Seurat object...")
saveRDS(SAB_slim, file = "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_slimmed_final.rds")

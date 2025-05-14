# Set options early (important for HPC!)
options(timeout = 300)

# Load required libraries
message("Checking and installing required packages...")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of required packages
packages <- c(
  "ggplot2", 
  "Seurat", 
  "reshape2", 
  "plyr", 
  "DESeq2", 
  "ggrepel", 
  "SingleCellExperiment", 
  "scmap", 
  "parallel", 
  "randomForest"
)

# Install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing package: ", pkg))
    if (pkg %in% c("DESeq2", "SingleCellExperiment", "scmap")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  } else {
    message(paste0("Package ", pkg, " is already installed."))
  }
}

# Load all packages
message("Loading all packages...")
lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE))
  message(paste0("Loaded package: ", pkg))
})

message("Finished installing and loading packages.")

# Start actual processing
message("Reading Healthy dataset...")
Healthy <- readRDS(url("https://ndownloader.figshare.com/files/28408638"))

# Control flag
reproduce <- FALSE

message("Reading merged dataset...")
merged <- readRDS("C:/Users/Robbe Fonteyn/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/Processed_SCC/PROCESSED/SAB002/SAB002_reseq/SAB002_reseq_seurat.rds")

# Extract counts
message("Extracting assay data...")
normalized_counts <- GetAssayData(merged, assay = "RNA", slot = "data")
counts <- GetAssayData(merged, assay = "RNA", slot = "counts")

# Prepare projection matrix
message("Preparing counts projection...")
counts.projection <- normalized_counts[intersect(rownames(normalized_counts), rownames(Healthy)),]

# Set up query data
message("Setting up SingleCellExperiment objects...")
sce_Culture <- SingleCellExperiment(assays = list(normcounts = as.matrix(counts.projection)))
logcounts(sce_Culture) <- normcounts(sce_Culture)
rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)

sce_All <- SingleCellExperiment(assays = list(normcounts = as.matrix(Healthy@assays$BOTH@data[rownames(counts.projection),])))
logcounts(sce_All) <- normcounts(sce_All)
rowData(sce_All)$feature_symbol <- rownames(sce_All)
sce_All <- sce_All[!duplicated(rownames(sce_All)), ]

# Index cells
message("Indexing cells for scmap...")
sce_Culture <- setFeatures(sce_Culture, features = rownames(sce_Culture))
sce_Culture <- indexCell(sce_Culture)

sce_All <- setFeatures(sce_All, features = rownames(sce_All))
sce_All <- indexCell(sce_All)

# Index the reference clusters
sce_All <- indexCluster(sce_All)

# Perform scmapCluster projection (MUCH faster than scmapCell)
message("Running scmapCluster projection...")
Culture_Map <- scmapCluster(
  projection = sce_Culture,
  index_list = list(
    sce_All = metadata(sce_All)$scmap_cluster_index
  )
)

# Define Calc function
message("Defining correlation calculation function...")

Calc <- function(id, cult) {
  u <- cult[, id]
  xcoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,1][u]
  ycoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[,2][u]
  x = mean(xcoords)
  y = mean(ycoords)
  ct.t = table(Idents(Healthy)[u])
  ct.t <- ct.t[order(ct.t, decreasing = TRUE)]
  nearest <- Healthy@assays$BOTH@data[rownames(sce_Culture), u]
  query <- normcounts(sce_Culture)[, id]
  cornn <- apply(nearest, 2, function(x) cor(x, query))
  
  data.frame(row.names = id, x = x, y = y, ct = names(ct.t)[1], score = mean(cornn))
}

# Run mapping in parallel
message("Running parallel mapping... (using 6 cores)")

# HPC Note: mc.cores should ideally be dynamically set
n_cores <- 6
mapped <- mclapply(colnames(Culture_Map$sce_All[[1]]), Calc, cult = Culture_Map$sce_All[[1]])

# Combine results
mapped <- do.call(rbind, mapped)

# Update metadata
message("Updating merged object metadata with projected cell types...")

# Add Triana_Celltypes to metadata
merged$Triana_Celltypes <- mapped$ct

# Save the newly annotated merged object
message("Saving the newly annotated merged object...")
saveRDS(merged, file = "merged_annotated.rds")

message("Annotated merged object saved as 'merged_annotated.rds'")

message("Pipeline finished successfully!")



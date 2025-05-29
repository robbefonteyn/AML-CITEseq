# Set options
options(timeout = 300)

# Load required libraries
message("Checking and installing required packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("ggplot2", "Seurat", "reshape2", "plyr", "DESeq2", "ggrepel", "SingleCellExperiment", "scmap", "parallel", "randomForest")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("DESeq2", "SingleCellExperiment", "scmap")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}
lapply(packages, function(pkg) suppressMessages(library(pkg, character.only = TRUE)))
message("Finished installing and loading packages.")

# Define directory
rds_directory <- "C:/Users/Robbe Fonteyn/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/Processed_SCC/PROCESSED"
output_file <- "Triana_celltype_assignments.csv"

# Load Triana Healthy reference
message("Loading Triana Healthy reference...")
Healthy <- readRDS(url("https://ndownloader.figshare.com/files/28408638"))

# Prepare reference
sce_All <- SingleCellExperiment(assays = list(normcounts = as.matrix(Healthy@assays$BOTH@data)))
logcounts(sce_All) <- normcounts(sce_All)
rowData(sce_All)$feature_symbol <- rownames(sce_All)
sce_All <- sce_All[!duplicated(rownames(sce_All)),]
sce_All <- setFeatures(sce_All, features = rownames(sce_All))
sce_All <- indexCell(sce_All)

# List RDS files
rds_files <- list.files(rds_directory, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

# Function to project and assign Triana celltypes
process_rds <- function(file_path) {
  message(paste("Processing", file_path))
  merged <- readRDS(file_path)
  
  # Normalize RNA
  DefaultAssay(merged) <- "RNA"
  merged <- NormalizeData(merged, assay = "RNA")
  merged <- FindVariableFeatures(merged, assay = "RNA")
  merged <- ScaleData(merged, assay = "RNA")
  merged <- RunPCA(merged, assay = "RNA", reduction.name = "pca")
  rna_data <- GetAssayData(merged, assay = "RNA", layer = "data")
  
  # Normalize ADT
  DefaultAssay(merged) <- "ADT"
  merged <- NormalizeData(merged, assay = "ADT", normalization.method = "CLR")
  merged <- ScaleData(merged, assay = "ADT")
  merged <- RunPCA(merged, assay = "ADT", reduction.name = "apca")
  adt_data <- GetAssayData(merged, assay = "ADT", layer = "data")
  
  # Multimodal integration
  merged <- FindMultiModalNeighbors(merged, reduction.list = list("pca", "apca"), dims.list = list(1:30, 1:18))
  merged <- RunUMAP(merged, nn.name = "weighted.nn", reduction.name = "wnn.umap")
  
  # Combine RNA and ADT into BOTH assay
  combined <- rbind(rna_data, adt_data)
  merged[["BOTH"]] <- CreateAssayObject(data = combined)
  
  # Use BOTH assay for projection
  DefaultAssay(merged) <- "BOTH"
  normalized_counts <- GetAssayData(merged, assay = "BOTH", layer = "data")
  common_genes <- intersect(rownames(normalized_counts), rownames(sce_All))
  counts.projection <- normalized_counts[common_genes,]
  
  sce_Culture <- SingleCellExperiment(assays = list(normcounts = as.matrix(counts.projection)))
  logcounts(sce_Culture) <- normcounts(sce_Culture)
  rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)
  sce_Culture <- setFeatures(sce_Culture, features = rownames(sce_Culture))
  sce_Culture <- indexCell(sce_Culture)
  
  Culture_Map <- scmapCell(
    projection = sce_Culture,
    index_list = list(
      sce_All = metadata(sce_All)$scmap_cell_index
    ),
    w = 5
  )
  
  Calc <- function(id, cult) {
    u <- cult[, id]
    xcoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[, 1][u]
    ycoords <- Healthy@reductions$MOFAUMAP@cell.embeddings[, 2][u]
    x = mean(xcoords)
    y = mean(ycoords)
    ct.t = table(Idents(Healthy)[u])
    ct.t <- ct.t[order(ct.t, decreasing = TRUE)]
    nearest <- Healthy@assays$BOTH@data[rownames(sce_Culture), u]
    query <- normcounts(sce_Culture)[, id]
    cornn <- apply(nearest, 2, function(x) cor(x, query))
    data.frame(row.names = id, Triana_celltype = names(ct.t)[1], score = mean(cornn))
  }
  
  mapped <- lapply(colnames(Culture_Map$sce_All[[1]]), Calc, cult = Culture_Map$sce_All[[1]])
  mapped <- do.call(rbind, mapped)
  
  result <- data.frame(
    Barcode = rownames(mapped),
    Triana_celltype = mapped$Triana_celltype,
    Original_identity = basename(file_path)
  )
  rownames(result) <- result$Barcode
  return(result)
}

# Run the pipeline
message("Starting batch processing...")
all_results <- lapply(rds_files, process_rds)
final_table <- do.call(rbind, all_results)

# Save output
write.csv(final_table, output_file, row.names = TRUE)
message("All done! Final table saved to:", output_file)

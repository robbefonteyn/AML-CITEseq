# Load necessary libraries
library(Seurat)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(presto)

# Define function to save plots with high resolution
save_high_res_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300, units = "in")
}

message(paste("Loading merged object..."))
path <- "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_Full_with_meta.rds"
SAB <- readRDS(path)

DefaultAssay(SAB)  <- "RNA"

message(paste("Subsetting Seurat object by CD34+"))
# Subset the Seurat object to only include cells with "CD34+" in the 'Sorting' column
SAB_CD34_AML_HC <- subset(SAB, subset = Sorting == "CD34+")

message(paste("Subsetting by AML and HC"))
unique_patients <- unique(SAB_CD34_AML_HC@meta.data$Patient)
aml_patients <- unique_patients[grepl("^AML\\d+", unique_patients)]
hc_patients <- unique_patients[grepl("^HC\\d+", unique_patients)]

message(paste("Normalizing subset"))
# Normalize the data
SAB_CD34_AML_HC <- NormalizeData(SAB_CD34_AML_HC)

message(paste("Identities setting for DE to patient column"))
# Set the identities for differential expression to the Patient column
Idents(SAB_CD34_AML_HC) <- "Patient"

message(paste("Finding markers"))
# Perform differential expression analysis between AML and HC
markers_AML_vs_HC <- FindMarkers(SAB_CD34_AML_HC, 
                                 ident.1 = aml_patients,
                                 ident.2 = hc_patients,
                                 test.use = "wilcox_limma",
                                 verbose = TRUE)

# Sort results by adjusted p-value
markers_AML_vs_HC <- markers_AML_vs_HC %>% arrange(p_val_adj)

# Save the DE results to a CSV file
write.csv(markers_AML_vs_HC, file = "DE_AML_vs_HC_CD34_wilcox_limma.csv")

message(paste("Generating Volcano Plot"))
# Generate a Volcano Plot
markers_AML_vs_HC$gene <- rownames(markers_AML_vs_HC)
volcano_plot <- ggplot(markers_AML_vs_HC, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes (AML vs HC)")

save_high_res_plot(volcano_plot, "Volcano_Plot_DE_AML_vs_HC_wilcox_limma.jpeg")

message(paste("Generating and saving individual violin plots"))
# Generate and save individual violin plots for each of the top 5 DE genes
top_genes <- rownames(markers_AML_vs_HC)[1:5]

for (gene in top_genes) {
  violin_plot <- VlnPlot(SAB_CD34_AML_HC, features = gene, group.by = "Patient", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  
  save_high_res_plot(violin_plot, paste0("Violin_", gene, ".jpeg"))
}

message(paste("Script successfully finished!"))

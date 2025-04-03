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
path <- "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_Full_with_patient.rds"
SAB <- readRDS(path)

DefaultAssay(SAB)  <- "RNA"

message(paste("Subsetting Seurat object by CD34+"))
# Subset the Seurat object to only include cells with "CD34+" in the 'Sorting' column
SAB_CD34 <- subset(SAB, subset = Sorting == "CD34+")

message(paste("Subsetting to cells with non-NA values in Genetics"))
# Subset to cells where Genetics is not NA
SAB_CD34_gen <- subset(SAB_CD34, subset = !is.na(Genetics))


message(paste("Normalizing subset"))
# Normalize the data
SAB_CD34_gen <- NormalizeData(SAB_CD34_gen)

message(paste("Setting identities to Genetics"))
# Set the identities for differential expression to the Genetics column
Idents(SAB_CD34_gen) <- "Genetics"

message(paste("Finding markers between 'mutated NPM1' and all others"))
# Perform differential expression analysis between 'mutated NPM1' and all other identities
markers_npm1 <- FindMarkers(SAB_CD34_gen, 
                            ident.1 = "mutated NPM1",
                            ident.2 = NULL,  # Compare against all others
                            test.use = "wilcox_limma",
                            verbose = TRUE)

# Sort results by adjusted p-value
markers_npm1 <- markers_npm1 %>% arrange(p_val_adj)

# Save the DE results to a CSV file
write.csv(markers_npm1, file = "DE_Genetics_MutatedNPM1_vs_Others_wilcox_limma.csv")

message(paste("Generating Volcano Plot"))
# Generate a Volcano Plot
markers_npm1$gene <- rownames(markers_npm1)
volcano_plot_npm1 <- ggplot(markers_npm1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of DE Genes (mutated NPM1 vs Others)")

save_high_res_plot(volcano_plot_npm1, "Volcano_Plot_MutatedNPM1_vs_Others.jpeg")

message(paste("Generating and saving individual violin plots"))
# Generate and save violin plots for the top 5 DE genes
top_genes_npm1 <- rownames(markers_npm1)[1:5]

for (gene in top_genes_npm1) {
  violin_plot_npm1 <- VlnPlot(SAB_CD34_gen, features = gene, group.by = "Genetics", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  
  save_high_res_plot(violin_plot_npm1, paste0("Violin_", gene, "_MutatedNPM1.jpeg"))
}

message(paste("Script successfully finished!"))

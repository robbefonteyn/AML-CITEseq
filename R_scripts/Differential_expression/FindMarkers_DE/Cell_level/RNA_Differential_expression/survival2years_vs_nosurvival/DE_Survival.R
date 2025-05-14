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

message(paste("Subsetting to cells with 'yes' or 'no' in Survival_two_years"))
# Subset to cells where Survival_two_years is "yes" or "no"
SAB_CD34_surv <- subset(SAB_CD34, subset = Survival_two_years %in% c("yes", "no"))

message(paste("Normalizing subset"))
# Normalize the data
SAB_CD34_surv <- NormalizeData(SAB_CD34_surv)

message(paste("Setting identities to Survival_two_years"))
# Set the identities for differential expression to the Survival_two_years column
Idents(SAB_CD34_surv) <- "Survival_two_years"
write.csv(SAB_CD34_surv@meta.data, file = "SAB_CD34_survival_subset_metadata.csv")

message(paste("Finding markers between Survival 'yes' and 'no'"))
# Perform differential expression analysis between Survival_two_years == "yes" and "no"
markers_surv <- FindMarkers(SAB_CD34_surv, 
                            ident.1 = "yes",
                            ident.2 = "no",
                            test.use = "wilcox_limma",
                            verbose = TRUE)

# Sort results by adjusted p-value
markers_surv <- markers_surv %>% arrange(p_val_adj)

# Save the DE results to a CSV file
write.csv(markers_surv, file = "DE_Survival_CD34_yes_vs_no_wilcox_limma.csv")

message(paste("Generating Volcano Plot"))
# Generate a Volcano Plot
markers_surv$gene <- rownames(markers_surv)
volcano_plot <- ggplot(markers_surv, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of DE Genes (Survival: Yes vs No)")

save_high_res_plot(volcano_plot, "Volcano_Plot_Survival_CD34_yes_vs_no.jpeg")

message(paste("Generating and saving individual violin plots"))
# Generate and save violin plots for the top 5 DE genes
top_genes <- rownames(markers_surv)[1:5]

for (gene in top_genes) {
  violin_plot <- VlnPlot(SAB_CD34_surv, features = gene, group.by = "Survival_two_years", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  
  save_high_res_plot(violin_plot, paste0("Violin_", gene, ".jpeg"))
}

message(paste("Script successfully finished!"))
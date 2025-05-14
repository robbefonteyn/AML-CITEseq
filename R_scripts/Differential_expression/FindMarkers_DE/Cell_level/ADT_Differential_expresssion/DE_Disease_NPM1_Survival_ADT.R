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

message("Loading merged object...")
path <- "/kyukon/home/gent/458/vsc45888/data/Master_dissertation/Processed_SCC/Merge_test/SAB_Full_with_patient.rds"
SAB <- readRDS(path)

DefaultAssay(SAB) <- "ADT"

message("Subsetting Seurat object by CD34+")
SAB_CD34 <- subset(SAB, subset = Sorting == "CD34+")

## 1. DE by Survival (Yes vs No)
message("Running DE analysis for Survival status (Yes vs No)...")
SAB_CD34_surv <- subset(SAB_CD34, subset = Survival_two_years %in% c("yes", "no"))
SAB_CD34_surv <- NormalizeData(SAB_CD34_surv)
Idents(SAB_CD34_surv) <- "Survival_two_years"

markers_surv <- FindMarkers(SAB_CD34_surv, ident.1 = "yes", ident.2 = "no", test.use = "wilcox_limma", verbose = TRUE)
markers_surv <- markers_surv %>% arrange(p_val_adj)
write.csv(markers_surv, file = "DE_Survival_CD34_ADT_yes_vs_no.csv")

markers_surv$gene <- rownames(markers_surv)
volcano_plot_surv <- ggplot(markers_surv, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of DE ADT Markers (Survival Yes vs No)")
save_high_res_plot(volcano_plot_surv, "Volcano_Survival_CD34_ADT_yes_vs_no.jpeg")

top_genes_surv <- rownames(markers_surv)[1:5]
for (gene in top_genes_surv) {
  vplot <- VlnPlot(SAB_CD34_surv, features = gene, group.by = "Survival_two_years", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  save_high_res_plot(vplot, paste0("Violin_Survival_", gene, "_ADT.jpeg"))
}

## 2. DE by AML vs HC
message("Running DE analysis for AML vs HC...")
unique_patients <- unique(SAB_CD34@meta.data$Patient)
aml_patients <- unique_patients[grepl("^AML\\d+", unique_patients)]
hc_patients <- unique_patients[grepl("^HC\\d+", unique_patients)]

SAB_CD34_AML_HC <- NormalizeData(SAB_CD34)
Idents(SAB_CD34_AML_HC) <- "Patient"
markers_AML_vs_HC <- FindMarkers(SAB_CD34_AML_HC, ident.1 = aml_patients, ident.2 = hc_patients, test.use = "wilcox_limma", verbose = TRUE)
markers_AML_vs_HC <- markers_AML_vs_HC %>% arrange(p_val_adj)
write.csv(markers_AML_vs_HC, file = "DE_AML_vs_HC_CD34_ADT.csv")

markers_AML_vs_HC$gene <- rownames(markers_AML_vs_HC)
volcano_plot_aml_hc <- ggplot(markers_AML_vs_HC, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of DE ADT markers (AML vs HC)")
save_high_res_plot(volcano_plot_aml_hc, "Volcano_Plot_DE_AML_vs_HC_ADT.jpeg")

top_genes_aml_hc <- rownames(markers_AML_vs_HC)[1:5]
for (gene in top_genes_aml_hc) {
  vplot <- VlnPlot(SAB_CD34_AML_HC, features = gene, group.by = "Patient", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  save_high_res_plot(vplot, paste0("Violin_AML_HC_", gene, "_ADT.jpeg"))
}

## 3. DE by NPM1 Mutation Status
message("Running DE analysis for NPM1 mutation status...")
message("Subsetting to non-NA Genetics")
non_na_genetics <- !is.na(SAB_CD34@meta.data$Genetics)
SAB_CD34_gen <- subset(SAB_CD34, cells = Cells(SAB_CD34)[non_na_genetics])

SAB_CD34_gen <- NormalizeData(SAB_CD34_gen)
Idents(SAB_CD34_gen) <- "Genetics"

markers_npm1 <- FindMarkers(SAB_CD34_gen, ident.1 = "mutated NPM1", ident.2 = NULL, test.use = "wilcox_limma", verbose = TRUE)
markers_npm1 <- markers_npm1 %>% arrange(p_val_adj)
write.csv(markers_npm1, file = "DE_Genetics_MutatedNPM1_vs_Others_ADT.csv")

markers_npm1$gene <- rownames(markers_npm1)
volcano_plot_npm1 <- ggplot(markers_npm1, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot of DE ADT markers (mutated NPM1 vs Others)")
save_high_res_plot(volcano_plot_npm1, "Volcano_MutatedNPM1_vs_Others_ADT.jpeg")

top_genes_npm1 <- rownames(markers_npm1)[1:5]
for (gene in top_genes_npm1) {
  vplot <- VlnPlot(SAB_CD34_gen, features = gene, group.by = "Genetics", pt.size = 0.1) +
    ggtitle(paste("Violin Plot for", gene))
  save_high_res_plot(vplot, paste0("Violin_MutatedNPM1_", gene, "_ADT.jpeg"))
}

message("Unified DE script using ADT assay completed.")

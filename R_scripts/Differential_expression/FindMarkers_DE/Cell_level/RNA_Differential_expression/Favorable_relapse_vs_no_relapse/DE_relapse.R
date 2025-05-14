# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Define function to save high resolution plots
save_high_res_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300, units = "in")
}

# Load Seurat object
message("Loading merged object...")
path <- "SAB_Full_with_patient.rds"
SAB <- readRDS(path)
DefaultAssay(SAB)  <- "RNA"

# Subset to L/D
message("Subsetting Seurat object by L/D")
SAB_LD <- subset(SAB, subset = Sorting == "L/D")
rm(SAB)

# Subset to cells with non-NA ELN and Relapse
message("Subsetting to non-NA ELN and Relapse")
non_na_eln_relapse <- !is.na(SAB_LD@meta.data$ELN) & !is.na(SAB_LD@meta.data$Relapse)
SAB_LD <- subset(SAB_LD, cells = Cells(SAB_LD)[non_na_eln_relapse])

# Define relapse labels
relapse_labels <- c("after chemo", "prim refr")

# Subset groups
message("Defining patient groups for contrast")
fav_no_relapse_cells <- rownames(SAB_LD@meta.data[
  SAB_LD@meta.data$ELN == "favorable" & !(SAB_LD@meta.data$Relapse %in% relapse_labels),
])
SAB_LD_fav_no_relapse <- subset(SAB_LD, cells = fav_no_relapse_cells)

favint_relapse_cells <- rownames(SAB_LD@meta.data[
  SAB_LD@meta.data$ELN %in% c("favorable", "intermediate") & SAB_LD@meta.data$Relapse %in% relapse_labels,
])
SAB_LD_favint_relapse <- subset(SAB_LD, cells = favint_relapse_cells)

# Merge subsets and define group
SAB_DE_contrast <- merge(SAB_LD_fav_no_relapse, y = SAB_LD_favint_relapse)

# Create contrast group label
SAB_DE_contrast$Group <- ifelse(
  SAB_DE_contrast@meta.data$ELN == "favorable" & !(SAB_DE_contrast@meta.data$Relapse %in% relapse_labels),
  "favorable_no_relapse",
  "favint_with_relapse"
)

# Assign identity for DE analysis
Idents(SAB_DE_contrast) <- "Group"

# Normalize and run DE
message("Normalizing and running differential expression")
SAB_DE_contrast <- NormalizeData(SAB_DE_contrast)
markers_relapse_contrast <- FindMarkers(SAB_DE_contrast,
                                        ident.1 = "favint_with_relapse",
                                        ident.2 = "favorable_no_relapse",
                                        test.use = "wilcox",
                                        verbose = TRUE)

# Sort and save results
markers_relapse_contrast <- markers_relapse_contrast %>% arrange(p_val_adj)
write.csv(markers_relapse_contrast, "DE_FavIntRelapse_vs_FavNoRelapse.csv")

# Volcano plot
message("Generating volcano plot")
markers_relapse_contrast$gene <- rownames(markers_relapse_contrast)
volcano_plot <- ggplot(markers_relapse_contrast, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = p_val_adj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("DE: Fav/Int Risk Relapse vs Fav No Relapse")

save_high_res_plot(volcano_plot, "Volcano_FavIntRelapse_vs_FavNoRelapse.jpeg")

message("Script successfully finished!")

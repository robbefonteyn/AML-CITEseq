setwd("C:/Users/Robbe Fonteyn/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/Annotation")
meta <- read.csv("Metadata_with_seurat_cellNr.csv", header = TRUE, sep = ",")

library(ggplot2)
library(scales)  # For rescaling function

# Define custom blue scale colors
color_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                "#a55194", "#393b79", "#8c6d31", "#e7ba52", "#843c39",
                "#d6616b", "#7b4173", "#fdae6b", "#31a354", "#636363")

#library(RColorBrewer)
#color_palette <- brewer.pal(n = 20, name = "Set3")  # or use "Set3", "Spectral", "Dark2"
#scale_fill_manual(values = color_palette)

library(viridis)
#color_palette <- viridis(20, option = "E")  # Options: "A", "B", "C", "D", "E"
scale_fill_manual(values = color_palette)

# Create the plot with improved label placement
ggplot(meta, aes(x = exp, y = Nr.cells.sorted, fill = Patient)) + 
  geom_bar(position = "fill", stat = "identity", color = "grey") + 
  
  # Labels: Vertical for most, Horizontal for small bars
  geom_text(aes(label = Nr.cells.sorted, 
                angle = ifelse(Nr.cells.sorted < quantile(Nr.cells.sorted, 0.1), 0, 90)),  
            position = position_fill(vjust = 0.5),  # Keeps text centered within bars
            size = 3, 
            color = "white", 
            check_overlap = TRUE) +  # Prevents excessive overlap
  theme_bw() +
  facet_wrap(~ Sorting, scales = "free_x")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_palette)
# Save plot as high-resolution PNG
ggsave("dataset_overview_sorted_cells.png", width = 10, height = 6, dpi = 300)

# Create the plot shaded by sample.type
ggplot(meta, aes(x = exp, y = Nr.cells.sorted, fill = sample.type)) + 
  geom_bar(position = "fill", stat = "identity", color = "grey") + 
  
  # Labels: Vertical for most, Horizontal for small bars
  geom_text(aes(label = Nr.cells.sorted, 
                angle = ifelse(Nr.cells.sorted < quantile(Nr.cells.sorted, 0.1), 0, 90)),  
            position = position_fill(vjust = 0.5),
            size = 3, 
            color = "white", 
            check_overlap = TRUE) +
  
  theme_bw() +
  facet_wrap(~ Sorting, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) #+
  
  # Optional custom colors for two sample types
  #scale_fill_manual(values = c("BM" = "#E41A1C", "PB" = "#377EB8", "AF" = "orange"))  # adjust names to match your values

# Save plot as high-resolution PNG
ggsave("dataset_overview_sorted_cells_by_sampletype.png", width = 10, height = 6, dpi = 300)


# Create the plot with improved label placement
ggplot(meta, aes(x = exp, y = Seurat_cells, fill = Patient)) + 
  geom_bar(position = "fill", stat = "identity", color = "grey") + 
  
  # Labels: Vertical for most, Horizontal for small bars
  geom_text(aes(label = Seurat_cells, 
                angle = ifelse(Seurat_cells < quantile(Seurat_cells, 0.1), 0, 90)),  
            position = position_fill(vjust = 0.5),  # Keeps text centered within bars
            size = 3, 
            color = "white", 
            check_overlap = TRUE) +  # Prevents excessive overlap
  theme_bw() + 
  facet_wrap(~ Sorting, scales = "free_x")  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_palette)
# Save plot as high-resolution PNG
ggsave("dataset_overview_seurat_cells.png", width = 10, height = 6, dpi = 300)


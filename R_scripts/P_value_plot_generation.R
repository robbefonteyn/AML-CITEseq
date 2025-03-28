DE <- read.csv("DE_Survival_CD34_yes_vs_no_wilcox_limma.csv", header = TRUE, sep = ",", row.names = 1)

# Move rownames into a column called "gene"
DE$gene <- rownames(DE)

# Now plot using ggplot2
library(ggplot2)

# Convert adjusted p-values to numeric (if needed)
DE$p_val_adj <- as.numeric(DE$p_val_adj)

# Count significant genes
n_signif <- sum(DE$p_val_adj < 0.05, na.rm = TRUE)

# Plot histogram with vertical line at 0.05
library(ggplot2)
ggplot(DE, aes(x = p_val_adj)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  xlab("Adjusted p-value") +
  ylab("Count") +
  ggtitle(paste0("Distribution of adj P values 2 year survival\nSignificant genes (padj < 0.05): ", n_signif))


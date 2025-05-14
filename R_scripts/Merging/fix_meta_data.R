library(Seurat)

SAB <- readRDS(file= "C:/Users/fonte/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/SAB_Full.rds")

meta <- SAB@meta.data
patientinfo <- read.csv2("C:/Users/fonte/OneDrive - UGent/AJ 2024-2025/Master's_Dissertation/AML-CITEseq/Data/Metadata_table_clean.csv")

# Add empty columns with NA values to meta table
meta$day <- NA
meta$exp <- NA
meta$patient <- NA
meta$sorting <- NA
meta$sample.ID <- NA
meta$sample.type <- NA
meta$sample.nr <- NA
meta$X..cells.sorted <- NA

# change the meta table to have the right columns
for (i in 1:nrow(meta)){
  
}

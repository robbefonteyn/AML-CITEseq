library(Seurat)
library(dplyr)
library(remotes)
library(SeuratObject)

path <- "SAB001_annotated.rds"
SAB <- readRDS(path)
message(paste("Loaded SAB001 to merge onto, from", path))


path <- "Processed_SCC/"
seurat_files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)

message(paste("Reading in RDS files from path:", path))
seurat_list <- lapply(seurat_files, readRDS)

message(paste("Done reading in files, now the vector is merged"))

sab_vector <- c("SAB001", "SAB002", "SAB003", "SAB004", "SAB005", "SAB006", "SAB007", "SAB008", "SAB009", "SAB010", "SAB011", "SAB012", "SAB013", "SAB014", "SAB015", "SAB016", "SAB017", "SAB018", "SAB019", "SAB020", "SAB021", "SAB022", "SAB023", "SAB024", "SAB025", "SAB026", "SAB027", "SAB028", "SAB029", "SAB030", "SAB031", "SAB032", "SAB033", "SAB035")

message(paste("Merging all the files together now, this will take some time"))
SAB <- merge(x= SAB, y= seurat_list, add.cell.ids = sab_vector, merge.data = TRUE, merge.dr=NA, project="AML")

message(paste("Saving new merged object:"))
saveRDS(SAB, "SAB_merged.rds")

message("Done with merging and saving")

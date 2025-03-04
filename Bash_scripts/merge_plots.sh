#!/usr/bin/env bash

#PBS -N merging_plots_script
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=350gb  

module load R/4.4.1-gfbf-2023b
#module load Seurat/5.1.0

cd /kyukon/home/gent/458/vsc45888/Master_dissertation/R_scripts/Merging

# Install missing dependencies if not available
# Rscript -e 'if (!requireNamespace("listenv", quietly = TRUE)) install.packages("listenv", repos="http://cran.us.r-project.org")'
# Rscript -e 'if (!requireNamespace("globals", quietly = TRUE)) install.packages("globals", repos="http://cran.us.r-project.org")'
# Rscript -e 'if (!requireNamespace("SeuratObject", quietly = TRUE)) install.packages("SeuratObject", repos="http://cran.us.r-project.org")'
# Rscript -e 'if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat", repos="http://cran.us.r-project.org")'

echo "Loading R script"
Rscript Merge_plots_v4.R

echo "merging_script finished"
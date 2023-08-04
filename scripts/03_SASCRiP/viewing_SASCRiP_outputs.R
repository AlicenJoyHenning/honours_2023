# VIEWING SASCRiP OUTPUTS

install.packages("readr")
library(readr)
install.packages("Matrix")
library(Matrix)

# viewing the seurat_matrix outputs of 2023 and 2020 index 

setwd("Alicen/work_files/ifnalpha/seurat_matrix/")

alpha_matrix_23 <- readMM("matrix.mtx.gz")

setwd("~/Alicen/work_files/tests/ifnalpha_darisia_index/seurat_matrix/")

alpha_matrix_20 <- readMM("matrix.mtx.gz")

# viewing the run_cqc outputs of 2023 and 2020  index 

setwd("~/Alicen/work_files/ifnalpha/run_cqc/")
alpha_run_cqc_2023 <- readRDS("ifnalpha_subset_seurat.rds")
alpha_run_cqc_2023
View(alpha_run_cqc_2023)

setwd("~/Alicen/work_files/tests/ifnalpha_darisia_index/run_cqc/")
alpha_run_cqc_2020 <- readRDS("ifnalpha_subset_seurat.rds")
alpha_run_cqc_2020 
View(alpha_run_cqc_2020)

# viewing the sctransform_normalize outputs 

setwd("~/Alicen/work_files/ifnalpha/sctransform_normalize/")
alpha_sctransform_normalize_2023 <- readRDS("alpha_test_normalised_seurat.rds")
View(alpha_sctransform_normalize_2023)
head(alpha_sctransform_normalize_2023)

setwd("~/Alicen/work_files/tests/ifnalpha_darisia_index/sctransform_normalize/")
alpha_sctransform_normalize_2020 <- readRDS("alpha_test_normalised_seurat.rds")
View(alpha_sctransform_normalize_2020)

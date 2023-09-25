# Quantifying PBMC populations across treatment 

##### [1] Load annotated dataset : ####
TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
DefaultAssay(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated) <- "RNA"

##### [2] View and store PBMC population names : ####
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated)
levels(TreatmentAnnotated)

##### [3] Quantify using for loop : ####
# Define the list of cell types
cell_types <- c(
  "monocytes", "CD4+_helper", "neutrophils", "T", "naive_CD8+ T",
  "B", "NKT", "mDCs", "cytotoxic_CD8+_T", "Tregs",
  "NK", "platelets", "unknown", "DCs", "CD4+_T", "pDCs"
)

# Initialize an empty list to store the results
results <- list()

# Loop through each cell type
for (cell_type in cell_types) {
  # Create a subset for each treatment and cell type combination
  T_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type)
  A_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "alpha")
  L_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "lambda")
  U_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "untreated")
  
  # Quantify the number of cells in each subset
  num_cells_T <- dim(T_subset)[2]
  num_cells_A <- dim(A_subset)[2]
  num_cells_L <- dim(L_subset)[2]
  num_cells_U <- dim(U_subset)[2]
  
  # Store the results in a list
  cat("\n Cell Type | ", cell_type, "\n", "Total:", num_cells_A)
  cat("\n Alpha :", num_cells_A)
  cat("\n Lambda :", num_cells_L)
  cat("\n Untreated :", num_cells_U, "\n")
}




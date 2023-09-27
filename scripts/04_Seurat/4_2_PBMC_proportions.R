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
  cat("\n Cell Type | ", cell_type, "\n", "Total:", num_cells_T)
  cat("\n Alpha :", num_cells_A)
  cat("\n Lambda :", num_cells_L)
  cat("\n Untreated :", num_cells_U, "\n")

}


##### [4] Alterations to the loop to account for errors when no subsetting results are found #####

# Loop through each PMBC pop to quantify cells in the pop per treatment type 

for (cell_type in cell_types) {
  # Create a subset for each treatment and cell type combination
  T_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type)
  
  # Initialize variables to store cell # for each treatment
  num_cells_A <- 0
  num_cells_L <- 0
  num_cells_U <- 0
  
  # Wrap each subseting calculation in a tryCatch block to prevent premature
  # loop exits when an 'error' is outputted from one subsetting result 
  
  # Alpha treatment 
  tryCatch({
    A_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "alpha") # create subset 
    num_cells_A <- dim(A_subset)[2] # view dimensions of resulting Seurat object (= # cells)
  }, error = function(e) { 
    num_cells_A <- 0  # Set num_cells_A to 0 if there was an error
  })
  
  # Lambda treatment 
  tryCatch({
    L_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "lambda")
    num_cells_L <- dim(L_subset)[2]
  }, error = function(e) {
    num_cells_L <- 0  
  })
  
  # Untreated control 
  tryCatch({
    U_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & treatment == "untreated")
    num_cells_U <- dim(U_subset)[2]
  }, error = function(e) {
    num_cells_U <- 0  
  })
  
  # Store the results in a list
  cat("\nCell Type | ", cell_type, "\n")  # \n to print results nicely on a separate line 
  cat("Total:", dim(T_subset)[2], "\n")
  cat("Alpha :", ifelse(num_cells_A > 0, num_cells_A, "No data"), "\n")
  cat("Lambda :", ifelse(num_cells_L > 0, num_cells_L, "No data"), "\n")
  cat("Untreated :", ifelse(num_cells_U > 0, num_cells_U, "No data"), "\n") 
}











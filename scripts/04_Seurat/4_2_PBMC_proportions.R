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
  "monocytes", "neutrophils", "DCs", "mDCs", "pDCs", "platelets", "T helper", "naive CD4 T", "naive CD8 T", "cytotoxic T", "NKT", "Tregs", "Tcm", "NK", "B"
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


##### [4] Alterations to the loop to account for errors when no subseting results are found #####
# Define the list of cell types
cell_types <- c(
  "monocytes", "neutrophils", "DCs", "mDCs", "pDCs", "platelets", "T helper", "naive CD4 T", "naive CD8 T", "cytotoxic T", "NKT", "Tregs", "Tcm", "NK", "B"
)

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
    A_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & sample == "alpha") # create subset 
    num_cells_A <- dim(A_subset)[2] # view dimensions of resulting Seurat object (= # cells)
  }, error = function(e) { 
    num_cells_A <- 0  # Set num_cells_A to 0 if there was an error
  })
  
  # Lambda treatment 
  tryCatch({
    L_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & sample == "lambda")
    num_cells_L <- dim(L_subset)[2]
  }, error = function(e) {
    num_cells_L <- 0  
  })
  
  # Untreated control 
  tryCatch({
    U_subset <- subset(TreatmentAnnotated, subset = celltype == cell_type & sample == "untreated")
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


##### [4] Present in graph #####
# Y values : 
proportion <- c(1674, 929, 680, # monocytes
          766, 618, 415, # neutrophils
          # 1209, 13, # myeloid 
          52, 21, 15,   # DCs
          388, 8, 9, # myeloid dericed dendritic cells
          17, 0, 0, # pDCs
          
          # 297, 51, # dendritic cells 
          
          101, 41, 27, # platelets 
          
          718, 326, 463, # T helper
          759, 613, 668, # naive CD4
          454, 126, 191, # naive cd8
          259, 71, 79, # cyto CD8
          259, 37, 38, # NKT
          
          112, 90, 91, # T regs 
          44, 24, 34, # Tcm
          185, 34, 42, # NK
          # 935, 33, # overall T 
          
          347, 54, 82 # B
)

# x groups : 
CellTypesProp <- c(rep("monocytes", 3), 
               rep("neutrophils", 3),
               # rep("myeloid", 2),
               rep("DCs", 3),
               rep("mDCs", 3),
               rep("pDCs", 3),
               
               rep("platelets", 3),
               
               rep("T helper", 3),
               rep("naive CD4 T", 3),
               rep("naive CD8 T", 3),
               rep("cytotoxic T", 3),
               rep("natural killer T", 3),
               rep("T regulatory", 3),
               rep("T central memory", 3),
               rep("natural killer", 3),
               # rep("all_T", 2),
               
               rep("B", 3))
# stacks : 
Treatment <- rep(c("alpha","lambda","untreated"), 15)

# create data frame :  
Prop <- data.frame(Treatment, CellTypesProp, proportion)

# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
Prop$CellTypesProp <- factor(Prop$CellTypesProp, levels = c(
  "monocytes", "neutrophils", "DCs","mDCs","pDCs","platelets",
  "T helper","naive CD4 T","naive CD8 T","cytotoxic T","natural killer T","T regulatory","T central memory","natural killer",
  "B"))
Prop$Treatment <- factor(Prop$Treatment, levels = c("alpha", "lambda", "untreated"))

colours <- c("lightgrey","#6ab5ba","darkgrey")


# Plot : 
# Plot : 
proportion <- ggplot(
  Prop,
  aes(fill = Treatment, y = CellTypesProp, x = proportion)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = colours) +
  labs(title = "", x = "Number of cells", y = "", fill = "IFN treatment") +
  theme(
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),legend.title = element_text(size = 14, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 14, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )
proportion






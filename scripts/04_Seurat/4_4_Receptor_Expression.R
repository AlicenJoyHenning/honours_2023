# Receptor expression across cell types 

# [1] Load integrated data set :
TreatmentAnnotated <- read_rds("honours/results/FinalIndex/TreatmentAnnotated.rds")
Idents(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated)

# [2] View receptor expression across cell types : 

receptorFP <- FeaturePlot(TreatmentAnnotated,
                          features = c("IFNAR1"),
                          cols = c("#e2e2e2", "black"))

receptorDP <- DotPlot(TreatmentAnnotated,
                          features = c("IFNAR1",  "IFNLR1"),
                          cols = c("black", "black")) +
  coord_flip()
  
                          
#blend to view alpha and lambda receptors!
?FeaturePlot()

# [3] Combining T cell populations 

ReceptorTreatment <- RenameIdents(TreatmentAnnotated, 
                                  'naive CD4 T' = "T cells",
                                  'T helper' = "T cells",
                                  'naive CD8 T' = "T cells",
                                  'cytotoxic T' = "T cells", 
                                  'NKT' = "T cells",
                                  'Tregs' = "T cells",
                                  'NK' = "T cells",
                                  'Tcm' = "T cells",
                                  'Tcm'= "T cells",
                                  "B" = 'B cells')


receptorDP <- DotPlot(ReceptorTreatment,
                      features = c("IFNAR1",  "IFNLR1"),
                      cols = c("black", "black")) +
  coord_flip()

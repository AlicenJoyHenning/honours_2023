# Venn Diagrams 
# load Venn diagram package
library("VennDiagram")

# move to new plotting page
grid.newpage()

# M!
draw.pairwise.venn(area1=904, area2=0,cross.area=0,
                   category=c("A","L"),
                   fill=c("Red","Yellow"))


# M2
draw.pairwise.venn(area1=519, area2=47,cross.area=2,
                   category=c("A","L"),
                   fill=c("Red","Yellow"))


# L1
L1 <- draw.pairwise.venn(area1=345, area2=3,cross.area=0,
                     category=c("A","L"),
                     fill=c("Red","Yellow"))
# L2 

L2 <- draw.pairwise.venn(area1=10, area2=10,cross.area=5,
                         category=c("A","L"),
                         fill=c("Red","Yellow"))
  
  
  
L6 <- draw.pairwise.venn(
  area1 = 10,               # Standard circle size
  area2 = 10,               # Standard circle size
  cross.area = 5,           # Standard intersection size
  category = c("316", "12"),  # Labels for circles
  fill = c("Red", "Yellow"),
  category.pos = c(0, 0),   # Position of category labels
  category.dist = 0.1,      # Distance of category labels from the center
  category.just = list(c(0.5, 0.5)),  # Justification of category labels
  cat.cex = 2              # Font size of category labels
)







  
  
  
  
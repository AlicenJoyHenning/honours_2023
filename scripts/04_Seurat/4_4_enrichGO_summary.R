
# Summarizing all enrich GO output in 2 bar graphs 
# normal bar graph quantifying related enrich GO terms for groups of cell types 

library(ggplot2)

##### Alpha #####
# Define alpha terms to be compared across cell types : 
terms <- c("mono", "neu", "T", "B")

# Define counts associated with these terms for each PBMC population : 
counts1 <- c(25, 25, 23, 25) # virus
counts2 <- c(25, 22, 22, 27) # IFN 
counts3 <- c(17, 17, 17, 17) # cytokine
counts4 <- c(7, 7, 6, 6) # innate 
counts5 <- c(10, 5, 4, 2) # inflam
  
df1 <- data.frame(terms, counts1)
df2 <- data.frame(terms, counts2)
df3 <- data.frame(terms, counts3)
df4 <- data.frame(terms, counts4)
df5 <- data.frame(terms, counts5)

alphaplot1 <- ggplot(df1,
                     aes(y = counts1, x = terms)) +
  geom_bar(fill = "lightgrey", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "1. Viral\n responses", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 10)),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

alphaplot2 <- ggplot(df2,
                     aes(y = counts2, x = terms)) +
  geom_bar(fill = "lightgrey", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "2. Inteferon\n signaling", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

alphaplot3 <- ggplot(df3,
                     aes(y = counts3, x = terms)) +
  geom_bar(fill = "lightgrey", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "3. Cytokine related\n function", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

alphaplot4 <- ggplot(df4,
                     aes(y = counts4, x = terms)) +
  geom_bar(fill = "lightgrey", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "4. Innate Immune\n Responses ", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

alphaplot5 <- ggplot(df5,
                     aes(y = counts5, x = terms)) +
  geom_bar(fill = "lightgrey", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "5. Inflammatory\n responses ", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    # axis.title.y = element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

alphaplots <- alphaplot1 | alphaplot2 | alphaplot3 | alphaplot4 | alphaplot5


##### Lambda #####
# Define alpha terms to be compared across cell types : 
terms <- c("mono", "neu", "T", "B")

# Define counts associated with these terms for each PBMC population : 
counts1 <- c(0, 0, 0, 25) # IFN
counts2 <- c(0, 0, 0, 22) # virus
counts3 <- c(0, 9, 0, 0) # vit d
counts4 <- c(6, 3, 0, 0) # bone 
counts5 <- c(0, 0, 4, 0) # RNA

df1 <- data.frame(terms, counts1)
df2 <- data.frame(terms, counts2)
df3 <- data.frame(terms, counts3)
df4 <- data.frame(terms, counts4)
df5 <- data.frame(terms, counts5)

lambdaplot1 <- ggplot(df1,
                     aes(y = counts1, x = terms)) +
  geom_bar(fill = "#6ab5ab", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "1. Inteferon\n signaling", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 10)),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

lambdaplot2 <- ggplot(df2,
                     aes(y = counts2, x = terms)) +
  geom_bar(fill = "#6ab5ab", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "2. Viral\n responses", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

lambdaplot3 <- ggplot(df3,
                     aes(y = counts3, x = terms)) +
  geom_bar(fill = "#6ab5ab", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "3. Vitamin D\n responses", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

lambdaplot4 <- ggplot(df4,
                     aes(y = counts4, x = terms)) +
  geom_bar(fill = "#6ab5ab", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "4. Bone\n mineralization", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

lambdaplot5 <- ggplot(df5,
                     aes(y = counts5, x = terms)) +
  geom_bar(fill = "#6ab5ab", position = "stack", stat = "identity") + 
  theme_minimal() +
  labs(title = "5. RNA\n processing ", x = "", y = "GO term counts") +
  theme(
    axis.text.x = element_text(size = 14, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    # axis.title.y = element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 14, face = "italic", hjust = 0.5),  # Center-align plot title and set its size
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

lambdaplots <- lambdaplot1 | lambdaplot2 | lambdaplot3 | lambdaplot4 | lambdaplot5


alphaplots / lambdaplots


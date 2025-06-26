# Installing necessary packages 
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("energy")
library(ggplot2)
library(dplyr)
library(energy)



labels <- read.csv("RNA_DATA/labels.csv", row.names = 1)
Y <- model.matrix(~ 0 + factor(labels[,1]))
X <- read.csv("RNA_DATA/data.csv", row.names = 1)

dc <- dcor(X, Y)
print(paste("Distance Correlation:", dc))



dcor_per_gene <- apply(X, 2, function(gene_expr) {dcor(gene_expr,Y)})
sorted_dcor <- sort(dcor_per_gene, decreasing = TRUE)
#print(sorted_dcor)
head(sorted_dcor, 10)

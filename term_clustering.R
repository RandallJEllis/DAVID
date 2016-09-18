library(GOSemSim)
library(pvclust)
library(RColorBrewer)

divergent <- brewer.pal(11, "RdBu")
highlight <- brewer.pal(3, "Set2")[1:2]

df <- read.csv("most_common_terms.csv", header=FALSE)
colnames(df) <- c("go", "term", "count")

df <- df[1:75, ]
davidmatrix <- termSim(df['go'], df['go'], method="Wang", ont = 'BP', organism = 'human')
# retain rows with less than 2 NAs
davidmatrix2 <- davidmatrix[rowSums(is.na(davidmatrix)) < 2, ]
# retain columns with less than 2 NAs
davidmatrix3 <- davidmatrix2[, colSums(is.na(davidmatrix2)) < 2]
# davidmatrix.nona <- na.omit(davidmatrix)
# davidmatrix.scaled <- scale(davidmatrix)

d <- dist(davidmatrix3, method = "euclidean")
d.nona <- na.omit(d)
fit <- hclust(d, method="ward")
plot(fit)
groups <- cutree(fit, k=6) # cut tree into 5 clusters
rect.hclust(fit, k=6, border="red")

identical(fit$labels, colnames(davidmatrix3))
# if the above line is TRUE, then equate fit$labels with df['term']. However, the removed row is still
# present in df, so you must remove it. 
df2 <- df[-59, ]
fit$labels <- df2['term']
plot(fit)

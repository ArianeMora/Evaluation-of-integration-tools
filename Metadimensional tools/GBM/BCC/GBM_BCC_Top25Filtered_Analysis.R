# BCC works only when further filtering of the data is done.
# 
## ------------------------------------------------------------------------
# Doesn't work for R = 4.0.3
library(bayesCC)
gbm.dat.me = readRDS("GBM_Methylation_Top25Filtered.rd")
gbm.dat.rna = readRDS("GBM_RNA_Top25Filtered.rd")
gbm.dat.protein = readRDS("GBM_Protein_Top25Filtered.rd")
data.labels = readRDS("GBM_DataLabels.rd")


# Data details
## ------------------------------------------------------------------------
mrna2 <- t(gbm.dat.rna)
me2 <- t(gbm.dat.me)
protein2 <- t(gbm.dat.protein)
X = list(mrna2, me2, protein2)
lapply(X,dim)
# [[1]]
# [1] 3011  100
# 
# [[2]]
# [1] 5744  100
# 
# [[3]]
# [1]  33 100

Results = list()
set.seed(7439285)
tic("BCC GBM")
Results[[1]] = BCC(X,K=4)  
toc()
# BCC GBM: 900.116 sec elapsed
K4 = Results[[1]]
Clusters = rep(1,100)
Clusters[K4$Cbest[,2]==1] = 2
Clusters[K4$Cbest[,3]==1] = 3
Clusters[K4$Cbest[,4]==1] = 4
clusters_label <- data.frame(data.labels$Label)
colnames(clusters_label) <- "Samples"
clusters_label$Samples <- as.character(clusters_label$Samples)
table(clusters_label$Samples)

# Classical Mesenchymal      Neural   Proneural 
# 28          24          17          31 

clusters_label$Truth1 <- clusters_label$Samples
clusters_label$Clusters <- Clusters
confusion.mat <- table(Truth = clusters_label$Truth1, Pred = clusters_label$Clusters)
confusion.mat
#              Pred
# Truth          1  2  3  4
# Classical      0  0  4 24
# Mesenchymal    1  0 20  3
# Neural         1  0  9  7
# Proneural      4  1 23  3

## CLUSTER 1= MESENCHYMAL; CLUSTER 2 = NEURAL; CLUSTER 3= PRONEURAL; CLUSTER 4 = CLASSICAL
clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth <- gsub("Classical", "4", clusters_label$Truth)
clusters_label$Truth <- gsub("Mesenchymal", "1", clusters_label$Truth)
clusters_label$Truth <- gsub("Neural", "2", clusters_label$Truth)
clusters_label$Truth <- gsub("Proneural", "3", clusters_label$Truth)
confusion.mat <- table(Truth = clusters_label$Truth, Pred = clusters_label$Clusters)
confusion.mat
#      Pred
# Truth  1  2  3  4
# 1      1  0 20  3
# 2      1  0  9  7
# 3      4  1 23  3
# 4      0  0  4 24
 
BCC_Accuracy = sum(clusters_label$Truth == clusters_label$Clusters)/length(clusters_label$Clusters)
round(BCC_Accuracy,4)
# [1] 0.48

# MESENCHYMAL Samples (Cluster 1)
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0

# NEURAL (Cluster 2) NO samples in cluster 2
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# PRONEURAL (Cluster 3)
Pg3 <- round(confusion.mat[4,4]/sum(confusion.mat[,4]),2) # Precision
Rg4 <- round(confusion.mat[4,4]/sum(confusion.mat[4,]),2) # Recall
Fg4 <- round(2*((Pg4*Rg4)/(Pg4+Rg4)),2) # F1-score
Fg4[is.nan(Fg4)] <- 0

# CLASSICAL samples (Cluster 4)
Pg4 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0

AvgF <- round((Fg1+Fg2+Fg3+Fg4)/4,2)
setwd("/home/anita/Benchmarking/GBM/BCC")
save.image(file = "GBM_BCC_Results.RData")


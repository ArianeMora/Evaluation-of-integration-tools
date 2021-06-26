## iClusterPlus multi-class classification analysis for GBM data 
# We how good the results are for 4 cluster solution. Hence k =3 as number of clusters=k+1

library(iClusterPlus)
library(lattice)
library(tictoc)

# Loading data
gbm.dat.me = readRDS("Data/GBM_Methylation_Top25Filtered.rd")
gbm.dat.rna = readRDS("Data/GBM_RNA_Top25Filtered.rd")
gbm.dat.protein = readRDS("Data/GBM_Protein_Top25Filtered.rd")
data.labels = readRDS("Data/GBM_DataLabels.rd")

# formatting data
mrna2 = as.matrix(gbm.dat.rna)
protein2 = as.matrix(gbm.dat.protein)
me2 = as.matrix(gbm.dat.me)

## Running iClusterPlus
set.seed(999)
k = 3 # number of clusters is k+1

tic("iClusterGBM")
cv.fit = tune.iClusterPlus(cpu = 36, dt1 = mrna2, dt2 = protein2,
                           dt3 = me2, type = c( "gaussian","gaussian", "gaussian"),
                           n.lambda = NULL, K=k)
toc()

output = alist()
output[[1]] <- cv.fit

nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}

best.fit=output[[1]]$fit[[which.min(BIC[,1])]]

# important features if required later
features = alist()
features[[1]] <- colnames(mrna2)
features[[2]] <- colnames(protein2)
features[[3]] <- colnames(me2)
sigfeatures = alist()
for (i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures) <- c("mRNA", "Protein", "Methylation")

clusters_label <- getClusters(output)
clusters_label <- data.frame(clusters_label)
colnames(clusters_label)[1] <- "Predicted"

dataname = "GBM"
save(best.fit, sigfeatures, file = paste(dataname,".best.cv.fit.k", k,".RData", sep = ""))
save.image(file = "iCP/GBM_iCPAnalysis_Top25Filtered.RData")

## Run the below after loading GBM_iCPAnalysis.RData
load("iCP/GBM_iCPAnalysis_Top25Filtered.RData")

clustersGBM = clusters_label
clustersGBM$Truth <- data.labels$Label
clustersGBM$Truth <- gsub("Classical", "1", clustersGBM$Truth)
clustersGBM$Truth <- gsub("Mesenchymal", "2", clustersGBM$Truth)
clustersGBM$Truth <- gsub("Neural", "3", clustersGBM$Truth)
clustersGBM$Truth <- gsub("Proneural", "4", clustersGBM$Truth)

confusion.mat =table(truth=clustersGBM$Truth, predicted = clustersGBM$Predicted)
confusion.mat
#         predicted
# truth  1  2  3  4
#   1   13  4 11  0
#   2    7 10  7  0
#   3    5  5  7  0
#   4    4  9  7 11

# Classical (Cluster 1)
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0

# Mesenchymal (Cluster 2)
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# Neural (Cluster 3)
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0

# Proneural (Cluster 4)
Pg4 <- round(confusion.mat[4,4]/sum(confusion.mat[,4]),2) # Precision
Rg4 <- round(confusion.mat[4,4]/sum(confusion.mat[4,]),2) # Recall
Fg4 <- round(2*((Pg4*Rg4)/(Pg4+Rg4)),2) # F1-score
Fg4[is.nan(Fg4)] <- 0

Accuracy <-  sum(diag(confusion.mat))/60
Avg.F1 <- (Fg1+Fg2+Fg3+Fg4)/4
save.image(file = "iCP/GBM_iCPAnalysis_withPerformanceScore.RData")


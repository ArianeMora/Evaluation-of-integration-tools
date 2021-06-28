# analysis of Melanoma dataset using iGC using Oncodrive Input
library(iGC)
# find_cna_driven_gene
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset_OncodriveCISInput.Rdata")

skcm.cnv[skcm.cnv==2] <- 1
skcm.cnv[skcm.cnv==-2] <- -1
skcm.cnv <- cbind(rownames(skcm.cnv), skcm.cnv)
colnames(skcm.cnv)[1] <- "GENE"
skcm.cnv[,1] <- as.character(skcm.cnv[,1])

skcm.ge <- cbind(rownames(skcm.ge), skcm.ge)
colnames(skcm.ge)[1] <- "GENE"
skcm.ge[,1] <- as.character(skcm.ge[,1])

# Change the data types to comply with the tool requirement
library(data.table)
setDT(skcm.cnv)
setDT(skcm.ge)
class(skcm.cnv)
class(skcm.ge)

# running iGC
library(tictoc)
tic("iGC")
cna_driven_genes <- find_cna_driven_gene(skcm.cnv, skcm.ge)
toc()

# Computing gain/loss sample proportion ...
# 
# Computing CNA gain driven gene records ... 
# 
# |================================| 100%
# Computing CNA loss driven gene records ... 
# 
# |================================| 100%
# > toc()
# iGC: 31.496 sec elapsed

# saving results
setwd("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/")
write.table(cna_driven_genes$gain_driven, file = "iGC_SKCM_GainDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_SKCM_LossDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_SKCM_BothDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )

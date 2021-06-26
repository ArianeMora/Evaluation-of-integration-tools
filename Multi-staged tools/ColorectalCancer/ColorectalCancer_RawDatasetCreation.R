## Creating colorectal cancer dataset
rm(list=ls())
library(caret)
### reading thresholded copy number downloaded from UCSC Xena browser 
# https://xenabrowser.net/datapages/?cohort=TCGA%20Colon%20Cancer%20(COAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# https://tcga-xena-hub.s3.us-east-1.amazonaws.com/latest/TCGA.COAD.sampleMap%2FSNP6_nocnv_genomicSegment.gz
# https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FSNP6_nocnv_genomicSegment&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

cn <- read.table(gzfile("COAD_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"), header = T, sep = "\t")
cn[1:5,1:5]
rownames(cn) <- cn[,1]
cn <- cn[,-1]

### reading gene expression file downloaded from UCSC Xena browser
ge.exp <- read.table(gzfile("COAD_HiSeqV2.gz"), header = T, sep = "\t")
ge.exp[1:5,1:5]
rownames(ge.exp) <- ge.exp[,1]
ge.exp <- ge.exp[,-1]

### reading methylation data downloaded from firebrowse
me.exp <- read.table("COAD.meth.by_mean.data.txt", header = T, sep = "\t")
me.exp[1:5,1:5]
me.exp <- me.exp[-1,]
rownames(me.exp) <- me.exp[,1]
me.exp <- me.exp[,-1]
me.exp <- na.omit(me.exp)

## identifying common genes and samples in all the datasets
a <- intersect(colnames(cn), colnames(ge.exp))
a1 <- intersect(colnames(me.exp), a)
b <- intersect(rownames(cn), rownames(ge.exp))
b1 <- intersect(rownames(me.exp), b)

ge.exp2 <- ge.exp[b1,a1]
cn2 <- cn[b1,a1]
me.exp2 <- me.exp[b1,a1]

## removing genes with zero and near zero variance in gene expression dataset and creating final dataset
# ge <-  Gene expression dataframe
# cnv <- copy number variation dataframe
# me <- methylation intensity dataframe
tmp <- t(ge.exp2)
rmcols <- nearZeroVar(tmp)
coad.ge <- t(tmp[, -rmcols])
coad.ge <- data.frame(coad.ge)

coad.cnv <- cn2[rownames(coad.ge),]
coad.me <- me.exp2[rownames(coad.ge),]

coad.cnv[1:5,1:5]
coad.ge[1:5,1:5]
coad.me[1:5,1:5]

# changing me columns as numeric
coad.me <- apply(coad.me, 2, as.numeric)
rownames(coad.me) <- row.names(coad.ge)

# final check
all(rownames(coad.ge) == rownames(coad.cnv))
all(rownames(coad.cnv) == rownames(coad.me))

all(colnames(coad.ge) == colnames(coad.cnv))
all(colnames(coad.cnv) == colnames(coad.me))

# saving the data
save(file = "ColorectalCancerRawDataset.Rdata", coad.ge, coad.cnv, coad.me)

## Creating pancreatic cancer dataset
rm(list=ls())
library(caret)

### reading thresholded copy number values downloaded from UCSC Xena browser
# https://xenabrowser.net/datapages/?dataset=TCGA.COAD.sampleMap%2FSNP6_nocnv_genomicSegment&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
cn <- read.table(gzfile("PAAD_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"), header = T, sep = "\t")
cn[1:5,1:5]
rownames(cn) <- cn[,1]
cn <- cn[,-1]

### reading gene expression file downloaded from UCSC Xena
ge.exp <- read.table(gzfile("PAAD_HiSeqV2.gz"), header = T, sep = "\t")
ge.exp[1:5,1:5]
rownames(ge.exp) <- ge.exp[,1]
ge.exp <- ge.exp[,-1]

### reading methylation data downloaded from firebrowse
me.exp <- read.table("PAAD.meth.by_mean.data.txt", header = T, sep = "\t")
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

cn2[1:5,1:5]
ge.exp2[1:5,1:5]
me.exp2[1:5,1:5]

## removing genes with zero and near zero variance in gene expression dataset and creating final dataset
# ge <-  Gene expression dataframe
# cnv <- copy number variation dataframe
# me <- methylation intensity dataframe
tmp <- t(ge.exp2)
rmcols <- nearZeroVar(tmp)
paad.ge <- t(tmp[, -rmcols])
paad.ge <- data.frame(paad.ge)

paad.cnv <- cn2[rownames(paad.ge),]
paad.me <- me.exp2[rownames(paad.ge),]

paad.cnv[1:5,1:5]
paad.ge[1:5,1:5]
paad.me[1:5,1:5]

# changing me columns as numeric
paad.me <- apply(paad.me, 2, as.numeric)
rownames(paad.me) <- row.names(paad.ge)

# final check
all(rownames(paad.ge) == rownames(paad.cnv))
all(rownames(paad.cnv) == rownames(paad.me))

all(colnames(paad.ge) == colnames(paad.cnv))
all(colnames(paad.cnv) == colnames(paad.me))

# save the dataset
save(file = "PancreaticCancerRawDataset.Rdata", paad.ge, paad.cnv, paad.me)

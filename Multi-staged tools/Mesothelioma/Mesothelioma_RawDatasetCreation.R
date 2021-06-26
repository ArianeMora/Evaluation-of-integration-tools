# -------------------------------------------------------------------------------------------
#            Dataset download for Mesothelioma (25/06/2021) Ariane Mora
#     
#      Main site: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Mesothelioma%20(MESO)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#
#     CNV: TCGA-MESO.gistic.tsv: copy number (gene-level) - GISTIC - focal score by gene: https://xenabrowser.net/datapages/?dataset=TCGA-MESO.gistic.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#     DNA Metylation: TCGA-MESO.methylation450.tsv: DNA methylation - Illumina Human Methylation 450: https://xenabrowser.net/datapages/?dataset=TCGA-MESO.methylation450.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#     RNAseq counts: TCGA-MESO.htseq_counts.tsv: gene expression RNAseq - HTSeq - Counts: https://xenabrowser.net/datapages/?dataset=TCGA-MESO.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# 
# -------------------------------------------------------------------------------------------

## Creating mesothelioma dataset
rm(list=ls())
library(caret)

### reading copy number values and seperating them into amplification and deletion calls
cn <- read.table(gzfile("TCGA-MESO.gistic.tsv.gz"), header = T, sep = "\t")
cn[1:5,1:5]
rownames(cn) <- cn[,1]
cn <- cn[,-1]

### reading gene expression file
ge.exp <- read.table(gzfile("TCGA-MESO.htseq_counts.tsv.gz"), header = T, sep = "\t")
ge.exp[1:5,1:5]
rownames(ge.exp) <- ge.exp[,1]
# Set Ensembl ID to be the main ID
ge.exp <- ge.exp[,2:length(colnames(ge.exp))]
# Remove the normalisation by log2 + 1
ge.exp <- (2^ge.exp) - 1 

### reading methylation data
me.exp <- read.table("TCGA-MESO.methylation450.tsv.gz", header = T, sep = "\t")
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
meso.ge <- t(tmp[, -rmcols])
meso.ge <- data.frame(meso.ge)

meso.cnv <- cn2[rownames(meso.ge),]
meso.me <- me.exp2[rownames(meso.ge),]

meso.cnv[1:5,1:5]
meso.ge[1:5,1:5]
meso.me[1:5,1:5]

# changing me columns as numeric
meso.me <- apply(meso.me, 2, as.numeric)
rownames(meso.me) <- row.names(meso.ge)

# final check
all(rownames(meso.ge) == rownames(meso.cnv))
all(rownames(meso.cnv) == rownames(meso.me))

all(colnames(meso.ge) == colnames(meso.cnv))
all(colnames(meso.cnv) == colnames(meso.me))

save(file = "MesotheliomaRawDataset.Rdata", meso.ge, meso.cnv, meso.me)

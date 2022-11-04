install.packages("BiocManager")
library("BiocManager")
BiocManager::install("BiocManager")
BiocManager::install("maSigPro")
BiocManager::install("edgeR")
BiocManager::install("GenomicFeatures")
BiocManager::install("rhdf5")
install.packages("tximport")
install.packages("readr")
install.packages("tidyr")
install.packages("dplyr")

library("tximport")
library("readr")
library("tidyr")
library("dplyr")
library("maSigPro")
library("edgeR")
library("GenomicFeatures")
library("rhdf5")

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/masigpro"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

#Create a text file with the names of the Kallisto folders and assign as "samples"
sample <-c('JMIK','JMIL','JMIM','JMIN','JMIP','JMGU','JMGW','JMGX','JMGY','JMGZ','JMHA','JMHB','JMIC','JMID','JMIE','JMIF','JMHH','JMHI','JMHJ','JMHK','JMHL','JMHM','JMHN','JMHP','JMHC','JMHD','JMHE','JMHG')

#Create file path for abundance.tsv file in each of the folders that are specified in "samples"
files <- file.path(datapath, sample,"abundance.tsv")

#pull info from gff
txdb<- makeTxDbFromGFF("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Dionaea_muscipula.gff")
seqinfo(txdb)
genes(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
tail(tx2gene)
tx2genedm <- as.vector(tx2gene)
head (tx2genedm)

#extracts counts data and builds matrix with geneID and counts
txi <- tximport(files, type="kallisto", tx2gene = tx2genedm, txOut = TRUE, geneIdCol =
                  TRUE, abundanceCol = FALSE, countsCol = TRUE)
head(txi)

#create variable with just the counts
counts <- txi$counts

#rename "sample1", "sample2", and so on to thea actual names of the samples
colnames(counts) <- c('LT_Time5_1',
'LT_Time5_2',
'LT_Time5_3',
'LT_Time5_4',
'LT_Time60_1',
'LT_Time60_2',
'LT_Time60_3',
'LT_Time60_4',
'LT_Time180_1',
'LT_Time180_2',
'LT_Time180_3',
'LT_Time180_4',
'LT_Time720_1',
'LT_Time720_2',
'LT_Time720_3',
'LT_Time720_4',
'LT_Time1440_1',
'LT_Time1440_2',
'LT_Time1440_3',
'LT_Time1440_4',
'LT_Time2880_1',
'LT_Time2880_2',
'LT_Time2880_3',
'LT_Time2880_4',
'LT_Time4320_1',
'LT_Time4320_2',
'LT_Time4320_3',
'LT_Time4320_4')


#now that we have a counts matrix we will normalize the counts to TMM using EdgeR
#first assign groups
groups<-c('LT_Time5_1',
'LT_Time5_2',
'LT_Time5_3',
'LT_Time5_4',
'LT_Time60_1',
'LT_Time60_2',
'LT_Time60_3',
'LT_Time60_4',
'LT_Time180_1',
'LT_Time180_2',
'LT_Time180_3',
'LT_Time180_4',
'LT_Time720_1',
'LT_Time720_2',
'LT_Time720_3',
'LT_Time720_4',
'LT_Time1440_1',
'LT_Time1440_2',
'LT_Time1440_3',
'LT_Time1440_4',
'LT_Time2880_1',
'LT_Time2880_2',
'LT_Time2880_3',
'LT_Time2880_4',
'LT_Time4320_1',
'LT_Time4320_2',
'LT_Time4320_3',
'LT_Time4320_4')

#"Creates a DGEList object from a table of counts, group indicator for each column"
ddm<-DGEList(counts, group=groups)

#"Calculate scaling factors to convert raw library sizes into effective library sizes"
ddm<-calcNormFactors(ddm, method="TMM")

#"Compute counts per million (CPM)"
normdm<-cpm(ddm, normalize.lib.sizes=TRUE)

#"Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples"
#to identify outliers
pdf(file="MaSigProPreySeriesMDS.pdf")
plotMDS(normdm, labels=groups)
dev.off()

#outlier identified -> Trap_prey_time1440_rep10_4, remove using the subset function
normdm <- subset(normdm, select = - c(LT_Time1440_4))

#Removes genes that have all zero read counts
normdm <- normdm[rowSums(normdm)!=0, ]

#experimental design
Time <- rep(c(5,60,180,720,1440,2880,4320), each = 4)
Replicate <- rep(c(1:7), each = 4)
Prey <- rep(c(1), each = 28)
edesign.VFT <- cbind(Time,Replicate,Prey)
edesign.VFT <- as.data.frame(edesign.VFT)
rownames(edesign.VFT) <- groups
edesign.VFT <- edesign.VFT[-c(20), ] #remove outlier

#make exp. design for masigpro
design <- make.design.matrix(edesign.VFT, degree = 6)
design$groups.vector

#using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normdm, design, counts=TRUE, min.obs=0) #
NBt<-T.fit(NBp)

#remove influential genes
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normdm<-normdm[!rownames(normdm) %in% inf.genenames, ]

#redo regression with influential genes removed
NBp<-p.vector(normdm, design, counts=TRUE)
NBt<-T.fit(NBp)

#see how many clusters
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
pdf(file="MaSigProWSSClusterRange.pdf")
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")
dev.off()

#get list of significant genes
sigs <- get.siggenes(NBt,rsq=0.6, vars="all")

#visualize clusters
pdf(file="MaSigProPreySeriesCluster6.pdf")
see.genes(sigs$sig.genes, k = 6)
dev.off()

k6vis<- see.genes(sigs$sig.genes, k = 6)

#get genes in each cluster
cluster1 <- names(which(k6vis$cut==1))
cluster2 <- names(which(k6vis$cut==2))
cluster3 <- names(which(k6vis$cut==3))
cluster4 <- names(which(k6vis$cut==4))
cluster5 <- names(which(k6vis$cut==5))
cluster6 <- names(which(k6vis$cut==6))


#make cluster csv's
write.csv(x = cluster1, file = "MaSigProPreyCluster1.csv")
write.csv(x = cluster2, file = "MaSigProPreyCluster2.csv")
write.csv(x = cluster3, file = "MaSigProPreyCluster3.csv")
write.csv(x = cluster4, file = "MaSigProPreyCluster4.csv")
write.csv(x = cluster5, file = "MaSigProPreyCluster5.csv")
write.csv(x = cluster6, file = "MaSigProPreyCluster6.csv")

chooseCRANmirror(ind = 80)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("BiocManager", quietly = TRUE))
    install.packages("edgeR")

install.packages("tximport")
install.packages("readr")
install.packages("tidyr")
install.packages("dplyr")

library("tximport")
library("readr")
library("tidyr")
library("dplyr")

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant/"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/masigpro"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

#Create a text file with the names of the Kallisto folders and assign as "samples"
sample <- c('JMIJ'
'JMHF'
'JMHQ'
'JMIB'
'JMIK'
'JMIL'
'JMIM'
'JMIN'
'JMIP'
'JMGU'
'JMGW'
'JMGX'
'JMGY'
'JMGZ'
'JMHA'
'JMHB'
'JMIC'
'JMID'
'JMIE'
'JMIF'
'JMHH'
'JMHI'
'JMHJ'
'JMHK'
'JMHL'
'JMHM'
'JMHN'
'JMHP'
'JMHC'
'JMHD'
'JMHE'
'JMHG')

#Create file path for abundance.tsv file in each of the folders that are specified in "samples"
files <- file.path(datapath, sample)

#extracts counts data and builds matrix with geneID and counts
txi <- tximport(files, type="kallisto", geneIdCol =TRUE, countsCol = TRUE)
head(txi)

#create variable with just the counts
counts <- txi$counts

#rename "sample1", "sample2", and so on to thea actual names of the samples
colnames(counts) <- c('LeafTrap_Time0_1',
'LT_Time0_2',
'LT_Time0_3',
'LT_Time0_4',
'LT_Time5_1',
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
groups<-c('LeafTrap_Time0_1',
'LT_Time0_2',
'LT_Time0_3',
'LT_Time0_4',
'LT_Time5_1',
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

#Removes genes that have all zero read counts
normdm<-normdm[ rowSums(normdm)!=0, ]

#importing experimental design
Time <- rep(c(0,5,60,180,720,1440,2880,4320), each = 4)
Replicates <- rep(c(1:8), each = 4)
edesign.VFT <- cbind(Time,Replicates,groups)
edesign.VFT <- as.data.frame(edesign.VFT)

#matches column names of counts matrix to edesign matrix
colnames(normdm) <- rownames(edesign.VFT)

#make exp. design for masigpro
design <- make.design.matrix(edesign.VFT, degree = 7)
design$groups.vector

#using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normdm, design, counts=TRUE) #
NBt<-T.fit(NBp)

#get list of significant genes
sigs <- get.siggenes(NBt)

#visualize clusters
k8vis<-see.genes(sigs$sig.genes, k = 8)
k6vis<-see.genes(sigs$sig.genes, k = 6)

pdf(file="MaSigProPreySeriesCluster8.pdf")
k8vis
dev.off()

pdf(file="MaSigProPreySeriesCluster8.pdf")
k6vis
dev.off()

#get genes in each cluster
cluster1 <- names(which(k8vis$cut==1))
cluster2 <- names(which(k8vis$cut==2))
cluster3 <- names(which(k8vis$cut==3))
cluster4 <- names(which(k8vis$cut==4))
cluster5 <- names(which(k8vis$cut==5))
cluster6 <- names(which(k8vis$cut==6))
cluster7 <- names(which(k8vis$cut==7))
cluster8 <- names(which(k8vis$cut==8))

#make cluster csv's
write.csv(x = cluster1, file = "MaSigProCluster1.csv")
write.csv(x = cluster2, file = "MaSigProCluster2.csv")
write.csv(x = cluster3, file = "MaSigProCluster3.csv")
write.csv(x = cluster4, file = "MaSigProCluster4.csv")
write.csv(x = cluster5, file = "MaSigProCluster5.csv")
write.csv(x = cluster6, file = "MaSigProCluster6.csv")
write.csv(x = cluster7, file = "MaSigProCluster7.csv")
write.csv(x = cluster8, file = "MaSigProCluster8.csv")

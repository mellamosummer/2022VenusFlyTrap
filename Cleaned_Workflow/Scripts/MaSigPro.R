suppressMessages({library("edgeR")
library("maSigPro")
library("tximport")
library("readr")
library("tidyr")
library("dplyr")})

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
str(normdm)
colnames(normdm)
rownames(edesign.VFT)
colnames(edesign.VFT)
rownames(normdm)

design <- make.design.matrix(edesign.VFT, degree = 7)
design$groups.vector

##using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normdm, design, counts=TRUE) #please choose dis. family with care. default for counts is neg. binomial
NBt<-T.fit(NBp)

sigs <- get.siggenes(NBt, rsq = 0.6, vars = "groups")

pdf(file="MaSigProPreySeriesCluster8.pdf")
see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
cluster.method="hclust" ,cluster.data = 1, k = 8)
dev.off()

pdf(file="MaSigProPreySeriesCluster6.pdf")
see.genes(sigs$sig.genes, show.fit = T, dis =design$dis,
cluster.method="hclust" ,cluster.data = 1, k = 6)
dev.off()

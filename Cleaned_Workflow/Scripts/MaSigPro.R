##################################
#  INSTALL PACKAGES
##################################

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

##################################
#  SET DIRECTORIES
##################################

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/masigpro"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

##################################
#  PREY TIME SERIES
##################################

##################################
#  BUILD COUNTS MATRIX
##################################

#Create a text file with the names of the Kallisto folders and assign as "sample"
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

##################################
#  NORMALIZE COUNTS MATRIX
##################################

#"Creates a DGEList object from a table of counts, group indicator for each column"
ddm<-DGEList(counts, group=groups)

#"Calculate scaling factors to convert raw library sizes into effective library sizes"
ddm<-calcNormFactors(ddm, method="TMM")

#"Compute counts per million (CPM)"
normdm<-cpm(ddm, normalize.lib.sizes=TRUE)

##################################
#  PLOT MDS OF SAMPLES TO QC
##################################

#"Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples"
#to identify outliers
pdf(file="MaSigProPreySeriesMDS.pdf")
plotMDS(normdm, labels=groups)
dev.off()

#outlier identified -> Trap_prey_time1440_rep10_4, remove using the subset function
normdm <- subset(normdm, select = - c(LT_Time1440_4))

#Removes genes that have all zero read counts
normdm <- normdm[rowSums(normdm)!=0, ]

##################################
#  MAKE MASIGPRO E. DESIGN
##################################

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

##################################
#  RUN TIME SERIES ANALYSIS
##################################

#using the negative binomial options in masigpro, calculate polynomial regressions for each gene
# p.vector generates list of FDR corrected significant genes
NBp<-p.vector(normdm, design, counts=TRUE, min.obs=0) #Global regression step: 18926 genes total -> 8577 genes
#T.fit selects the best regression model for each gene using stepwise regression.
NBt<-T.fit(NBp) #fits 8577 genes, 1463 with influential data

#remove influential genes
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normdm<-normdm[!rownames(normdm) %in% inf.genenames, ]

#redo regression with influential genes removed
NBp<-p.vector(normdm, design, counts=TRUE) #Global regression step: 17463 genes total -> 6924 genes
NBt<-T.fit(NBp)

#plot within sum of squares to determine cluster #
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
pdf(file="MaSigProWSSPreyClusterRange.pdf")
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")
dev.off()

#get list of significant genes
sigs <- get.siggenes(NBt,rsq=0.6, vars="all")
write.csv(x = sigs, file = "MaSigProPreySignificantGenes.csv")

##################################
#  CLUSTER VISUALIZATION
##################################

pdf(file="MaSigProPreySeriesCluster6.pdf")
see.genes(sigs$sig.genes, k = 6)
dev.off()

MaSigProPreySigProfiles<-sigs$sig.genes$sig.profiles
write.csv(x = MaSigProPreySigProfiles, file = "MaSigProPreySigProfiles.csv")
MaSigProPreySigPValues<-sigs$sig.genes$sig.pvalues
write.csv(x = MaSigProPreySigPValues, file = "MaSigProPreySigPValues.csv")
MaSigProPreySigCoefficients<-sigs$sig.genes$coefficients
write.csv(x = MaSigProPreySigCoefficients, file = "MaSigProPreySigCoefficients.csv")

##################################
#  CLUSTER OUTPUT
##################################

#get genes in each cluster
k6vis<- see.genes(sigs$sig.genes, k = 6)
cluster1 <- names(which(k6vis$cut==1))
cluster2 <- names(which(k6vis$cut==2))
cluster3 <- names(which(k6vis$cut==3))
cluster4 <- names(which(k6vis$cut==4))
cluster5 <- names(which(k6vis$cut==5))
cluster6 <- names(which(k6vis$cut==6))

#make cluster csv's
write.csv(x = cluster1, file = "MaSigProPreyCluster1.csv") #794 genes
write.csv(x = cluster2, file = "MaSigProPreyCluster2.csv") #393 genes
write.csv(x = cluster3, file = "MaSigProPreyCluster3.csv") #201 genes
write.csv(x = cluster4, file = "MaSigProPreyCluster4.csv") #123 genes
write.csv(x = cluster5, file = "MaSigProPreyCluster5.csv") #140 genes
write.csv(x = cluster6, file = "MaSigProPreyCluster6.csv") #149 genes

##################################
#  NO PREY TIME SERIES
##################################

##################################
#  BUILD COUNTS MATRIX
##################################

#Create a text file with the names of the Kallisto folders and assign as "sample"
sample <-c('JMIJ','JMHF','JMHQ','JMIB','JMGR','JMGT','JMHR','JMHS','JMGQ','JMHY','JMHZ','JMIA','JMHT','JMHU','JMHW','JMHX')

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
colnames(counts) <- c('LT_Time0_1',
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
'LT_Time1440_1',
'LT_Time1440_2',
'LT_Time1440_3',
'LT_Time1440_4')


#now that we have a counts matrix we will normalize the counts to TMM using EdgeR
#first assign groups
groups<-c('LT_Time0_1',
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
'LT_Time1440_1',
'LT_Time1440_2',
'LT_Time1440_3',
'LT_Time1440_4')

##################################
#  NORMALIZE COUNTS MATRIX
##################################

#"Creates a DGEList object from a table of counts, group indicator for each column"
ddm<-DGEList(counts, group=groups)

#"Calculate scaling factors to convert raw library sizes into effective library sizes"
ddm<-calcNormFactors(ddm, method="TMM")

#"Compute counts per million (CPM)"
normdm<-cpm(ddm, normalize.lib.sizes=TRUE)

##################################
#  PLOT MDS OF SAMPLES TO QC
##################################

#"Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples"
#to identify outliers
pdf(file="MaSigProNoPreySeriesMDS.pdf")
plotMDS(normdm, labels=groups)
dev.off()

#Removes genes that have all zero read counts
normdm <- normdm[rowSums(normdm)!=0, ]

##################################
#  MAKE MASIGPRO E. DESIGN
##################################

#experimental design
Time <- rep(c(0,5,60,1440), each = 4)
Replicate <- rep(c(1:4), each = 4)
NoPrey <- rep(c(1), each = 16)
edesign.VFT <- cbind(Time,Replicate,NoPrey)
edesign.VFT <- as.data.frame(edesign.VFT)
rownames(edesign.VFT) <- groups

#make exp. design for masigpro
design <- make.design.matrix(edesign.VFT, degree = 3)
design$groups.vector

##################################
#  RUN TIME SERIES ANALYSIS
##################################

#using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normdm, design, counts=TRUE, min.obs=0) #Global regression step: 18592 genes total -> 2821 genes
NBt<-T.fit(NBp) #240 genes with influential data

#remove influential genes
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normdm<-normdm[!rownames(normdm) %in% inf.genenames, ]

#redo regression with influential genes removed
NBp<-p.vector(normdm, design, counts=TRUE) #Global regression step: 18352 genes <- 2475 genes
NBt<-T.fit(NBp)

#plot within sum of squares to determine cluster #
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
pdf(file="MaSigProWSSNoPreyClusterRange.pdf")
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
plot(1:15, wss, type="b")
dev.off()

#get list of significant genes
sigs <- get.siggenes(NBt,rsq=0.6, vars="all")

##################################
#  CLUSTER VISUALIZATION
##################################

pdf(file="MaSigNoProPreySeriesCluster2.pdf")
see.genes(sigs$sig.genes, k = 2)
dev.off()

MaSigProNoPreySigProfiles<-sigs$sig.genes$sig.profiles
write.csv(x = MaSigProNoPreySigProfiles, file = "MaSigProNoPreySigProfiles.csv")
MaSigProNoPreySigPValues<-sigs$sig.genes$sig.pvalues
write.csv(x = MaSigProNoPreySigPValues, file = "MaSigProNoPreySigPValues.csv")
MaSigProNoPreySigCoefficients<-sigs$sig.genes$coefficients
write.csv(x = MaSigProNoPreySigCoefficients, file = "MaSigProNoPreySigCoefficients.csv")


##################################
#  CLUSTER OUTPUT
##################################

k2vis<- see.genes(sigs$sig.genes, k = 2)

#get genes in each cluster
cluster1 <- names(which(k2vis$cut==1)) #279 genes
cluster2 <- names(which(k2vis$cut==2)) #434 genes

#make cluster csv's
write.csv(x = cluster1, file = "MaSigProNoPreyCluster1.csv")
write.csv(x = cluster2, file = "MaSigProNoPreyCluster2.csv")

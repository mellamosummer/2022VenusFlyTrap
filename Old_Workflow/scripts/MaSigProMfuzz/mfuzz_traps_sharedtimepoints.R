library(edgeR)
library(maSigPro)
library(Mfuzz)
library(tximport)
library(readr)
library(tidyr)
library(GenomicFeatures)
library(dplyr)

#Create a text file with the names of the Kallisto folders and assign as "samples"
setwd("/Users/srblanco/Desktop/UGA/Dissertation/Leebens Mack Rotation/VFTAnalysis2/")
samples <- read.table("traps_prey_noprey_sametimepoint.txt", header = FALSE) #samples.txt is a txt file with a list of folders to be imported
samples
samples <- samples[-c(25),] 

#Create file path for abundance.tsv file in each of the folders that are specified in "samples"
files <- file.path("/Users/srblanco/Desktop/UGA/Dissertation/Leebens Mack Rotation/VFTAnalysis2/Kallisto_quant_results", samples, "abundance.tsv")
files
names(files) <- paste0("sample", 1:16) #name each sample
head(files)
file.exists(files) #all need to be TRUE

#Pulls annotation info from the gff for Dionaea muscipula 
txdb<- makeTxDbFromGFF("/Users/srblanco/Desktop/UGA/Dissertation/Leebens Mack Rotation/VFTAnalysis2/Dionaea_muscipula.gff")
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
colnames(counts) <- c("Trap_no_prey_time0005_rep3_1", "Trap_no_prey_time0005_rep3_2", "Trap_no_prey_time0005_rep3_3", "Trap_no_prey_time0005_rep3_4",
                      "Trap_no_prey_time0060_rep4_1", "Trap_no_prey_time0060_rep4_2", "Trap_no_prey_time0060_rep4_3", "Trap_no_prey_time0060_rep4_4",
                      "Trap_no_prey_time1440_rep5_1", "Trap_no_prey_time1440_rep5_2", "Trap_no_prey_time1440_rep5_3", "Trap_no_prey_time1440_rep5_4",
                      "Trap_prey_time0005_rep6_1",    "Trap_prey_time0005_rep6_2",    "Trap_prey_time0005_rep6_3",    "Trap_prey_time0005_rep6_4",  
                      "Trap_prey_time0060_rep7_1",    "Trap_prey_time0060_rep7_2",    "Trap_prey_time0060_rep7_3",    "Trap_prey_time0060_rep7_4",  
                      "Trap_prey_time1440_rep10_1",   "Trap_prey_time1440_rep10_2",   "Trap_prey_time1440_rep10_3",   "Trap_prey_time1440_rep10_4")
                      
#now that we have a counts matrix we will normalize the counts to TMM using EdgeR

setwd("/Users/srblanco/Desktop/UGA/Dissertation/Leebens Mack Rotation/VFTAnalysis2/")

#first assign groups 

groups<-c("Trap_no_prey_time0005_rep3_1", "Trap_no_prey_time0005_rep3_2", "Trap_no_prey_time0005_rep3_3", "Trap_no_prey_time0005_rep3_4",
          "Trap_no_prey_time0060_rep4_1", "Trap_no_prey_time0060_rep4_2", "Trap_no_prey_time0060_rep4_3", "Trap_no_prey_time0060_rep4_4",
          "Trap_no_prey_time1440_rep5_1", "Trap_no_prey_time1440_rep5_2", "Trap_no_prey_time1440_rep5_3", "Trap_no_prey_time1440_rep5_4",
          "Trap_prey_time0005_rep6_1",    "Trap_prey_time0005_rep6_2",    "Trap_prey_time0005_rep6_3",    "Trap_prey_time0005_rep6_4",  
          "Trap_prey_time0060_rep7_1",    "Trap_prey_time0060_rep7_2",    "Trap_prey_time0060_rep7_3",    "Trap_prey_time0060_rep7_4",  
          "Trap_prey_time1440_rep10_1",   "Trap_prey_time1440_rep10_2",   "Trap_prey_time1440_rep10_3",   "Trap_prey_time1440_rep10_4")

#"Creates a DGEList object from a table of counts, group indicator for each column"
ddm<-DGEList(counts, group=groups)

#"Calculate scaling factors to convert raw library sizes into effective library sizes"
ddm<-calcNormFactors(ddm, method="TMM")

#"Compute counts per million (CPM)"
normdm<-cpm(ddm, normalize.lib.sizes=TRUE)

#"Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples"
#to identify outliers
pdf("MDS_trapssharedtimepoints.pdf")
plotMDS(normdm, labels=groups)
dev.off()

#outlier identified -> Trap_prey_time1440_rep10_4, remove using the subset function
normdm <- subset(normdm, select = - c(Trap_prey_time1440_rep10_4))

#Removes genes that have all zero read counts
normdm<-normdm[ rowSums(normdm)!=0, ]


#importing experimental design
edesign.VFT <- read.csv("traps_prey_noprey_sametimepoint.csv")
edesign.VFT <- as.data.frame(edesign.VFT)
edesign.VFT
edesign.VFT <- edesign.VFT[-c(24), ] #remove outlier
edesign.VFT<- na.omit(edesign.VFT) 

names(edesign.VFT)[1] <- "row_names"
rownames(edesign.VFT) <- edesign.VFT$row_names
edesign.VFT<-edesign.VFT[-c(0,1)]
edesign.VFT

str(edesign.VFT)

#matches column names of counts matrix to edesign matrix
colnames(normdm) <- rownames(edesign.VFT) 
str(normdm)
colnames(normdm)
rownames(edesign.VFT)
colnames(edesign.VFT) 
rownames(normdm)

design <- make.design.matrix(edesign.VFT, degree = 2)
design$groups.vector


##using the negative binomial options in masigpro, calculate polynomial regressions for each gene
NBp<-p.vector(normdm, design, counts=TRUE) #please choose dis. family with care. default for counts is neg. binomial
NBt<-T.fit(NBp, epsilon=0.00001)

#Removes influential genes
influential<-NBt$influ.info
inf.genenames<-colnames(influential)
normdm<-normdm[!rownames(normdm) %in% inf.genenames, ]

##pick k
NBp<-p.vector(normdm, design, counts=TRUE)
wss<-(nrow(NBp$SELEC)-1)*sum(apply(NBp$SELEC,2,var))
for (i in 2:15) wss[i]<- sum(kmeans(NBp$SELEC, centers=i, iter.max=20)$withinss)
pdf("WSSclusterssametimepoints.pdf")
plot(1:15, wss, type="b")


NBt<-T.fit(NBp, epsilon=0.00001)
sigs <- get.siggenes(NBt, rsq = 0.5, vars = "all")

## Average replicates (required by MFuzz) and make into ExpressionSet
fit_averaged <- data.frame(row.names = row.names(NBp$SELEC))
fit_averaged$T1 <- rowMeans(NBp$SELEC[ , c(1:4,13:16)], na.rm=TRUE)
fit_averaged$T2 <- rowMeans(NBp$SELEC[ , c(5:8, 17:20)], na.rm=TRUE)
fit_averaged$T3 <- rowMeans(NBp$SELEC[ , c(9:12, 21:23)], na.rm=TRUE)
fit_averaged <- as.matrix(fit_averaged)
fit_averaged

##estimate m for mfuzz
genes <- ExpressionSet(assayData=fit_averaged)
genes.f <- filter.std(genes, min.std=0)
genes.s <- standardise(genes.f)

## Estimate m for mfuzz
m <- mestimate(genes.s)
m

## Find cluster number with minimum centroid distance
pdf("mfuzz_sharedtimepointstraps_Dmin_clusters_06032022.pdf")
dmin <- Dmin(genes.s, m=m, crange=seq(2,22,1), visu=TRUE)
dev.off()


## Do soft clustering

#2 clusters
clusters <- mfuzz(genes.s, c=2, m=m)

genesinclusters2 <- clusters$cluster
genesinclusters2 <- as.data.frame(genesinclusters2)

write.table(genesinclusters2, "geneclusters2.csv", sep=",")

numberofgenespercluster <- genesinclusters2 %>%
  group_by(genesinclusters2) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster2.csv")

pdf("mfuzz_sharedtimepointstraps_2clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()


#8 clusters
clusters <- mfuzz(genes.s, c=8, m=m)

genesinclusters8 <- clusters$cluster
genesinclusters8 <- as.data.frame(genesinclusters8)

write.table(genesinclusters8, "geneclusters8.csv", sep=",")

numberofgenespercluster <- genesinclusters8 %>%
  group_by(genesinclusters8) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster8.csv")

pdf("mfuzz_sharedtimepointstraps_8clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()

#9 clusters
clusters <- mfuzz(genes.s, c=9, m=m)

genesinclusters9 <- clusters$cluster
genesinclusters9 <- as.data.frame(genesinclusters9)

write.table(genesinclusters9, "geneclusters9.csv", sep=",")

numberofgenespercluster <- genesinclusters9 %>%
  group_by(genesinclusters9) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster9.csv")

pdf("mfuzz_sharedtimepointstraps_9clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()

#10 clusters
clusters <- mfuzz(genes.s, c=10, m=m)

genesinclusters10 <- clusters$cluster
genesinclusters10 <- as.data.frame(genesinclusters10)

write.table(genesinclusters10, "geneclusters10.csv", sep=",")

numberofgenespercluster <- genesinclusters10 %>%
  group_by(genesinclusters10) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster10.csv")

pdf("mfuzz_sharedtimepointstraps_10clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()

#11 clusters
clusters <- mfuzz(genes.s, c=11, m=m)

genesinclusters11 <- clusters$cluster
genesinclusters11 <- as.data.frame(genesinclusters11)

write.table(genesinclusters11, "geneclusters11.csv", sep=",")

numberofgenespercluster <- genesinclusters11 %>%
  group_by(genesinclusters11) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster11.csv")

pdf("mfuzz_sharedtimepointstraps_11clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()

#12 clusters
clusters <- mfuzz(genes.s, c=12, m=m)

genesinclusters12 <- clusters$cluster
genesinclusters12 <- as.data.frame(genesinclusters12)

write.table(genesinclusters12, "geneclusters12.csv", sep=",")

numberofgenespercluster <- genesinclusters12 %>%
  group_by(genesinclusters12) %>%
  summarise(count=n())

numberofgenespercluster

write_csv(numberofgenespercluster, "numberofgenespercluster12.csv")

pdf("mfuzz_sharedtimepointstraps_12clusters_06032022.pdf")

mfuzz.plot2(genes.s, cl=clusters, mfrow=c(4,4), time.labels=c(1,2,3),
            centre = TRUE, centre.lwd = 2,
            centre.col = "white", x11=FALSE)
dev.off()

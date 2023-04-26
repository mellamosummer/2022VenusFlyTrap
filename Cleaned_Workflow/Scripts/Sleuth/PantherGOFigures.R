library(tidyverse)
library(patchwork)
library(ggplot2)
options(scipen = 999)

##############################
#GENE COEXPRESSION GO PLOTS
##############################

##############################
#GO terms Module 3
##############################

setwd("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/7_GeneCoexpressionAnalysis/AllTrapsAnalysis/NetworkModules/PantherGOClassification")


Module3GoBP<- read.delim(header=FALSE, "Module3/pantherBiologicalProcess.txt")
Module3GoBP <- Module3GoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
 filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

Module3GoBPGGplot <-ggplot(Module3GoBP) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess))

Module3GoCC<- read.delim(header=FALSE, "Module3/pantherCellularComponent.txt")
Module3GoCC <- Module3GoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

Module3GoCCGGplot <- ggplot(Module3GoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

Module3GoMF<- read.delim(header=FALSE, "Module3/pantherMolecularFunction.txt")
Module3GoMF <- Module3GoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
 filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

Module3GoMFGGplot <- ggplot(Module3GoMF) +
  geom_point(aes(y = MolecularFunction , x = NumberOfGenes ))

Module3Go <- wrap_plots(Module3GoBPGGplot, Module3GoCCGGplot, Module3GoMFGGplot, nrow = 1)

ggsave("Module3GO_unclassifiedremoved.pdf", height = 5, width = 20)

Module3GoPC<- read.delim(header=FALSE, "Module3/pantherProteinClass.txt")
Module3GoPC <- Module3GoPC %>% select(V1:V5) %>% 
  rename(number = V1, ProteinClass = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!ProteinClass=='No PANTHER category is assigned (UNCLASSIFIED)')

Module3GoPCGGplot <- ggplot(Module3GoPC) +
  geom_point(aes(y = ProteinClass , x = NumberOfGenes ))

ggsave("Module3ProteinClass.pdf")

Module3GoP<- read.delim(header=FALSE, "Module3/pantherPathway.txt")
Module3GoP <- Module3GoP %>% select(V1:V5) %>% 
  rename(number = V1, Pathway = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!Pathway=='No PANTHER category is assigned (UNCLASSIFIED)')

Module3GoPGGplot <- ggplot(Module3GoP) +
  geom_point(aes(y = Pathway , x = NumberOfGenes))

ggsave("Module3Pathway.pdf", height = 20)

##############################
#GO terms Module 4
##############################

Module4GoBP<- read.delim(header=FALSE, "Module4/pantherBiologicalProcess.txt")
Module4GoBP <- Module4GoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

Module4GoBPGGplot <-ggplot(Module4GoBP) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess))

Module4GoCC<- read.delim(header=FALSE, "Module4/pantherCellularComponent.txt")
Module4GoCC <- Module4GoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

Module4GoCCGGplot <- ggplot(Module4GoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

Module4GoMF<- read.delim(header=FALSE, "Module4/pantherMolecularFunction.txt")
Module4GoMF <- Module4GoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

Module4GoMFGGplot <- ggplot(Module4GoMF) +
  geom_point(aes(y = MolecularFunction , x = NumberOfGenes ))

Module4Go <- wrap_plots(Module4GoBPGGplot, Module4GoCCGGplot, Module4GoMFGGplot, nrow = 1)

ggsave("Module4GO_unclassifiedremoved.pdf", height = 5, width = 20)

Module4GoPC<- read.delim(header=FALSE, "Module4/pantherProteinClass.txt")
Module4GoPC <- Module4GoPC %>% select(V1:V5) %>% 
  rename(number = V1, ProteinClass = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!ProteinClass=='No PANTHER category is assigned (UNCLASSIFIED)')

Module4GoPCGGplot <- ggplot(Module4GoPC) +
  geom_point(aes(y = ProteinClass , x = NumberOfGenes ))

ggsave("Module4ProteinClass.pdf")

Module4GoP<- read.delim(header=FALSE, "Module4/pantherPathway.txt")
Module4GoP <- Module4GoP %>% select(V1:V5) %>% 
  rename(number = V1, Pathway = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!Pathway=='No PANTHER category is assigned (UNCLASSIFIED)')

Module4GoPGGplot <- ggplot(Module4GoP) +
  geom_point(aes(y = Pathway , x = NumberOfGenes))

ggsave("Module4Pathway.pdf", height = 20)

##############################
#GO terms Module 10
##############################

Module10GoBP<- read.delim(header=FALSE, "Module10/pantherBiologicalProcess.txt")
Module10GoBP <- Module10GoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

Module10GoBPGGplot <-ggplot(Module10GoBP) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess))

Module10GoCC<- read.delim(header=FALSE, "Module10/pantherCellularComponent.txt")
Module10GoCC <- Module10GoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

Module10GoCCGGplot <- ggplot(Module10GoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

Module10GoMF<- read.delim(header=FALSE, "Module10/pantherMolecularFunction.txt")
Module10GoMF <- Module10GoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

Module10GoMFGGplot <- ggplot(Module10GoMF) +
  geom_point(aes(y = MolecularFunction , x = NumberOfGenes ))

Module10Go <- wrap_plots(Module10GoBPGGplot, Module10GoCCGGplot, Module10GoMFGGplot, nrow = 1)

ggsave("Module10GO_unclassifiedremoved.pdf", height = 5, width = 20)

Module10GoPC<- read.delim(header=FALSE, "Module10/pantherProteinClass.txt")
Module10GoPC <- Module10GoPC %>% select(V1:V5) %>% 
  rename(number = V1, ProteinClass = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!ProteinClass=='No PANTHER category is assigned (UNCLASSIFIED)')

Module10GoPCGGplot <- ggplot(Module10GoPC) +
  geom_point(aes(y = ProteinClass , x = NumberOfGenes ))

ggsave("Module10ProteinClass.pdf")

Module10GoP<- read.delim(header=FALSE, "Module10/pantherPathway.txt")
Module10GoP <- Module10GoP %>% select(V1:V5) %>% 
  rename(number = V1, Pathway = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) %>% 
  filter(!Pathway=='No PANTHER category is assigned (UNCLASSIFIED)')

Module10GoPGGplot <- ggplot(Module10GoP) +
  geom_point(aes(y = Pathway , x = NumberOfGenes))

ggsave("Module10Pathway.pdf", height = 20)


##############################
#SLEUTH DEG GO PLOTS
##############################

setwd("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/3_sleuth")

##############################
#GO terms traps vs. petioles
##############################
SleuthPetiolesVTrapsGoBP<- read.delim(header=FALSE, "PetiolesVsTraps/PantherGOClassification/pantherBiologicalProcess.txt")
SleuthPetiolesVTrapsGoBP <- SleuthPetiolesVTrapsGoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
 # filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

SleuthPetiolesVTrapsGoBPGGplot <-ggplot(SleuthPetiolesVTrapsGoBP) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess))

SleuthPetiolesVTrapsGoCC<- read.delim(header=FALSE, "PetiolesVsTraps/PantherGOClassification/pantherCellularComponent.txt")
SleuthPetiolesVTrapsGoCC <- SleuthPetiolesVTrapsGoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
  #filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

SleuthPetiolesVTrapsGoCCGGplot <- ggplot(SleuthPetiolesVTrapsGoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

SleuthPetiolesVTrapsGoMF<- read.delim(header=FALSE, "PetiolesVsTraps/PantherGOClassification/pantherMolecularFunction.txt")
SleuthPetiolesVTrapsGoMF <- SleuthPetiolesVTrapsGoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5)# %>% 
 # filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

SleuthPetiolesVTrapsGoMFGGplot <- ggplot(SleuthPetiolesVTrapsGoMF) +
  geom_point(aes(x = NumberOfGenes, y = MolecularFunction))

SleuthPetiolesVTrapsGo <- wrap_plots(SleuthPetiolesVTrapsGoBPGGplot, SleuthPetiolesVTrapsGoCCGGplot, SleuthPetiolesVTrapsGoMFGGplot, nrow = 1)

ggsave("SleuthSleuthPetiolesVTrapsGOPlots.pdf", height = 3, width = 15)

#####################
#GO terms 1 hr- all
#####################

Sleuth1hrGoBP<- read.delim(header=FALSE, "PreyVsNoPrey_1hr/PantherGOClassification/pantherBiologicalProcess.txt")
Sleuth1hrGoBP <- Sleuth1hrGoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
  #filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth1hrGoBP$Time <- c("1hr")

Sleuth1hrGoBPGGplot <-ggplot(Sleuth1hrGoBP) +
  geom_point(aes(x = PerecentGe, y = BiologicalProcess))

Sleuth1hrGoCC<- read.delim(header=FALSE, "PreyVsNoPrey_1hr/PantherGOClassification/pantherCellularComponent.txt")
Sleuth1hrGoCC <- Sleuth1hrGoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
  #filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth1hrGoCC$Time <- c("1hr")

Sleuth1hrGoCCGGplot <- ggplot(Sleuth1hrGoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

Sleuth1hrGoMF<- read.delim(header=FALSE, "PreyVsNoPrey_1hr/PantherGOClassification/pantherMolecularFunction.txt")
Sleuth1hrGoMF <- Sleuth1hrGoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
 # filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth1hrGoMF$Time <- c("1hr")

Sleuth1hrGoMFGGplot <- ggplot(Sleuth1hrGoMF) +
  geom_point(aes(x = NumberOfGenes, y = MolecularFunction))

Sleuth1hrGo <- wrap_plots(Sleuth1hrGoBPGGplot, Sleuth1hrGoCCGGplot, Sleuth1hrGoMFGGplot, nrow = 1)

ggsave("Sleuth1HrGOPlots.pdf", height = 3, width = 15)

#####################
#GO terms 24 hr- all
#####################

Sleuth24hrGoBP<- read.delim(header=FALSE, "PreyVsNoPrey_24hr/PantherGOClassification/pantherBiologicalProcess.txt")
Sleuth24hrGoBP <- Sleuth24hrGoBP %>% select(V1:V5) %>% 
  rename(number = V1, BiologicalProcess = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5)# %>% 
 # filter(!BiologicalProcess=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth24hrGoBP$Time <- c("24hr")

Sleuth24hrGoBPGGplot <-ggplot(Sleuth24hrGoBP) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess))

Sleuth24hrGoCC<- read.delim(header=FALSE, "PreyVsNoPrey_24hr/PantherGOClassification/pantherCellularComponent.txt")
Sleuth24hrGoCC <- Sleuth24hrGoCC %>% select(V1:V5) %>% 
  rename(number = V1, CellularComponent = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
#  filter(!CellularComponent=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth24hrGoCC$Time <- c("24hr")

Sleuth24hrGoCCGGplot <- ggplot(Sleuth24hrGoCC) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent))

Sleuth24hrGoMF<- read.delim(header=FALSE, "PreyVsNoPrey_24hr/PantherGOClassification/pantherMolecularFunction.txt")
Sleuth24hrGoMF <- Sleuth24hrGoMF %>% select(V1:V5) %>% 
  rename(number = V1, MolecularFunction = V2, NumberOfGenes = V3, PercentGeneHitAgainstTotalGenes = V4, PercentGeneHitAgainstTotalHits = V5) #%>% 
#  filter(!MolecularFunction=='No PANTHER category is assigned (UNCLASSIFIED)')

Sleuth24hrGoMF$Time <- c("24hr")

Sleuth24hrGoMFGGplot <- ggplot(Sleuth24hrGoMF) +
  geom_point(aes(x = NumberOfGenes, y = MolecularFunction))

Sleuth24hrGo <- wrap_plots(Sleuth24hrGoBPGGplot, Sleuth24hrGoCCGGplot, Sleuth24hrGoMFGGplot, nrow = 1)

ggsave("Sleuth24HrGOPlots.pdf", height = 3, width = 15)


####################
#DEG # plot
####################

HigherExpressionIn <- c("PreyTrigger","MechanicalTrigger", "PreyTrigger","MechanicalTrigger")
Timepoint<- c("1hr","1hr", "24hr","24hr")
NumberOfGenes <- c(103, 71, 151, 2)

DEGs <- data.frame(Timepoint, Direction, HigherExpressionIn)


ggplot(data=DEGs, aes(x=Timepoint, y=NumberOfGenes, fill=HigherExpressionIn)) +
  geom_bar(stat="identity")

ggsave("SleuthPairwiseDEGs.pdf")


##########

SleuthBPCombined<- rbind(Sleuth1hrGoBP, Sleuth24hrGoBP)

SleuthMFCombined<- rbind(Sleuth1hrGoMF, Sleuth24hrGoMF)

SleuthCCCombined<- rbind(Sleuth1hrGoCC, Sleuth24hrGoCC)

SleuthCombinedBPGGplot <- ggplot(SleuthBPCombined) +
  geom_point(aes(x = NumberOfGenes, y = BiologicalProcess, color=Time), alpha = 0.7)

SleuthCombinedBPGGplot

SleuthCombinedCCGGplot <- ggplot(SleuthCCCombined) +
  geom_point(aes(x = NumberOfGenes, y = CellularComponent, color=Time), alpha = 0.7)

SleuthCombinedCCGGplot

SleuthCombinedMFGGplot <- ggplot(SleuthMFCombined) +
  geom_point(aes(x = NumberOfGenes, y = MolecularFunction, color=Time), alpha = 0.7)

SleuthCombinedMFGGplot

SleuthTimePointGOComparison <- wrap_plots(SleuthCombinedBPGGplot, SleuthCombinedMFGGplot, SleuthCombinedCCGGplot, nrow = 1)
SleuthTimePointGOComparison
ggsave("SleuthCombinedGOplot.pdf", height=5, width=25)


####################
# DEGs in Modules Plot
####################

Sleuth1HrListModules<- read_csv("PreyVsNoPrey_1hr/DEGcsvs/SleuthList1HrAnnotated_modulesandorthogroups.csv")
Sleuth1HrListModules$module <- (as.character(Sleuth1HrListModules$module))

ggplot(data=Sleuth1HrListModules, aes(x=factor(module, level=c('2', '3', '4', '5','6','7','9','10','16','NA')), fill=DE_Direction)) +
        geom_histogram(stat="count") +
        xlab("Modules") + ylab("NumberOfDEGs")

ggsave("Sleuth1hrDEGsinModules.pdf")

Sleuth24HrListModules<- read_csv("PreyVsNoPrey_24hr/DEGcsvs/SleuthList24HrAnnotated_modulesandorthogroups.csv")
Sleuth24HrListModules$module <- (as.character(Sleuth24HrListModules$module))

ggplot(data=Sleuth24HrListModules, aes(x=factor(module, level=c('2', '3', '4','6','10','NA')), fill=DE_Direction)) +
  geom_histogram(stat="count") +
  xlab("Modules") + ylab("NumberOfDEGs")

ggsave("Sleuth24hrDEGsinModules.pdf")

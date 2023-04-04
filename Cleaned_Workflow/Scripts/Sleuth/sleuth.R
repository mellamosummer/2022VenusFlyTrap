install.packages("sleuth")
library("sleuth")

##################################
# All Sample Check
##################################

#set input and output dirs
datapath <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/quant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/sleuth"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)


#create a sample-to-condition metadata object
  
sample <- c('JMGQ',
               'JMGR',
               'JMGS',
               'JMGT',
               'JMGU',
               'JMGW',
               'JMGX',
               'JMGY',
               'JMGZ',
               'JMHA',
               'JMHB',
               'JMHC',
               'JMHD',
               'JMHE',
               'JMHF',
               'JMHG',
               'JMHH',
               'JMHI',
               'JMHJ',
               'JMHK',
               'JMHL',
               'JMHM',
               'JMHN',
               'JMHP',
               'JMHQ',
               'JMHR',
               'JMHS',
               'JMHT',
               'JMHU',
               'JMHW',
               'JMHX',
               'JMHY',
               'JMHZ',
               'JMIA',
               'JMIB',
               'JMIC',
               'JMID',
               'JMIE',
               'JMIF',
               'JMIG',
               'JMIH',
               'JMII',
               'JMIJ',
               'JMIK',
               'JMIL',
               'JMIM',
               'JMIN',
               'JMIP')

kallisto_dirs <- file.path(datapath, sample) #create vector of paths to kallisto output directories

condition <- c('LeafTrap_NoPrey_0060',
               'LeafTrap_NoPrey_0005',
               'Petiole_NoPrey_0000',
               'LeafTrap_NoPrey_0005',
               'LeafTrap_Prey_0060',
               'LeafTrap_Prey_0060',
               'LeafTrap_Prey_0060',
               'LeafTrap_Prey_0180',
               'LeafTrap_Prey_0180',
               'LeafTrap_Prey_0180',
               'LeafTrap_Prey_0180',
               'LeafTrap_Prey_4320',
               'LeafTrap_Prey_4320',
               'LeafTrap_Prey_4320',
               'LeafTrap_NoPrey_0000',
               'LeafTrap_Prey_4320',
               'LeafTrap_Prey_1440',
               'LeafTrap_Prey_1440',
               'LeafTrap_Prey_1440',
               'LeafTrap_Prey_1440',
               'LeafTrap_Prey_2880',
               'LeafTrap_Prey_2880',
               'LeafTrap_Prey_2880',
               'LeafTrap_Prey_2880',
               'LeafTrap_NoPrey_0000',
               'LeafTrap_NoPrey_0005',
               'LeafTrap_NoPrey_0005',
               'LeafTrap_NoPrey_1440',
               'LeafTrap_NoPrey_1440',
               'LeafTrap_NoPrey_1440',
               'LeafTrap_NoPrey_1440',
               'LeafTrap_NoPrey_0060',
               'LeafTrap_NoPrey_0060',
               'LeafTrap_NoPrey_0060',
               'LeafTrap_NoPrey_0000',
               'LeafTrap_Prey_0720',
               'LeafTrap_Prey_0720',
               'LeafTrap_Prey_0720',
               'LeafTrap_Prey_0720',
               'Petiole_NoPrey_0000',
               'Petiole_NoPrey_0000',
               'Petiole_NoPrey_0000',
               'LeafTrap_NoPrey_0000',
               'LeafTrap_Prey_0005',
               'LeafTrap_Prey_0005',
               'LeafTrap_Prey_0005',
               'LeafTrap_Prey_0005',
               'LeafTrap_Prey_0060')

#create vector of treatments in same order as sample IDs
samples_to_conditions <- data.frame(sample,condition) #create dataframe associating treatments to sample IDs
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs) #add kallisto output paths to dataframe

# check that directories and metadata object are OK
print(kallisto_dirs)
print(samples_to_conditions)

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#plot PCA
pdf(file="SleuthPCAAllSamples.pdf")
plot_pca(sleuth_object, color_by='condition')
dev.off()

pdf(file="SleuthPCALoadingsAllSamples.pdf")
plot_loadings(sleuth_object)
dev.off()

##################################
# LeafTrap Vs. Petiole
##################################

#create a sample-to-condition metadata object
sample <- c('JMGQ','JMGR','JMGS','JMGT','JMGU','JMGW','JMGX','JMGY','JMGZ','JMHA','JMHB','JMHC','JMHD','JMHE','JMHF','JMHG','JMHH','JMHI','JMHJ','JMHK','JMHL','JMHM','JMHN','JMHP','JMHQ','JMHR','JMHS','JMHT','JMHU','JMHW','JMHX','JMHY','JMHZ','JMIA','JMIB','JMIC','JMID','JMIE','JMIF','JMIG','JMIH','JMII','JMIJ','JMIK','JMIL','JMIM','JMIN','JMIP')

#create vector of sample IDs
kallisto_dirs <- file.path(datapath, sample) #create vector of paths to kallisto output directories

condition <- c('LeafTrap','LeafTrap','Petiole','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','Petiole','Petiole','Petiole','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap','LeafTrap')

#create vector of treatments in same order as sample IDs
samples_to_conditions <- data.frame(sample,condition) #create dataframe associating treatments to sample IDs
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs) #add kallisto output paths to dataframe

# check that directories and metadata object are OK
# print(kallisto_dirs)
# print(samples_to_conditions)

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#summarize the sleuth results for DE genes with q-values < 0.05
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#print the summary table for DE genes with q-values < 0.05
write.csv(x = sleuth_significant, file = "petioleVstraps_sleuth_q_0.05.csv", row.names = FALSE)

#visualize results for the 10 most significant DE genes
pdf(file="SleuthResults.pdf")
for(i in sleuth_significant$target_id[1:1623]) {
  p1 <- plot_bootstrap(sleuth_object, i, units = "tpm", color_by = "condition")
  print(p1)
}
dev.off()

#plot PCA
pdf(file="SleuthPCALeafVsTraps.pdf")
plot_pca(sleuth_object, color_by='condition')
dev.off()

pdf(file="SleuthPCALoadingsLeafVsTraps.pdf")
plot_loadings(sleuth_object)
dev.off()

rm(list = ls())

##################################
# Prey vs. No Prey 5 min
##################################

#NO DEGS!

##################################
# Prey vs. No Prey 60 min
##################################

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant/"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/sleuth"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)


#create a sample-to-condition metadata object
sample <- c('JMGU','JMGW','JMGX','JMIP','JMGQ','JMIA','JMHY','JMHZ')
#create vector of sample IDs
kallisto_dirs <- file.path(datapath, sample) #create vector of paths to kallisto output directories
condition <- c('Prey','Prey','Prey','Prey','NoPrey','NoPrey','NoPrey','NoPrey')

 #create vector of treatments in same order as sample IDs
samples_to_conditions <- data.frame(sample,condition) #create dataframe associating treatments to sample IDs
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs) #add kallisto output paths to dataframe

# check that directories and metadata object are OK
# print(kallisto_dirs)
# print(samples_to_conditions)

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#plot PCA
pdf(file="SleuthPCAPreyVsNoPrey1hr.pdf")
plot_pca(sleuth_object, color_by='condition')
dev.off()

pdf(file="SleuthPCALoadingsPreyVsNoPrey1hr.pdf")
plot_loadings(sleuth_object)
dev.off()

#summarize the sleuth results for DE genes with q-values < 0.05
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#print the summary table for DE genes with q-values < 0.05
write.csv(x = sleuth_significant, file = "SleuthPreyNoPrey60minResults_q_0.05.csv", row.names = FALSE)

#visualize results for the 10 most significant DE genes
pdf(file="SleuthPreyNoPrey60minResults.pdf")
for(i in sleuth_significant$target_id[1:174]) {
  p1 <- plot_bootstrap(sleuth_object, i, units = "tpm", color_by = "condition")
  print(p1)
}
dev.off()

rm(list = ls())

##################################
# Prey vs. No Prey 1440 min
##################################

#set input and output dirs
datapath <- "/scratch/srb67793/2022VenusFlyTrap/kallisto/quant/"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/scratch/srb67793/2022VenusFlyTrap/sleuth"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

#create a sample-to-condition metadata object
sample <- c('JMHH','JMHI','JMHJ','JMHK','JMHT','JMHU','JMHW','JMHX')

#create vector of sample IDs
kallisto_dirs <- file.path(datapath, sample) #create vector of paths to kallisto output directories
condition <- c('Prey','Prey','Prey','Prey','NoPrey','NoPrey','NoPrey','NoPrey')

 #create vector of treatments in same order as sample IDs
samples_to_conditions <- data.frame(sample,condition) #create dataframe associating treatments to sample IDs
samples_to_conditions <- dplyr::mutate(samples_to_conditions, path = kallisto_dirs) #add kallisto output paths to dataframe

# check that directories and metadata object are OK
# print(kallisto_dirs)
# print(samples_to_conditions)

# read data into sleuth_object, import bootstrap summaries and TPMs, and perform normalization/filtering steps
sleuth_object <- sleuth_prep(samples_to_conditions, extra_bootstrap_summary = TRUE, read_bootstrap_tpm=TRUE)

# estimate parameters for the full linear model that includes the conditions as factors
sleuth_object <- sleuth_fit(sleuth_object, ~condition, 'full')

# estimate parameters for the reduced linear model that assumes equal transcript abundances in both conditions
sleuth_object <- sleuth_fit(sleuth_object, ~1, 'reduced')

# perform likelihood ratio test to identify transcripts whose fit is significantly better under full model relative to reduced model
sleuth_object <- sleuth_lrt(sleuth_object, 'reduced', 'full')

# check that sleuth object is OK
models(sleuth_object)

#plot PCA
pdf(file="SleuthPCAPreyVsNoPrey24hr.pdf")
plot_pca(sleuth_object, color_by='condition')
dev.off()

pdf(file="SleuthPCALoadingsPreyVsNoPrey24hr.pdf")
plot_loadings(sleuth_object)
dev.off()

#summarize the sleuth results for DE genes with q-values < 0.05
sleuth_table <- sleuth_results(sleuth_object, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#print the summary table for DE genes with q-values < 0.05
write.csv(x = sleuth_significant, file = "SleuthPreyNoPrey1440min_q_0.05.csv", row.names = FALSE)

#visualize results for the 10 most significant DE genes
pdf(file="SleuthPreyNoPrey1440minResults.pdf")
for(i in sleuth_significant$target_id[1:151]) {
  p1 <- plot_bootstrap(sleuth_object, i, units = "tpm", color_by = "condition")
  print(p1)
}
dev.off()

#Quit R
quit(save="no")


SleuthList1Hr<- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/sleuth/PreyVsNoPrey_1hr/SleuthPreyNoPrey60minResults_q_0.05.csv")

SleuthList1HrAnnotated <- SleuthList1Hr %>% 
  left_join(funct_anno, by = "target_id") %>% select(target_id,AT_gene_id,annotation)

SleuthList24Hr<- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/sleuth/PreyVsNoPrey_24hr/SleuthPreyNoPrey1440min_q_0.05.csv")

SleuthList24HrAnnotated <- SleuthList24Hr %>% 
  left_join(funct_anno, by = "target_id") %>% select(target_id,AT_gene_id,annotation)

TrapsVsPetioles<- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/sleuth/PetiolesVsTraps/SleuthPetioleVstraps_q_0.05.csv")

TrapsVsPetiolesAnnotated <- TrapsVsPetioles %>% 
  left_join(funct_anno, by = "target_id") %>% select(target_id,AT_gene_id,annotation)

write.csv(x = SleuthList1HrAnnotated, file = "SleuthList1HrAnnotated.csv")

write.csv(x = SleuthList24HrAnnotated, file = "SleuthList24HrAnnotated.csv")

write.csv(x = TrapsVsPetiolesAnnotated, file = "SleuthTrapsVsPetiolesAnnotated.csv")

SleuthPreyamylase1hr <- SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('amylase', ignore_case=TRUE)))

SleuthPreychitinase1hr <- SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('chitinase', ignore_case=TRUE)))

SleuthPreyendopeptidase1hr <- SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('endopeptidase', ignore_case=TRUE)))

SleuthPreyEsterases1hr <- SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('esterase', ignore_case=TRUE)))

SleuthPreyglucanase1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucanase', ignore_case=TRUE)))

SleuthPreyGlucosidases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucosidase', ignore_case=TRUE)))

SleuthPreyglucosaminidases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucosaminidase', ignore_case=TRUE)))

SleuthPreyLipases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('lipase', ignore_case=TRUE)))

SleuthPreyNucleases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('nuclease', ignore_case=TRUE)))

SleuthPreyPeroxidases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('peroxidase', ignore_case=TRUE)))

SleuthPreyPhosphatases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphatase', ignore_case=TRUE)))

SleuthPreyphosphoaminidase1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphoaminidase', ignore_case=TRUE)))

SleuthPreyproteases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('protease', ignore_case=TRUE)))

SleuthPreyribonuclease1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('ribonuclease', ignore_case=TRUE)))

SleuthPreyureases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('urease', ignore_case=TRUE)))

SleuthPreyXylosidases1hr<-SleuthList1HrAnnotated %>% 
  filter(str_detect(annotation, fixed('xylosidase', ignore_case=TRUE)))

SleuthPreyamylase24hr <- SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('amylase', ignore_case=TRUE)))

SleuthPreychitinase24hr<- SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('chitinase', ignore_case=TRUE)))

SleuthPreyendopeptidase24hr <- SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('endopeptidase', ignore_case=TRUE)))

SleuthPreyEsterases24hr <- SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('esterase', ignore_case=TRUE)))

SleuthPreyglucanase24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucanase', ignore_case=TRUE)))

SleuthPreyGlucosidases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucosidase', ignore_case=TRUE)))

SleuthPreyglucosaminidases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('glucosaminidase', ignore_case=TRUE)))

SleuthPreyLipases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('lipase', ignore_case=TRUE)))

SleuthPreyNucleases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('nuclease', ignore_case=TRUE)))

SleuthPreyPeroxidases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('peroxidase', ignore_case=TRUE)))

SleuthPreyPhosphatases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphatase', ignore_case=TRUE)))

SleuthPreyphosphoaminidase24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphoaminidase', ignore_case=TRUE)))

SleuthPreyproteases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('protease', ignore_case=TRUE)))

SleuthPreyribonuclease24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('ribonuclease', ignore_case=TRUE)))

SleuthPreyureases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('urease', ignore_case=TRUE)))

SleuthPreyXylosidases24hr<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('xylosidase', ignore_case=TRUE)))

TrapsVsPetiolesamylase <- TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('amylase', ignore_case=TRUE)))

TrapsVsPetioleschitinase<- TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('chitinase', ignore_case=TRUE)))

TrapsVsPetiolesendopeptidase <- TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('endopeptidase', ignore_case=TRUE)))

TrapsVsPetiolesEsterases <- TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('esterase', ignore_case=TRUE)))

TrapsVsPetiolesglucanase<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('glucanase', ignore_case=TRUE)))

TrapsVsPetiolesGlucosidases<-SleuthList24HrAnnotated %>% 
  filter(str_detect(annotation, fixed('TrapsVsPetiolesAnnotated', ignore_case=TRUE)))

TrapsVsPetiolesglucosaminidases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('glucosaminidase', ignore_case=TRUE)))

TrapsVsPetiolesLipases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('lipase', ignore_case=TRUE)))

TrapsVsPetiolesNucleases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('nuclease', ignore_case=TRUE)))

TrapsVsPetioleseroxidases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('peroxidase', ignore_case=TRUE)))

TrapsVsPetiolesPhosphatases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphatase', ignore_case=TRUE)))

STrapsVsPetiolesphosphoaminidase<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('phosphoaminidase', ignore_case=TRUE)))

TrapsVsPetiolesproteases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('protease', ignore_case=TRUE)))

TrapsVsPetiolesribonuclease<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('ribonuclease', ignore_case=TRUE)))

TrapsVsPetiolesureases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('urease', ignore_case=TRUE)))

TrapsVsPetiolesXylosidases<-TrapsVsPetiolesAnnotated %>% 
  filter(str_detect(annotation, fixed('xylosidase', ignore_case=TRUE)))

write.csv(x = SleuthPreyamylase1hr, file = "SleuthPreyamylase1hr.csv")
write.csv(x = SleuthPreychitinase1hr, file = "SleuthPreychitinase1hr.csv")
write.csv(x = SleuthPreyendopeptidase1hr, file = "SleuthPreyendopeptidase1hr.csv")
write.csv(x = SleuthPreyEsterases1hr, file = "SleuthPreyEsterases1hr.csv")
write.csv(x = SleuthPreyglucanase1hr, file = "SleuthPreyglucanase1hr.csv")
write.csv(x = SleuthPreyGlucosidases1hr, file = "SleuthPreyGlucosidases1hr.csv")
write.csv(x = SleuthPreyglucosaminidases1hr, file = "SleuthPreyglucosaminidases1hr.csv")
write.csv(x = SleuthPreyLipases1hr, file = "SleuthPreyLipases1hr.csv")
write.csv(x = SleuthPreyNucleases1hr, file = "SleuthPreyNucleases1hr.csv")
write.csv(x = SleuthPreyPeroxidases1hr, file = "SleuthPreyPeroxidases1hr.csv")
write.csv(x = SleuthPreyPhosphatases1hr, file = "SleuthPreyPhosphatases1hr.csv")
write.csv(x = SleuthPreyphosphoaminidase1hr, file = "SleuthPreyphosphoaminidase1hr.csv")
write.csv(x = SleuthPreyproteases1hr, file = "SleuthPreyproteases1hr.csv")
write.csv(x = SleuthPreyribonuclease1hr, file = "SleuthPreyribonuclease1hr.csv")
write.csv(x = SleuthPreyureases1hr, file = "SleuthPreyureases1hr.csv")
write.csv(x = SleuthPreyXylosidases1hr, file = "SleuthPreyXylosidases1hr.csv")

write.csv(x = TrapsVsPetiolesamylase, file = "TrapsVsPetiolesamylase.csv")
write.csv(x = TrapsVsPetiolesendopeptidase, file = "TrapsVsPetiolesendopeptidase.csv")
write.csv(x = TrapsVsPetioleseroxidases, file = "TrapsVsPetiolesPeroxidases.csv")
write.csv(x = TrapsVsPetiolesEsterases, file = "TrapsVsPetiolesEsterases.csv")
write.csv(x = TrapsVsPetiolesLipases, file = "TrapsVsPetiolesLipases.csv")
write.csv(x = TrapsVsPetiolesNucleases, file = "TrapsVsPetiolesNucleases.csv")
write.csv(x = TrapsVsPetiolesPhosphatases, file = "TrapsVsPetiolesPhosphatases.csv")
write.csv(x = TrapsVsPetiolesproteases, file = "TrapsVsPetiolesproteases.csv")
write.csv(x = TrapsVsPetiolesribonuclease, file = "TrapsVsPetiolesribonuclease.csv")



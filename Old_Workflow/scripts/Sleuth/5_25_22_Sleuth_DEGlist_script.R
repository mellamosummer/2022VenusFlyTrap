library(dplyr)
library(tidyverse)
library(plyr)

#prepping the data-------------------
#getting the data from the CSV file and cleaning
setwd("/Users/srblanco/Desktop/UGA/Dissertation/Leebens Mack Rotation/VFTAnalysis2/Kallisto_quant_results")

#determine a matrix of samples informations
samples <- read.table("samples.txt", header = FALSE)
samples


#building design matrix with path, treatment (prey/noprey), time

samples <- samples %>%
  add_column(condition = "Value")

samples <- samples %>%
  add_column(path = "Value")

samples <- samples %>%
  rename(sample = V1)

samples <- rename(samples, time = path)

samples

samples[1:4,2] <- "Petiole"
samples[5:20, 2] <- "No_prey"
samples[21:48, 2] <- "Prey"
samples


#changed for each time treatment
samples[1:4, 3] <- "0000"
samples[5:8, 3] <- "0005"
samples[9:12, 3] <- "0060"
samples[13:16, 3] <- "0000"
samples[17:20, 3] <- "1440"
samples[21:24, 3] <- "0005"
samples[25:28, 3] <- "0060"
samples[29:32, 3] <- "0180"
samples[33:36, 3] <- "0720"
samples[37:40, 3] <- "1440"
samples[41:44, 3] <- "2880"
samples[45:48, 3] <- "4320"
samples


samples<- as.vector(samples)
samples

#generating pairwise comparisons for sleuth

pairwise_time_5 <- filter(samples, time == "0005")
pairwise_time_5

pairwise_time_60 <- filter(samples, time == "0060")
pairwise_time_60

pairwise_time_60_outlier_removed <- pairwise_time_60[-c(1), ]
pairwise_time_60_outlier_removed

pairwise_time_1440 <- filter(samples, time == "1440")
pairwise_time_1440

pairwise_time_1440_outlier_removed <- pairwise_time_1440[-c(8), ]
pairwise_time_1440_outlier_removed

#generating individual condition dataframes to use to build pairwise comparison with petiole tissue
samples

traps_noprey_time0 <- samples %>% 
  filter(condition == "No_prey") %>% 
  filter(time == "0000")
traps_noprey_time0

traps_noprey_time5 <- samples %>% 
  filter(condition == "No_prey") %>% 
  filter(time == "0005")
traps_noprey_time5

traps_noprey_time60 <- samples %>% 
  filter(condition == "No_prey") %>% 
  filter(time == "0060")
traps_noprey_time60

traps_noprey_time1440 <- samples %>% 
  filter(condition == "No_prey") %>% 
  filter(time == "1440")
traps_noprey_time1440

traps_prey_time5 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "0005")
traps_prey_time5

traps_prey_time60 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "0060")
traps_prey_time60

traps_prey_time180 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "0180")
traps_prey_time180

traps_prey_time720 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "0720")
traps_prey_time720

traps_prey_time1440 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "1440")
traps_prey_time1440

traps_prey_time2880 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "2880")
traps_prey_time2880

traps_prey_time4320 <- samples %>% 
  filter(condition == "Prey") %>% 
  filter(time == "4320")
traps_prey_time4320

petioles <- samples %>% 
  filter(condition == "Petiole") 
petioles

#create pairwise petiole comparisons

petiole_vs_traps_noprey_time0 <- bind_rows(petioles, traps_noprey_time0)
petiole_vs_traps_noprey_time0

petiole_vs_traps_noprey_time5 <- bind_rows(petioles, traps_noprey_time5)
petiole_vs_traps_noprey_time5

petiole_vs_traps_noprey_time60 <- bind_rows(petioles, traps_noprey_time60)
petiole_vs_traps_noprey_time60

petiole_vs_traps_noprey_time1440 <- bind_rows(petioles, traps_noprey_time1440)
petiole_vs_traps_noprey_time1440

petiole_vs_traps_prey_time5 <- bind_rows(petioles, traps_prey_time5)
petiole_vs_traps_prey_time5

petiole_vs_traps_prey_time60 <- bind_rows(petioles, traps_prey_time60)
petiole_vs_traps_prey_time60

petiole_vs_traps_prey_time180 <- bind_rows(petioles, traps_prey_time180)
petiole_vs_traps_prey_time180

petiole_vs_traps_prey_time720 <- bind_rows(petioles, traps_prey_time720)
petiole_vs_traps_prey_time720

petiole_vs_traps_prey_time1440 <- bind_rows(petioles, traps_prey_time1440)
petiole_vs_traps_prey_time1440

petiole_vs_traps_prey_time2880 <- bind_rows(petioles, traps_prey_time2880)
petiole_vs_traps_prey_time2880

petiole_vs_traps_prey_time4320 <- bind_rows(petioles, traps_prey_time4320)
petiole_vs_traps_prey_time4320

#prey with traps only
traps_prey <- samples %>% 
  filter(condition == "Prey")

traps_prey

traps_prey_no_outliers <- traps_prey[-c(20, 22), ]

traps_noprey <- samples %>% 
  filter(condition == "No_prey")

traps_noprey
traps_noprey_outliers_removed <- traps_noprey[-c(11, 14), ]


#beginning of sleuth code , generates data frame with sample, condition, and path---------------
sample_id <- subset(pairwise_time_1440_outlier_removed, select = c(condition, sample, sample)) 
###I changed "pairwise_time_1440_outlier_removed" for each of the pairwise comparisons to generate the output files.
sample_id



#rename one of the samples to file path 
colnames(sample_id)[3] <- "path"
sample_id

#Sleuth Prep and model fitting
library("sleuth")
so <- sleuth_prep(sample_id, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

#Plot single gene for all conditions
plot_bootstrap(so, " Dm_00000118-RA", units = "est_counts", color_by = "condition")

#Plot PCA
jpeg("PCA_traps_noprey_outliers_removed.jpg")
plot_pca(so, color_by = 'condition')
dev.off()

#Plot group density
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates), "sample"), offset = 1)
#Shiny Packages                                                                                                    
install.packages("shiny")
library("shiny")
install.packages("gridExtra")
library("gridExtra")

#interactive data explorer thing
sleuth_live(so)

write_csv(sleuth_significant, "5_24_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time1440_outlier_removed.csv")

#combining DEG lists
#getting the data from the CSV file and cleaning
setwd("/Users/srblanco/Desktop/UGA/Leebens Mack Rotation/VFTAnalysis2/Sleuth_Analysis")

#prey vs no prey traps pairwise comparison
pairwiseanalysis_traps_vs_no_traps_time_1440 <- read_csv("5_24_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time1440.csv")
pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed <- read_csv("5_24_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time1440_outlier_removed.csv")

#prey vs no prey traps pairwise comparison
pairwiseanalysis_traps_vs_no_traps_time_60 <- read_csv("5_24_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time60.csv")
pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed <- read_csv("5_24_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time60_outlier_removed.csv")

#traps vs petiole  pairwise comparison
#prey
pairwiseanalysis_petiole_vs_traps_prey_time5 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time5.csv")
pairwiseanalysis_petiole_vs_traps_prey_time60 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time60.csv")
pairwiseanalysis_petiole_vs_traps_prey_time180 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time180.csv")
pairwiseanalysis_petiole_vs_traps_prey_time720 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time720.csv")
pairwiseanalysis_petiole_vs_traps_prey_time1440 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time1440.csv")
pairwiseanalysis_petiole_vs_traps_prey_time2880 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time2880.csv")
pairwiseanalysis_petiole_vs_traps_prey_time4320 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_prey_time4320.csv")
#noprey
pairwiseanalysis_petiole_vs_traps_noprey_time0 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time0.csv")
pairwiseanalysis_petiole_vs_traps_noprey_time5 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time5.csv")
pairwiseanalysis_petiole_vs_traps_noprey_time60 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time60.csv")
pairwiseanalysis_petiole_vs_traps_noprey_time1440 <- read_csv("2_3_2022_sleuth_pairwiseanalysis_petiole_vs_traps_noprey_time1440.csv")

#makes csvs into lists with gene names

list_pairwiseanalysis_traps_vs_no_traps_time_1440 <- subset(pairwiseanalysis_traps_vs_no_traps_time_1440, select = c(target_id))
list_pairwiseanalysis_traps_vs_no_traps_time_1440 <- dplyr::rename(list_pairwiseanalysis_traps_vs_no_traps_time_1440, list_pairwiseanalysis_traps_vs_no_traps_time_1440 = target_id)
write.table(list_pairwiseanalysis_traps_vs_no_traps_time_1440, file = "list_pairwiseanalysis_traps_vs_no_traps_time_1440.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed <- subset(pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed, select = c(target_id))
list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed <- dplyr::rename(list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed, list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed = target_id)
write.table(list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed, file = "list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_traps_vs_no_traps_time_60 <- subset(pairwiseanalysis_traps_vs_no_traps_time_60, select = c(target_id))
list_pairwiseanalysis_traps_vs_no_traps_time_60 <- dplyr::rename(list_pairwiseanalysis_traps_vs_no_traps_time_60, list_pairwiseanalysis_traps_vs_no_traps_time_60 = target_id)
write.table(list_pairwiseanalysis_traps_vs_no_traps_time_60, file = "list_pairwiseanalysis_traps_vs_no_traps_time_60.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed <- subset(pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed, select = c(target_id))
list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed <- dplyr::rename(list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed, list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed = target_id)
write.table(list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed, file = "list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time5 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time5, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time5 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time5, list_pairwiseanalysis_petiole_vs_traps_prey_time5 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time5, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time5.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time60 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time60, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time60 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time60, list_pairwiseanalysis_petiole_vs_traps_prey_time60 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time60, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time60.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time180 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time180, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time180 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time180, list_pairwiseanalysis_petiole_vs_traps_prey_time180 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time180, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time180.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time720 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time720, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time720 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time720, list_pairwiseanalysis_petiole_vs_traps_prey_time720 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time720, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time720.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time1440 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time1440, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time1440 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time1440, list_pairwiseanalysis_petiole_vs_traps_prey_time1440 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time1440, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time1440.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_petiole_vs_traps_prey_time2880 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time2880, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time2880 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time2880, list_pairwiseanalysis_petiole_vs_traps_prey_time2880 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time2880, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time2880.lst", col.names = F, row.names = F, quote = FALSE)

list_pairwiseanalysis_petiole_vs_traps_prey_time4320 <- subset(pairwiseanalysis_petiole_vs_traps_prey_time4320, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_prey_time4320 <- rename(list_pairwiseanalysis_petiole_vs_traps_prey_time4320, list_pairwiseanalysis_petiole_vs_traps_prey_time4320 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_prey_time4320, file = "list_pairwiseanalysis_petiole_vs_traps_prey_time4320.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_petiole_vs_traps_noprey_time0 <- subset(pairwiseanalysis_petiole_vs_traps_noprey_time0, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_noprey_time0 <- rename(list_pairwiseanalysis_petiole_vs_traps_noprey_time0, list_pairwiseanalysis_petiole_vs_traps_noprey_time0 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_noprey_time0, file = "list_pairwiseanalysis_petiole_vs_traps_noprey_time0.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_petiole_vs_traps_noprey_time5 <- subset(pairwiseanalysis_petiole_vs_traps_noprey_time5, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_noprey_time5 <- rename(list_pairwiseanalysis_petiole_vs_traps_noprey_time5, list_pairwiseanalysis_petiole_vs_traps_noprey_time5 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_noprey_time5, file = "list_pairwiseanalysis_petiole_vs_traps_noprey_time5.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_petiole_vs_traps_noprey_time60 <- subset(pairwiseanalysis_petiole_vs_traps_noprey_time60, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_noprey_time60 <- rename(list_pairwiseanalysis_petiole_vs_traps_noprey_time60, list_pairwiseanalysis_petiole_vs_traps_noprey_time60 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_noprey_time60, file = "list_pairwiseanalysis_petiole_vs_traps_noprey_time60.lst", col.names = F, row.names = F, quote = FALSE)


list_pairwiseanalysis_petiole_vs_traps_noprey_time1440 <- subset(pairwiseanalysis_petiole_vs_traps_noprey_time1440, select = c(target_id))
list_pairwiseanalysis_petiole_vs_traps_noprey_time1440 <- rename(list_pairwiseanalysis_petiole_vs_traps_noprey_time1440, list_pairwiseanalysis_petiole_vs_traps_noprey_time1440 = target_id)
write.table(list_pairwiseanalysis_petiole_vs_traps_noprey_time1440, file = "list_pairwiseanalysis_petiole_vs_traps_noprey_time1440.lst", col.names = F, row.names = F, quote = FALSE)


#combine gene lists
install.packages("gdata")
library(gdata)
combined_pairwise_petioles_vs_traps_prey <- cbindX(list_pairwiseanalysis_petiole_vs_traps_prey_time5, list_pairwiseanalysis_petiole_vs_traps_prey_time60, list_pairwiseanalysis_petiole_vs_traps_prey_time180, list_pairwiseanalysis_petiole_vs_traps_prey_time720, list_pairwiseanalysis_petiole_vs_traps_prey_time1440, list_pairwiseanalysis_petiole_vs_traps_prey_time2880, list_pairwiseanalysis_petiole_vs_traps_prey_time4320)
combined_pairwise_petioles_vs_traps_prey

combined_pairwise_petioles_vs_traps_noprey<- cbindX(list_pairwiseanalysis_petiole_vs_traps_noprey_time0, list_pairwiseanalysis_petiole_vs_traps_noprey_time5, list_pairwiseanalysis_petiole_vs_traps_noprey_time60, list_pairwiseanalysis_petiole_vs_traps_noprey_time1440)
combined_pairwise_petioles_vs_traps_noprey

combined_pairwise_traps <-cbindX(list_pairwiseanalysis_traps_vs_no_traps_time_60, list_pairwiseanalysis_traps_vs_no_traps_time_60_outlier_removed, list_pairwiseanalysis_traps_vs_no_traps_time_1440, list_pairwiseanalysis_traps_vs_no_traps_time_1440_outlier_removed)
combined_pairwise_traps

combined_petioles_vs_alltraps_prey_and_no_prey <- cbindX(combined_pairwise_petioles_vs_traps_prey, combined_pairwise_petioles_vs_traps_noprey)
combined_petioles_vs_alltraps_prey_and_no_prey

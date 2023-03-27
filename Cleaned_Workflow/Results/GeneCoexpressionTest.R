library(tidyverse)
library(igraph)
library(ggraph)

library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)

set.seed(666)

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
datapath <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/quant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

##################################
#  PREY TIME SERIES
##################################

##################################
#  BUILD COUNTS MATRIX
##################################

#Create a text file with the names of the Kallisto folders and assign as "sample"
sample <-c('JMGS',
           'JMIG',
           'JMIH',
           'JMII',
           'JMGR',
           'JMGT',
           'JMHR',
           'JMHS',
           'JMGQ',
           'JMHY',
           'JMHZ',
           'JMIA',
           'JMHT',
           'JMHU',
           'JMHW',
           'JMHX',
           'JMHF',
           'JMHQ',
           'JMIB',
           'JMIJ',
           'JMIK',
           'JMIL',
           'JMIM',
           'JMIN',
           'JMGU',
           'JMGW',
           'JMGX',
           'JMIP',
           'JMGY',
           'JMGZ',
           'JMHA',
           'JMHB',
           'JMIC',
           'JMID',
           'JMIE',
           'JMIF',
           'JMHH',
           'JMHI',
           'JMHJ',
           'JMHK',
           'JMHL',
           'JMHM',
           'JMHN',
           'JMHP',
           'JMHC',
           'JMHD',
           'JMHE',
           'JMHG')

#Create file path for abundance.tsv file in each of the folders that are specified in "samples"
files <- file.path(datapath, sample,"abundance.tsv")

#extracts counts data and builds matrix with geneID and counts
txi <- tximport(files, type="kallisto", txOut = TRUE, geneIdCol = TRUE, abundanceCol = TRUE, countsCol = FALSE)
head(txi)

#create variable with just the tpm
tpm <- txi$abundance

#rename "sample1", "sample2", and so on to thea actual names of the samples
colnames(tpm) <- sample

write.csv(x = tpm, file = "Blanco_tpm_all.csv")

Exp_table <- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/Blanco_tpm_all.csv", col_types = cols())
Metadata <-  read_excel("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/metadata.xlsx")

Exp_table_long <- Exp_table %>% 
  rename(gene_ID = `...1`) %>% 
  pivot_longer(cols = !gene_ID, names_to = "LibraryName", values_to = "tpm") %>% 
  mutate(logTPM = log10(tpm + 1)) 

Exp_table_log_wide <- Exp_table_long %>% 
  select(gene_ID, LibraryName, logTPM) %>% 
  pivot_wider(names_from = LibraryName, values_from = logTPM)

head(Exp_table_log_wide)
dim(Exp_table_log_wide)

my_pca <- prcomp(t(Exp_table_log_wide[, -1]))
pc_importance <- as.data.frame(t(summary(my_pca)$importance))
head(pc_importance, 20)

PCA_coord <- my_pca$x[, 1:10] %>% 
  as.data.frame() %>% 
  mutate(LibraryName = row.names(.)) %>% 
  full_join(Metadata %>% 
              select(LibraryName, Tissue, Time, Treatment, SampleName), by = "LibraryName")

PCA_coord

PCA_by_time<- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Time), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_time

PCA_by_tissue<- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Tissue), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_tissue

PCA_by_treatment<- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_treatment

wrap_plots(PCA_by_tissue, PCA_by_treatment, PCA_by_time, nrow = 1)
ggsave("PCA1_2.png", height = 3.5, width = 15, bg = "white")
ggsave("PCA1_2.png", height = 3.5, width = 15, bg = "white")


PCA_by_time<- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = Time), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_time

PCA_by_tissue<- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = Tissue), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_tissue

PCA_by_treatment<- PCA_coord %>% 
  ggplot(aes(x = PC2, y = PC3)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(n = 9, "Set1")) +
  labs(x = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC3 (", pc_importance[3, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_treatment

wrap_plots(PCA_by_tissue, PCA_by_treatment, PCA_by_time, nrow = 1)
ggsave("PCA2_3.png", height = 3.5, width = 15, bg = "white")
ggsave("PCA2_3.png", height = 3.5, width = 15, bg = "white")

Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord %>% 
              select(LibraryName, SampleName, Tissue, Treatment, Time), 
            by = "LibraryName") %>% 
  filter(Treatment == "Prey") %>% 
  group_by(gene_ID, SampleName, Time) %>% 
  summarise(mean.logTPM = mean(logTPM)) %>% 
  ungroup()  

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(gene_ID) %>% 
  mutate(z.score = (mean.logTPM - mean(mean.logTPM))/sd(mean.logTPM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)

high_var_genes <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var, 0.667))

head(high_var_genes)
dim(high_var_genes)

Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% high_var_genes$gene_ID)

head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  group_by(gene_ID) %>% 
  count() %>% 
  nrow()

z_score_wide <- Exp_table_long_averaged_z_high_var %>% 
  select(gene_ID, SampleName, z.score) %>% 
  pivot_wider(names_from = SampleName, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
head(z_score_wide)

cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)

number_of_Treatment_Time <- ncol(z_score_wide) - 1
number_of_Treatment_Time

PCA_coord %>% 
  filter(Tissue == "Trap" | Treatment == "Prey") %>% 
  group_by(Time) %>% 
  count() %>% 
  nrow()

cor_matrix_upper_tri <- cor_matrix
cor_matrix_upper_tri[lower.tri(cor_matrix_upper_tri)] <- NA

edge_table <- cor_matrix_upper_tri %>% 
  as.data.frame() %>% 
  mutate(from = row.names(cor_matrix)) %>% 
  pivot_longer(cols = !from, names_to = "to", values_to = "r") %>% 
  filter(is.na(r) == F) %>% 
  filter(from != to) %>% 
  mutate(t = r*sqrt((number_of_Treatment_Time-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_Treatment_Time-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_Treatment_Time-2, lower.tail = T)
  )) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr")) 

head(edge_table)

edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.05) %>% 
  slice_min(order_by = abs(r), n = 10)

edge_table %>% 
  filter(r > 0) %>% 
  filter(FDR < 0.01) %>% 
  slice_min(order_by = abs(r), n = 10)

edge_table %>% 
  slice_sample(n = 20000) %>% 
  ggplot(aes(x = r)) +
  geom_histogram(color = "white", bins = 100) +
  geom_vline(xintercept = 0.7, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

ggsave("../Results/r_histogram.png", height = 3.5, width = 5, bg = "white")

edge_table_select <- edge_table %>% 
  filter(r >= 0.7)

node_table <- data.frame(
  gene_ID = c(edge_table_select$from, edge_table_select$to) %>% unique()
)

my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F
)

modules <- cluster_leiden(my_network, resolution_parameter = 2, 
                          objective_function = "modularity")

optimize_resolution <- function(network, resolution){
  modules = network %>% 
    cluster_leiden(resolution_parameter = resolution,
                   objective_function = "modularity")
  
  parsed_modules = data.frame(
    gene_ID = names(membership(modules)),
    module = as.vector(membership(modules)) 
  )
  
  num_module_5 = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    nrow() %>% 
    as.numeric()
  
  num_genes_contained = parsed_modules %>% 
    group_by(module) %>% 
    count() %>% 
    arrange(-n) %>% 
    filter(n >= 5) %>% 
    ungroup() %>% 
    summarise(sum = sum(n)) %>% 
    as.numeric()
  
  c(num_module_5, num_genes_contained)
  
}

optimization_results <- purrr::map_dfc(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  t() %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() %>% 
  rename(num_module = V1,
         num_contained_gene = V2)

Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module)) +
  geom_line(size = 1.1, alpha = 0.8, color = "dodgerblue2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. modules\nw/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

Optimize_num_gene <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_contained_gene)) +
  geom_line(size = 1.1, alpha = 0.8, color = "violetred2") +
  geom_point(size = 3, alpha = 0.7) +
  geom_vline(xintercept = 2, size = 1, linetype = 4) +
  labs(x = "resolution parameter",
       y = "num. genes in\nmodules w/ >=5 genes") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

wrap_plots(Optimize_num_module, Optimize_num_gene, nrow = 2)

ggsave("../Results/Optimize_resolution.svg", height = 5, width = 3.2, bg ="white")
ggsave("../Results/Optimize_resolution.png", height = 5, width = 3.2, bg ="white")

my_network_modules <- data.frame(
  gene_ID = names(membership(modules)),
  module = as.vector(membership(modules)) 
) %>% 
  inner_join(node_table, by = "gene_ID")

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  ungroup() %>% 
  summarise(sum = sum(n))

module_5 <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_modules <- my_network_modules %>% 
  filter(module %in% module_5$module)

head(my_network_modules)

Exp_table_long_averaged_z_high_var_modules <- Exp_table_long_averaged_z_high_var %>% 
  inner_join(my_network_modules, by = "gene_ID")

head(Exp_table_long_averaged_z_high_var_modules)

modules_mean_z <- Exp_table_long_averaged_z_high_var_modules %>% 
  group_by(module, dev_stage, tissue, `Sample Name`) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)

module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp

module_line_plot <- Exp_table_long_averaged_z_high_var_modules %>% 
  mutate(order_x = case_when(
    str_detect(dev_stage, "5") ~ 1,
    str_detect(dev_stage, "10") ~ 2,
    str_detect(dev_stage, "20") ~ 3,
    str_detect(dev_stage, "30") ~ 4,
    str_detect(dev_stage, "MG") ~ 5,
    str_detect(dev_stage, "Br") ~ 6,
    str_detect(dev_stage, "Pk") ~ 7,
    str_detect(dev_stage, "LR") ~ 8,
    str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  filter(module == "1" |
           module == "5") %>% 
  ggplot(aes(x = dev_stage, y = z.score)) +
  facet_grid(module ~ tissue) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "1" |
               module == "5") %>% 
      mutate(order_x = case_when(
        str_detect(dev_stage, "5") ~ 1,
        str_detect(dev_stage, "10") ~ 2,
        str_detect(dev_stage, "20") ~ 3,
        str_detect(dev_stage, "30") ~ 4,
        str_detect(dev_stage, "MG") ~ 5,
        str_detect(dev_stage, "Br") ~ 6,
        str_detect(dev_stage, "Pk") ~ 7,
        str_detect(dev_stage, "LR") ~ 8,
        str_detect(dev_stage, "RR") ~ 9
      )) %>% 
      mutate(dev_stage = reorder(dev_stage, order_x)),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = NULL,
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_lines_color_strip <- expand.grid(
  Tissue = unique(Metadata$Treatment),
  Time = unique(Metadata$Time), 
  stringsAsFactors = F
) %>% 
  filter(dev_stage != "Anthesis") %>% 
  filter(str_detect(tissue, "epider|chyma|Vasc") == F) %>% 
  mutate(order_x = case_when(
    str_detect(dev_stage, "5") ~ 1,
    str_detect(dev_stage, "10") ~ 2,
    str_detect(dev_stage, "20") ~ 3,
    str_detect(dev_stage, "30") ~ 4,
    str_detect(dev_stage, "MG") ~ 5,
    str_detect(dev_stage, "Br") ~ 6,
    str_detect(dev_stage, "Pk") ~ 7,
    str_detect(dev_stage, "LR") ~ 8,
    str_detect(dev_stage, "RR") ~ 9
  )) %>% 
  mutate(stage = case_when(
    str_detect(dev_stage, "MG|Br|Pk") ~ str_sub(dev_stage, start = 1, end = 2),
    T ~ dev_stage
  )) %>% 
  mutate(stage = factor(stage, levels = c(
    "5 DPA",
    "10 DPA",
    "20 DPA",
    "30 DPA",
    "MG",
    "Br",
    "Pk",
    "LR",
    "RR"
  ))) %>% 
  mutate(dev_stage = reorder(dev_stage, order_x)) %>% 
  ggplot(aes(x = dev_stage, y = 1)) +
  facet_grid(. ~ tissue) +
  geom_tile(aes(fill = stage)) +
  scale_fill_manual(values = viridis(9, option = "D")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_blank(),
    text = element_text(size = 14),
    panel.spacing = unit(1, "lines")
  )

wrap_plots(module_line_plot, module_lines_color_strip,
           nrow = 2, heights = c(1, 0.08))

ggsave("../Results/module_line_plots.svg", height = 4, width = 8.2, bg = "white")
ggsave("../Results/module_line_plots.png", height = 4, width = 8.2, bg = "white")

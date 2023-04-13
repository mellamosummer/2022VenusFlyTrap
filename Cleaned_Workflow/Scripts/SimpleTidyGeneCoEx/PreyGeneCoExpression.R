##################################
#  INSTALL PACKAGES
##################################

library(tidyverse)
library(igraph)
library(ggraph)
library(readxl)
library(patchwork)
library(RColorBrewer)
library(viridis)
library("tximport")

set.seed(666)

##################################
#  SET DIRECTORIES
##################################

#set input and output dirs
datapath <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/2_KallistoQuant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

##################################
#  BUILD TPM MATRIX
##################################

#Create a text file with the names of the Kallisto folders and assign as "sample"
sample <-c('JMHF',
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

#create variable with just the tpm
tpm <- txi$abundance

#rename "sample1", "sample2", and so on to thea actual names of the samples
colnames(tpm) <- sample

head(tpm)
dim(tpm)

Exp_table <- tpm
Exp_table <- as.data.frame(Exp_table)

write.csv(x = tpm, file = "Blanco_tpm_all.csv")

Exp_table <- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/7_GeneCoexpressionAnalysis/Data/Blanco_tpm_all.csv")
Metadata <-  read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/7_GeneCoexpressionAnalysis/Data/Metadata.csv", col_types = cols())

#check

Metadata %>% 
  group_by(Time) %>% 
  count()

Metadata %>% 
  group_by(Tissue) %>% 
  count()

Metadata %>% 
  group_by(Treatment) %>% 
  count()

##################################
# PCA
##################################

Exp_table <- Exp_table %>% 
  rownames_to_column(var = "gene_ID")


Exp_table_long<- Exp_table %>% 
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
              select(LibraryName, Time, Treatment), by = "LibraryName")

PCA_coord <- PCA_coord %>% filter(Treatment == "prey" | Time == "0")
PCA_coord$Time <- as.numeric(PCA_coord$Time)
PCA_coord

PCA_by_time<- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Time), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )


PCA_by_time

wrap_plots(PCA_by_time, PCA_by_treatment, nrow = 1)

ggsave("PCA_Prey_by_time_tissue.png", height = 3.5, width = 8.5, bg = "white")

##################################
# Average up the reps
##################################

Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord, by = "LibraryName") %>% 
  group_by(gene_ID, Time) %>% 
  summarise(mean.logTPM = mean(logTPM)) %>% 
  ungroup()

head(Exp_table_long_averaged)
##################################
# z score
##################################

Exp_table_long_averaged_z <- Exp_table_long_averaged %>% 
  group_by(gene_ID) %>% 
  mutate(z.score = (mean.logTPM - mean(mean.logTPM))/sd(mean.logTPM)) %>% 
  ungroup()

head(Exp_table_long_averaged_z)


##################################
# Gene selection -- variance
##################################

Expressed_genes <- Exp_table_long %>% 
  filter(tpm > 5) %>% 
  group_by(gene_ID) %>% 
  count() %>% 
  filter(n >= 2)

dim(Expressed_genes)

high_var_genes <- Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% Expressed_genes$gene_ID) %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  filter(var > quantile(var,0.5))

dim(high_var_genes)

high_var_genes %>% 
  filter(gene_ID %in% Baits$gene_ID)

all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 

Exp_table_long_averaged_z_high_var <- Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% high_var_genes$gene_ID)

head(Exp_table_long_averaged_z_high_var)

Exp_table_long_averaged_z_high_var %>% 
  group_by(gene_ID) %>% 
  count() %>% 
  nrow()


##################################
# Gene selection -- F stats
##################################

compute_F <- function(data, gene){
  anova_table = lm(logTPM ~ Time, data %>% 
                     filter(gene_ID == gene)) %>% 
    anova()
  
  cbind(anova_table$`F value`[1], anova_table$`Pr(>F)`[1]) %>% 
    as.data.frame()
}

ANOVA_results <- purrr::map_dfr(
  .x = Expressed_genes$gene_ID,
  .f = compute_F,
  data =  Exp_table_long %>% 
    full_join(PCA_coord, by = "LibraryName")
) %>% 
  cbind(gene_ID = Expressed_genes$gene_ID) %>% 
  as.data.frame() %>% 
  rename(
    F_stat = V1,
    p.value = V2
  ) %>% 
  mutate(FDR = p.adjust(p.value, method = "fdr"))

F_high_var_comparison <- ANOVA_results %>% 
  left_join(all_var_and_ranks, by = "gene_ID") %>%
  mutate(high_F = case_when(
    F_stat >= 2 ~ "high F",
    T ~ "low F"
  )) %>% 
  mutate(high_var = case_when(
    gene_ID %in% high_var_genes$gene_ID ~ "high var",
    T ~ "low var"
  )) %>% 
  mutate(type = case_when(
    high_F == "high F" & 
      high_var == "high var" ~ "both",
    high_F == "high F" &
      high_var == "low var" ~ "high F only",
    high_F == "low F" & 
      high_var == "high var" ~ "high var only",
    T ~ "neither"
  ))

F_var_scatter <- F_high_var_comparison %>% 
  ggplot(aes(x = F_stat, y = sqrt(var))) +
  geom_point(aes(color = type), alpha = 0.5) +
  scale_color_manual(values = brewer.pal(n = 8, "Set2")) +
  labs(x = "F stat",
       y = "sd(log10(TPM))",
       color = NULL) +
  theme_classic() +
  theme(
    legend.position = c(0.75, 0.75),
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

F_var_scatter

F_dist <- F_high_var_comparison %>% 
  ggplot(aes(x = F_stat)) +
  geom_histogram(bins = 100, fill = brewer.pal(8, "Set2")[2]) +
  labs(x = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )


F_dist

var_dist <- F_high_var_comparison %>% 
  ggplot(aes(x = sqrt(var))) +
  geom_histogram(bins = 100, fill = brewer.pal(8, "Set2")[3]) +
  labs(x = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  ) +
  coord_flip()

var_dist

blank <- F_high_var_comparison %>% 
  ggplot(aes(x = -log10(FDR), y = sqrt(var))) +
  theme_void()

wrap_plots(
  F_dist, blank,
  F_var_scatter, var_dist,
  nrow = 2, ncol = 2, 
  heights = c(0.35, 1),
  widths = c(1, 0.35)
)

ggsave("var_F_scatter_prey.png", width = 5.5, height = 5.5)

F_high_var_comparison %>% 
  group_by(type) %>% 
  count()

F_high_var_comparison %>% 
  group_by(high_var) %>% 
  count()

F_high_var_comparison %>% 
  group_by(high_F) %>% 
  count()

high_F_only_examples <- F_high_var_comparison %>% 
  filter(type == "high F only") %>% 
  slice_max(order_by = F_stat, n = 1) %>% 
  inner_join(Exp_table_long, by = "gene_ID") %>% 
  inner_join(PCA_coord, by = "LibraryName")

high_var_only_examples <- F_high_var_comparison %>% 
  filter(type == "high var only") %>% 
  slice_max(order_by = var, n = 1) %>% 
  inner_join(Exp_table_long, by = "gene_ID") %>% 
  inner_join(PCA_coord, by = "LibraryName")

high_var_high_F_examples <- F_high_var_comparison %>% 
  filter(type == "both") %>% 
  slice_max(order_by = F_stat, n = 1) %>% 
  inner_join(Exp_table_long, by = "gene_ID") %>% 
  inner_join(PCA_coord, by = "LibraryName")

neither_examples <- F_high_var_comparison %>% 
  filter(type == "neither") %>% 
  slice_max(order_by = var, n = 1) %>% 
  inner_join(Exp_table_long, by = "gene_ID") %>% 
  inner_join(PCA_coord, by = "LibraryName")

high_F_only_graphs <- high_F_only_examples %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  geom_point(aes(fill = Treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(TPM)",
       title = paste0("high F only example\n",
                      high_F_only_examples$gene_ID)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

high_F_only_graphs

high_var_only_graphs <- high_var_only_examples %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  geom_point(aes(fill = Treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(TPM)",
       title = paste0("high var only example\n",
                      high_var_only_examples$gene_ID)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

high_var_high_F_graphs <- high_var_high_F_examples %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  geom_point(aes(fill = Treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(TPM)",
       title = paste0("high F high var example\n",
                      high_var_high_F_examples$gene_ID)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

low_var_low_F_graphs <- neither_examples %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  geom_point(aes(fill = Treatment), shape = 21, color = "white",
             size = 3, alpha = 0.8, position = position_jitter(0.2, seed = 666)) +
  scale_fill_manual(values = c("grey70", "tomato1")) +
  labs(x = "time point",
       y = "log10(TPM)",
       title = paste0("low F low var example\n",
                      neither_examples$gene_ID)) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12)
  )

wrap_plots(
  high_var_only_graphs,
  high_var_high_F_graphs,
  low_var_low_F_graphs,
  high_F_only_graphs,
  nrow = 2, guides = "collect"
) & theme(legend.position = "bottom")

high_var_or_high_F_genes <- F_high_var_comparison %>% 
  filter(type != "neither")

dim(high_var_or_high_F_genes)

Exp_table_long_averaged_z_high_var_or_high_F <- Exp_table_long_averaged_z %>% 
  filter(gene_ID %in% high_var_or_high_F_genes$gene_ID)

head(Exp_table_long_averaged_z_high_var_or_high_F)

##################################
# Gene wise correlation
##################################

z_score_wide <- Exp_table_long_averaged_z_high_var_or_high_F %>% 
  mutate(tag = paste(Time, sep = "-")) %>% 
  select(gene_ID, tag, z.score) %>% 
  pivot_wider(names_from = tag, values_from = z.score) %>% 
  as.data.frame()

row.names(z_score_wide) <- z_score_wide$gene_ID
head(z_score_wide)

cor_matrix <- cor(t(z_score_wide[, -1]))
dim(cor_matrix)

##################################
# Edge selection
##################################

number_of_time_treatment <- ncol(z_score_wide) - 1
number_of_time_treatment

PCA_coord %>% 
  group_by(Time, Treatment) %>% 
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
  mutate(t = r*sqrt((number_of_time_treatment-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_time_treatment-2, lower.tail = F),
    t <=0 ~ pt(t, df = number_of_time_treatment-2, lower.tail = T)
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
  geom_vline(xintercept = 0.75, color = "tomato1", size = 1.2) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black")
  )

##################################
# Module detection
##################################

edge_table_select <- edge_table %>% 
  filter(r >= 0.75)

dim(edge_table_select)

##################################
# Build graph object
##################################

funct_anno <- read_delim("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/5_BLAST/AllDmProteinsBLASTresults/DmProteinsTairdbBLASTconcat.txt", 
                         delim = "\t", col_names = F, col_types = cols())

funct_anno <- funct_anno %>% select(X1:X3) %>% 
  rename(gene_ID = X1, AT_gene_id = X2, annotation = X3)

node_table <- data.frame(
  gene_ID = c(edge_table_select$from, edge_table_select$to) %>% unique()
) %>% 
  left_join(funct_anno, by = "gene_ID")


head(node_table)
dim(node_table)

my_network <- graph_from_data_frame(
  edge_table_select,
  vertices = node_table,
  directed = F
)


##################################
# optimize clustering resolution
##################################

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
  
  cbind(num_module_5, num_genes_contained) %>% 
    as.data.frame()
  
}

optimization_results <- purrr::map_dfr(
  .x = seq(from = 0.25, to = 5, by = 0.25),
  .f = optimize_resolution, 
  network = my_network
) %>% 
  cbind(
    resolution = seq(from = 0.25, to = 5, by = 0.25)
  ) %>% 
  as.data.frame() 

head(optimization_results)

Optimize_num_module <- optimization_results %>% 
  ggplot(aes(x = resolution, y = num_module_5)) +
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
  ggplot(aes(x = resolution, y = num_genes_contained)) +
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

ggsave("module_optimization_prey.png", height = 10, width = , bg = "white")

##################################
# Graph based clustering
##################################

modules <- cluster_leiden(my_network, resolution_parameter = 2, 
                          objective_function = "modularity")

head(modules)
my_network_modules <- data.frame(
  gene_ID = names(membership(modules)),
  module = as.vector(membership(modules)) 
)

my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5) %>% 
  print(n = 21)

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

my_network_modules <- my_network_modules %>% 
  left_join(funct_anno, by = "gene_ID")

#serine endopeptidases
my_network_modules %>% 
  filter(gene_ID == "Dm_00000773-RA" |
           gene_ID == "Dm_00001789-RA" |
           gene_ID == "Dm_00004253-RA" |
           gene_ID == "Dm_00015236-RA")

#gene_ID == "Dm_00000773-RA" module 13
#gene_ID == "Dm_00001789-RA" module 6
#gene_ID == "Dm_00004253-RA" not in module
#gene_ID == "Dm_00015236-RA" module 9 


#cysteine endopeptidases

my_network_modules %>% 
  filter(gene_ID == "Dm_00010272-RA" |
           gene_ID == "Dm_00017138-RA") #not in modules

#metallo endopeptidases

my_network_modules %>% 
  filter(gene_ID == "Dm_00012406-RA" |
           gene_ID == "Dm_00002722-RA" |
           gene_ID == "Dm_00001146-RA" |
           gene_ID == "Dm_00002723-RA" |
           gene_ID == "Dm_00012588-RA")


#gene_ID == "Dm_00012406-RA" not in modules
#gene_ID == "Dm_00002722-RA" not in modules
#gene_ID == "Dm_00001146-RA" not in modules
#gene_ID == "Dm_00002723-RA" module 2
#gene_ID == "Dm_00012588-RA" module 2

#chitinases

my_network_modules %>% 
  filter(gene_ID == "Dm_00013373-RA" |
           gene_ID == "Dm_00019697-RA" |
           gene_ID == "Dm_00004094-RA" |
           gene_ID == "Dm_00000939-RA" |
           gene_ID == "Dm_00020493-RA")

#all in module 3

#lipases

my_network_modules %>% 
  filter(gene_ID == "Dm_00003116-RA" |
           gene_ID == "Dm_00002298-RA" |
           gene_ID == "Dm_00004094-RA" |
           gene_ID == "Dm_00002448-RA" |
           gene_ID == "Dm_00010368-RA")

#gene_ID == "Dm_00003116-RA" module 3
#gene_ID == "Dm_00002298-RA" module 4
#gene_ID == "Dm_00004094-RA" not in modules
#gene_ID == "Dm_00002448-RA" module 4
#gene_ID == "Dm_00010368-RA" not in modules

##################################
# Module Treatment correspondence
##################################

Exp_table_long_averaged_z_high_var_or_high_F_modules <- Exp_table_long_averaged_z_high_var_or_high_F %>% 
  inner_join(my_network_modules, by = "gene_ID")

head(Exp_table_long_averaged_z_high_var_or_high_F_modules)

modules_mean_z <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  group_by(module, Time) %>% 
  summarise(mean.z = mean(z.score)) %>% 
  ungroup()

head(modules_mean_z)

module_peak_exp <- modules_mean_z %>% 
  group_by(module) %>% 
  slice_max(order_by = mean.z, n = 1)

module_peak_exp

##################################
# Module QC
##################################

Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z, aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

ggsave("module_line_plots_ordered_prey.png", height = 10, width = 40, bg = "white")

module_line_plot_3 <-Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "3") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "3"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_3

ggsave("module3_line_plot_prey.png", height = 10, width = , bg = "white")

module_line_plot_4 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "4") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "4"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_4

ggsave("module_line_plot_4prey.png", height = 10, width = 8, bg = "white")

module_line_plot_5 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "5") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "5"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_5

ggsave("module_line_plot_5prey.png", height = 10, width = 8, bg = "white")

module_line_plot_6 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "6") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "6"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_6

ggsave("module_line_plot_6.png", height = 10, width = 8, bg = "white")


module_line_plot_7<- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "7") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "7"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_7

ggsave("module_line_plot_7prey.png", height = 10, width = 8, bg = "white")

module_line_plot_8 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "8") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "8"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_8

ggsave("module_line_plot_8.png", height = 10, width = 8, bg = "white")

module_line_plot_12 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "12") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "12"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_12

ggsave("module_line_plot_12prey.png", height = 10, width = 8, bg = "white")

module_line_plot_14 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "14") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "14"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_14

ggsave("module_line_plot_14prey.png", height = 10, width = 8, bg = "white")

module_line_plot_19 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "19") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "19"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_19

ggsave("module_line_plot_19prey.png", height = 10, width = 8, bg = "white")

module_line_plot_27 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "27") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "27"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_27

ggsave("module_line_plot_27prey.png", height = 10, width = 8, bg = "white")

module_line_plot_32 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "32") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "32"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_32

ggsave("module_line_plot_32prey.png", height = 10, width = 8, bg = "white")

module_line_plot_58 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "58") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "58"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_58

ggsave("module_line_plot_58prey.png", height = 10, width = 8, bg = "white")

module_line_plot_445 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "445") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "445"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_445

ggsave("module_line_plot_445prey.png", height = 10, width = 8, bg = "white")

module_line_plot_463 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "463") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "463"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_463

ggsave("module_line_plot_463prey.png", height = 10, width = 8, bg = "white")

module_line_plot_774 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  filter(module == "774") %>%
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") + 
  geom_line(data = modules_mean_z %>% 
              filter(module == "774"), aes(y = mean.z, group = module), size = 1.1, alpha = 0.8) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_774

ggsave("module_line_plot_774prey.png", height = 10, width = 8, bg = "white")


##################################
# Heatmap
##################################

modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, c(0.05, 0.95))

modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 2 ~ 2,
    mean.z < -2 ~ -2,
    T ~ mean.z
  ))

modules_mean_z_reordered <- modules_mean_z %>% 
  full_join(module_peak_exp %>% 
              rename(peak_time = Time) %>% 
              select(module, peak_time), by = "module") %>% 
  mutate(module = reorder(module, -peak_time))

head(modules_mean_z_reordered)

modules_mean_z_reordered %>% 
  ggplot(aes(x = as.factor(Time), y = as.factor(module))) +
  geom_tile(aes(fill = mean.z)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-2, 2),
                       breaks = c(-2, 0, 2), labels = c("< 2", "0", "2")) +
  labs(x = "Time point",
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  )

ggsave("prey_module_heatmap.png", height = 10, width = 8, bg = "white")

##################################
# Gene co-expression graphs
##################################

#chitinases

neighbors_of_bait_chitinase <- c(
  neighbors(my_network, v = "Dm_00013373-RA"),
  neighbors(my_network, v = "Dm_00019697-RA"),
  neighbors(my_network, v = "Dm_00004094-RA"), 
  neighbors(my_network, v = "Dm_00000939-RA"), 
  neighbors(my_network, v = "Dm_00020493-RA")
) %>% 
  unique()  

length(neighbors_of_bait_chitinase)

chitinase_transcription_factors <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait_chitinase)) %>% 
  filter(str_detect(annotation, "transcription factor"))

nrow(chitinase_transcription_factors)
head(chitinase_transcription_factors)

chitinase_subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait_chitinase) &
           to %in% names(neighbors_of_bait_chitinase)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

length(chitinase_subnetwork_edges)
dim(chitinase_subnetwork_edges)

chitinase_subnetwork_genes <- c(chitinase_subnetwork_edges$from, chitinase_subnetwork_edges$to) %>% unique()

chitinase_subnetwork_nodes <- node_table %>% 
  filter(gene_ID %in% chitinase_subnetwork_genes) %>% 
  left_join(my_network_modules, by = "gene_ID") %>% 
  left_join(module_peak_exp, by = "module") 

dim(chitinase_subnetwork_nodes)

chitinase_subnetwork <- graph_from_data_frame(chitinase_subnetwork_edges,
                                              vertices = chitinase_subnetwork_nodes,
                                              directed = F)

chitinase_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = Time)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 8), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("chitinase_subnetwork_graph.png", height = 10, width = 12, bg = "white")

#metalloendopeptidases

neighbors_of_bait_metallo <- c(
  neighbors(my_network, v = "Dm_00002723-RA"),
  neighbors(my_network, v = "Dm_00012588-RA")
) %>% 
  unique()  

length(neighbors_of_bait_metallo)

metallo_transcription_factors <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait_metallo)) %>% 
  filter(str_detect(annotation, "transcription factor"))

nrow(metallo_transcription_factors)
head(metallo_transcription_factors)

metallo_subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait_metallo) &
           to %in% names(neighbors_of_bait_metallo)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

length(metallo_subnetwork_edges)
dim(metallo_subnetwork_edges)

metallo_subnetwork_genes <- c(metallo_subnetwork_edges$from, metallo_subnetwork_edges$to) %>% unique()

metallo_subnetwork_nodes <- node_table %>% 
  filter(gene_ID %in% metallo_subnetwork_genes) %>% 
  left_join(my_network_modules, by = "gene_ID") %>% 
  left_join(module_peak_exp, by = "module") 

dim(metallo_subnetwork_nodes)

metallo_subnetwork <- graph_from_data_frame(metallo_subnetwork_edges,
                                            vertices = metallo_subnetwork_nodes,
                                            directed = F)

metallo_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = Time)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 8), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("metallo_subnetwork_graph.png", height = 10, width = 12, bg = "white")

#lipases
#gene_ID == "Dm_00003116-RA" module 3
#gene_ID == "Dm_00002298-RA" module 4
#gene_ID == "Dm_00002448-RA" module 4

neighbors_of_bait_lipase <- c(
  neighbors(my_network, v = "Dm_00003116-RA"),
  neighbors(my_network, v = "Dm_00002298-RA"),
  neighbors(my_network, v = "Dm_00002448-RA")
) %>% 
  unique()  

length(neighbors_of_bait_lipase)

lipase_transcription_factors <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait_lipase)) %>% 
  filter(str_detect(annotation, "transcription factor"))

nrow(lipase_transcription_factors)
head(lipase_transcription_factors)

lipase_subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait_lipase) &
           to %in% names(neighbors_of_bait_lipase)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

length(lipase_subnetwork_edges)
dim(lipase_subnetwork_edges)

lipase_subnetwork_genes <- c(lipase_subnetwork_edges$from, lipase_subnetwork_edges$to) %>% unique()

lipase_subnetwork_nodes <- node_table %>% 
  filter(gene_ID %in% lipase_subnetwork_genes) %>% 
  left_join(my_network_modules, by = "gene_ID") %>% 
  left_join(module_peak_exp, by = "module") 

dim(lipase_subnetwork_nodes)

lipase_subnetwork <- graph_from_data_frame(lipase_subnetwork_edges,
                                           vertices = lipase_subnetwork_nodes,
                                           directed = F)

lipase_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = Time)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 8), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("lipase_subnetwork_graph.png", height = 10, width = 12, bg = "white")


#serine endo
#gene_ID == "Dm_00000773-RA" module 13
#gene_ID == "Dm_00001789-RA" module 6
#gene_ID == "Dm_00015236-RA" module 9 

neighbors_of_bait_serine <- c(
  neighbors(my_network, v = "Dm_00000773-RA"),
  neighbors(my_network, v = "Dm_00001789-RA"),
  neighbors(my_network, v = "Dm_00015236-RA")
) %>% 
  unique()  

length(neighbors_of_bait_serine)

serine_transcription_factors <- my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait_serine)) %>% 
  filter(str_detect(annotation, "transcription factor"))

nrow(serine_transcription_factors)
head(serine_transcription_factors)

my_network_modules %>% 
  filter(gene_ID %in% names(neighbors_of_bait_serine))

nrow(my_protein_kinase)

serine_subnetwork_edges <- edge_table_select %>% 
  filter(from %in% names(neighbors_of_bait_serine) &
           to %in% names(neighbors_of_bait_serine)) %>% 
  group_by(from) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup() %>% 
  group_by(to) %>% 
  slice_max(order_by = r, n = 5) %>% 
  ungroup()

length(serine_subnetwork_edges)
dim(serine_subnetwork_edges)

serine_subnetwork_genes <- c(serine_subnetwork_edges$from, serine_subnetwork_edges$to) %>% unique()

serine_subnetwork_nodes <- node_table %>% 
  filter(gene_ID %in% serine_subnetwork_genes) %>% 
  left_join(my_network_modules, by = "gene_ID") %>% 
  left_join(module_peak_exp, by = "module") 

dim(serine_subnetwork_nodes)

serine_subnetwork <- graph_from_data_frame(serine_subnetwork_edges,
                                           vertices = serine_subnetwork_nodes,
                                           directed = F)

serine_subnetwork %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = Time)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 8), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )

ggsave("serine_subnetwork_graph.png", height = 10, width = 12, bg = "white")

coex <- graph_from_data_frame(edge_table_select,
                              directed = F)

my_network %>% 
  ggraph(layout = "kk", circular = F) +
  geom_edge_diagonal(color = "grey70", width = 0.5, alpha = 0.5) +
  geom_node_point(alpha = 0.8, color = "grey30", shape = 21, size = 2,
                  aes(fill = Time)) + 
  scale_fill_gradientn(colors = viridis(10, option = "A")) +
  labs(fill = "Peak time point") +
  guides(size = "none",
         fill = guide_colorsteps()) +
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = "bottom",
    legend.justification = 1,
    title = element_text(size = 12)
  )


neighbors_of_bait_serinecsv <- funct_anno %>% 
  filter(gene_ID %in% names(neighbors_of_bait_serine)) %>% 
  select(gene_ID, AT_gene_id, annotation)

write_excel_csv(neighbors_of_bait_serinecsv, "serineendopeptidase_neighbors.csv", col_names = T)

neighbors_of_bait_chitinaseecsv <- funct_anno %>% 
  filter(gene_ID %in% names(neighbors_of_bait_chitinase)) %>% 
  select(gene_ID, AT_gene_id, annotation)

write_excel_csv(neighbors_of_bait_chitinaseecsv, "chitinase_neighbors.csv", col_names = T)

neighbors_of_bait_lipasecsv <- funct_anno %>% 
  filter(gene_ID %in% names(neighbors_of_bait_lipase)) %>% 
  select(gene_ID, AT_gene_id, annotation)

write_excel_csv(neighbors_of_bait_lipasecsv, "lipase_neighbors.csv", col_names = T)

neighbors_of_bait_metallocsv <- funct_anno %>% 
  filter(gene_ID %in% names(neighbors_of_bait_metallo)) %>% 
  select(gene_ID, AT_gene_id, annotation)

write_excel_csv(neighbors_of_bait_metallocsv, "metallo_neighbors.csv", col_names = T)


##################################
# Annotate modules
##################################

write.csv(x = my_network_modules, file = "network_modules_prey.csv")

Annotations <- read_tsv("BLAST/DmProteinsTairdbBLASTconcat.txt", col_names = FALSE)
head(Annotations)

Annotation_subset <- Annotations %>% 
  select(1:3)

colnames(Annotation_subset) <- c("gene_ID","ATgene_ID", "TopBLASTHit")

head(Annotation_subset)

my_network_moduless_annotated <- my_network_modules %>% 
  full_join(Annotation_subset, by = "gene_ID")

my_network_modules %>% group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_moduless_annotated %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

my_network_moduless_annotated <- my_network_moduless_annotated %>% 
  filter(!module == "NA")

my_network_moduless_annotated %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)


write.csv(x = my_network_moduless_annotated, file = "network_modules_annotated.csv")

###

my_network_modules <- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/GeneCoexpressionAnalysis/NetworkModules/network_modules_annotated.csv")


genesinmodules <- my_network_modules %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n) %>% 
  filter(n >= 5)

write.csv(x = genesinmodules, file = "genesinmodulesprey.csv")

phosphatasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('phosphatase', ignore_case=TRUE)))

phosphatasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

proteasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('protease', ignore_case=TRUE)))

proteasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

chitinasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('chitinase', ignore_case=TRUE)))

chitinasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)


glucanasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('glucanase', ignore_case=TRUE)))

glucanasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)


esterasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('esterase', ignore_case=TRUE)))

esterasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

peroxidasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('peroxidase', ignore_case=TRUE)))

peroxidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

nucleasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('nuclease', ignore_case=TRUE)))

nucleasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

glucosaminidasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('glucosaminidase', ignore_case=TRUE)))

glucosaminidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

glucosidasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('glucosidase', ignore_case=TRUE)))

glucosidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

amylasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('amylase', ignore_case=TRUE)))

amylasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

lipasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('lipase', ignore_case=TRUE)))

lipasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

ribonucleasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('ribonuclease', ignore_case=TRUE)))

ribonucleasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

phosphoamidasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('phosphoamidase', ignore_case=TRUE)))

phosphoamidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)


xylosidasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('xylosidase', ignore_case=TRUE)))

xylosidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

ureasesinnetwork <- my_network_modules %>% 
  filter(str_detect(annotation, fixed('urease', ignore_case=TRUE)))

ureasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

cysteineendopeptidasesinnetwork<- my_network_modules %>% 
  filter(str_detect(annotation, fixed('cysteine endopeptidase', ignore_case=TRUE)))

cysteineendopeptidasesinnetwork %>% 
  group_by(module) %>% 
  count() %>% 
  arrange(-n)

write.csv(x = cysteineendopeptidasesinnetwork, file = "cysteineendopeptidasesinnetworkprey.csv")
write.csv(x = ureasesinnetwork, file = "ureasesinnetworkprey.csv")
write.csv(x = xylosidasesinnetwork, file = "xylosidasesinnetworkprey.csv")
write.csv(x = phosphoamidasesinnetwork, file = "phosphoamidasesinnetworkprey.csv")
write.csv(x = ribonucleasesinnetwork, file = "ribonucleasesinnetworkprey.csv")
write.csv(x = lipasesinnetwork, file = "lipasesinnetworkprey.csv")
write.csv(x = amylasesinnetwork, file = "amylasesinnetworkprey.csv")
write.csv(x = glucosidasesinnetwork, file = "glucosidasesinnetworkprey.csv")
write.csv(x = glucosaminidasesinnetwork, file = "glucosaminidasesinnetworkprey.csv")
write.csv(x = nucleasesinnetwork, file = "nucleasesinnetworkprey.csv")
write.csv(x = peroxidasesinnetwork, file = "peroxidasesinnetworkprey.csv")
write.csv(x = esterasesinnetwork, file = "esterasesinnetworkprey.csv")
write.csv(x = glucanasesinnetwork, file = "glucanasesinnetworkprey.csv")
write.csv(x = chitinasesinnetwork, file = "chitinasesinnetworkprey.csv")
write.csv(x = proteasesinnetwork, file = "proteasesinnetworkprey.csv")
write.csv(x = phosphatasesinnetwork, file = "phosphatasesinnetwork.csv")


Exp_table_long %>% 
  filter(gene_ID %in% chitinasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = "log10(TPM)")) +
  facet_grid(gene_ID ~ ., scales = "free_y") +
  geom_point(aes(fill = 'black'), color = "black", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = Time), fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = "Time",
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    strip.background = element_blank()
  )

chitinases <- Exp_table_long %>% 
  filter(gene_ID %in% chitinasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

chitinases + scale_x_continuous(limit = c(0,4320))

ggsave("chitinases_prey.png", height = 20, width = 3, bg = "white")


xylosidases <- Exp_table_long %>% 
  filter(gene_ID %in% xylosidasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )


xylosidases + scale_x_continuous(limit = c(0,4320))

ggsave("xylosidases_graph_truncated.png", height = 12, width = 5, bg = "white")

ureases <- Exp_table_long %>% 
  filter(gene_ID %in% ureasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

ureases + scale_x_continuous(limit = c(0,4320))

ggsave("ureases_graph_truncated.png", height = 12, width = 5, bg = "white")

cysteineendopeptidases <- Exp_table_long %>% 
  filter(gene_ID %in% cysteineendopeptidasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

cysteineendopeptidases + scale_x_continuous(limit = c(0,4320))

ggsave("cysteineendopeptidases_graph_truncated.png", height = 12, width = 5, bg = "white")

ribonucleases <- Exp_table_long %>% 
  filter(gene_ID %in% ribonucleasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

ribonucleases + scale_x_continuous(limit = c(0,4320))

ggsave("ribonucleases_graph_truncated.png", height = 15, width = 3, bg = "white")

lipases <- Exp_table_long %>% 
  filter(gene_ID %in% lipasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

lipases + scale_x_continuous(limit = c(0,4320))

ggsave("lipases_graph_truncated.png", height = 47, width = 3, bg = "white")

esterases <- Exp_table_long %>% 
  filter(gene_ID %in% esterasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

esterases + scale_x_continuous(limit = c(0,4320))

ggsave("esterases_graph_truncated.png", height = 55, width = 3, bg = "white", limitsize = FALSE)


glucanases <- Exp_table_long %>% 
  filter(gene_ID %in% glucanasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )
glucanases + scale_x_continuous(limit = c(0,4320))

ggsave("glucanases_prey.png", height = 8, width = 5, bg = "white", limitsize = FALSE)

glucosidase <- Exp_table_long %>% 
  filter(gene_ID %in% glucosidasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

glucosidase + scale_x_continuous(limit = c(0,4320))

ggsave("glucosidase_graph_truncated.png", height = 20, width = 3, bg = "white", limitsize = FALSE)

nucleases <- Exp_table_long %>% 
  filter(gene_ID %in% nucleasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )
nucleases + scale_x_continuous(limit = c(0,4320))

ggsave("nucleases_graph_truncated.png", height =40, width = 3, bg = "white", limitsize = FALSE)

peroxidases <- Exp_table_long %>% 
  filter(gene_ID %in% peroxidasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

peroxidases + scale_x_continuous(limit = c(0,4320))

ggsave("peroxidases_prey.png", height =30, width = 3, bg = "white", limitsize = FALSE)

phosphatases <- Exp_table_long %>% 
  filter(gene_ID %in% phosphatasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

phosphatases + scale_x_continuous(limit = c(0,4320))

ggsave("phosphatases_prey.png", height =130, width = 3, bg = "white", limitsize = FALSE)

proteases <- Exp_table_long %>% 
  filter(gene_ID %in% proteasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

proteases + scale_x_continuous(limit = c(0,4320))

ggsave("proteases_prey.png", height =50, width = 3, bg = "white", limitsize = FALSE)

amylases <- Exp_table_long %>% 
  filter(gene_ID %in% amylasesinnetwork$gene_ID) %>% 
  inner_join(PCA_coord, by = "LibraryName") %>% 
  ggplot(aes(x = Time, y = logTPM)) +
  facet_grid(gene_ID ~., scales = "free_y") +
  geom_point(aes(fill = logTPM), color = "white", size = 2, 
             alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  stat_summary(geom = "line", aes(color = logTPM), 
               fun = mean, alpha = 0.8, size = 1.1) +
  labs(x = NULL,
       y = "log10(TPM)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    strip.background = element_blank()
  )

amylases + scale_x_continuous(limit = c(0,4320))

ggsave("amylases_prey.png", height =25, width = 3, bg = "white", limitsize = FALSE)

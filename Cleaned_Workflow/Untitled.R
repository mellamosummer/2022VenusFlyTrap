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
datapath <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/quant"  # you need to modify this line to match the path made by your BASH script
resultdir <- "/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/"   # you need to modify this line to match the path made by your BASH script
setwd(resultdir)

##################################
#  BUILD TPM MATRIX
##################################

#Create a text file with the names of the Kallisto folders and assign as "sample"
sample <-c('JMGR',
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

head(tpm)
dim(tpm)

write.csv(x = tpm, file = "Blanco_tpm_all.csv")

Exp_table <- read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/Blanco_tpm_all.csv", col_types = cols())

Metadata <-  read_csv("/Users/summerblanco/Desktop/Github/2022VenusFlyTrap/Cleaned_Workflow/Results/metadata.csv", col_types = cols())

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

PCA_coord <- PCA_coord %>% filter(Tissue == "Trap") 
PCA_coord$Time <- as.numeric(PCA_coord$Time)
PCA_coord

PCA_by_tissue

PCA_by_treatment<- PCA_coord %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Treatment), color = "grey20", shape = 21, size = 3, alpha = 0.8) +
  labs(x = paste("PC1 (", pc_importance[1, 2] %>% signif(3)*100, "% of Variance)", sep = ""), 
       y = paste("PC2 (", pc_importance[2, 2] %>% signif(3)*100, "% of Variance)", "  ", sep = ""),
       fill = NULL) +  
  theme_bw() +
  theme(
    text = element_text(size= 14),
    axis.text = element_text(color = "black")
  )

PCA_by_treatment

ggsave("PCA_by_treatment.png", height = 3.5, width = 15, bg = "white")


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

ggsave("PCA_by_time_tissue.png", height = 3.5, width = 8.5, bg = "white")

##################################
# Average up the reps
##################################

Exp_table_long_averaged <- Exp_table_long %>% 
  full_join(PCA_coord, by = "LibraryName") %>% 
  group_by(gene_ID, Time, Treatment) %>% 
  summarise(mean.logTPM = mean(logTPM)) %>% 
  ungroup()

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

all_var_and_ranks <- Exp_table_long_averaged_z %>% 
  group_by(gene_ID) %>% 
  summarise(var = var(mean.logTPM)) %>% 
  ungroup() %>% 
  mutate(rank = rank(var, ties.method = "average")) 

all_var_and_ranks %>% 
  ggplot(aes(x = var, y = rank)) +
  geom_rect( 
    xmax = max(high_var_genes$var), 
    xmin = min(high_var_genes$var),
    ymax = nrow(all_var_and_ranks),
    ymin = nrow(all_var_and_ranks) - nrow(high_var_genes),
    fill = "dodgerblue2", alpha = 0.2
  ) +
  geom_line(size = 1.1) +
  labs(y = "rank",
       x = "var(log10(FPKM))",
       caption = "Blue box = high var expressed genes.") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.caption = element_text(hjust = 0)
  )


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
  anova_table = lm(logTPM ~ Time:Treatment, data %>% 
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

ggsave("var_F_scatter.png", width = 5.5, height = 5.5)

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
  mutate(tag = paste(Time, Treatment, sep = "-")) %>% 
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

my_network <- graph_from_data_frame(
  edge_table_select,
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



##################################
# Module Treatment correspondence
##################################

Exp_table_long_averaged_z_high_var_or_high_F_modules <- Exp_table_long_averaged_z_high_var_or_high_F %>% 
  inner_join(my_network_modules, by = "gene_ID")

head(Exp_table_long_averaged_z_high_var_or_high_F_modules)

modules_mean_z <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  group_by(module, Time, Treatment) %>% 
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

module_line_plot_3 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "3") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "3") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
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

ggsave("module_line_plot_3.png", height = 10, width = 8, bg = "white")

module_line_plot_2 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "2") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "2") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_2

ggsave("module_line_plot_2.png", height = 10, width = 8, bg = "white")

module_line_plot_4 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "4") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "4") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
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

ggsave("module_line_plot_4.png", height = 10, width = 8, bg = "white")


module_line_plot_10 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "10") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "10") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_10

ggsave("module_line_plot_10.png", height = 10, width = 8, bg = "white")

module_line_plot_7 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "7") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "7") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
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

ggsave("module_line_plot_7.png", height = 10, width = 8, bg = "white")

module_line_plot_9 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "9") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "9") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_9

ggsave("module_line_plot_9.png", height = 10, width = 8, bg = "white")

module_line_plot_6 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "6") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "6") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
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

module_line_plot_16 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "16") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "16") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_16

ggsave("module_line_plot_16.png", height = 10, width = 8, bg = "white")

module_line_plot_13 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "13") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "13") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_13

ggsave("module_line_plot_13.png", height = 10, width = 8, bg = "white")

module_line_plot_5 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "5") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "5") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
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

ggsave("module_line_plot_5.png", height = 10, width = 8, bg = "white")

module_line_plot_24 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "24") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "24") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_24

ggsave("module_line_plot_24.png", height = 10, width = 8, bg = "white")

module_line_plot_168 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "168") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "168") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_168

ggsave("module_line_plot_168.png", height = 10, width = 8, bg = "white")

module_line_plot_142 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "142") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "142") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_142

ggsave("module_line_plot_142.png", height = 10, width = 8, bg = "white")

module_line_plot_460 <- Exp_table_long_averaged_z_high_var_or_high_F_modules %>% 
  mutate(treatment = factor(Treatment, levels = c("prey", "no prey"))) %>% 
  filter(module == "460") %>%  
  ggplot(aes(x = Time, y = z.score)) +
  facet_grid(Treatment ~ module) +
  geom_line(aes(group = gene_ID), alpha = 0.3, color = "grey70") +
  geom_line(
    data = modules_mean_z %>% 
      filter(module == "460") %>%  
      mutate(Treatment = factor(Treatment, levels = c("prey", "no prey"))),
    aes(y = mean.z, group = module), 
    size = 1.1, alpha = 0.8
  ) +
  labs(x = "time point",
       y = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_blank(),
    panel.spacing = unit(1, "line")
  )

module_line_plot_460

ggsave("module_line_plot_460.png", height = 10, width = 8, bg = "white")

##################################
# Heatmap
##################################

modules_mean_z$mean.z %>% summary()
quantile(modules_mean_z$mean.z, c(0.05, 0.95))

modules_mean_z <- modules_mean_z %>% 
  mutate(mean.z.clipped = case_when(
    mean.z > 1.5 ~ 1.5,
    mean.z < -1.5 ~ -1.5,
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
  facet_grid(. ~ Treatment, scales = "free", space = "free") +
  geom_tile(aes(fill = mean.z)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1.5, 1.5),
                       breaks = c(-1.5, 0, 1.5), labels = c("< -1.5", "0", "> 1.5")) +
  labs(x = "Time point",
       y = "Module",
       fill = "z score") +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    legend.position = "top"
  )

ggsave("module_heatmap.png", height = 10, width = 8, bg = "white")

##################################
# Gene co-expression graphs
##################################

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


##################################
# Annotate modules
##################################

write.csv(x = my_network_modules, file = "network_modules.csv")

Annotations <- read_tsv("BLAST/DmProteinsTairdbBLASTconcat.txt")
colnames(Annotations) <- c(1:13)

head(Annotations)

Annotation_subset <- Annotations %>% 
  select(1:3)

head(Annotation_subset)

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

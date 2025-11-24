library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggpubr)




########################### FungiGut ################################


# Load data
Celiac_raw_FUN <- read_rds("Celiacos/phyloseq_celiacos_0.99_raw.rds")
Meta_raw_FUN <- read_rds("Metacardis/phyloseq_meta_0.99_raw.rds")
View(otu_table(Celiac_raw_FUN))
View(otu_table(Meta_raw_FUN))

sample_names(Meta_raw_FUN) <- gsub("\\.fastq$", "", sample_names(Meta_raw_FUN))

ps_fun_raw <- merge_phyloseq(Celiac_raw_FUN, Meta_raw_FUN)

# Build metadata
metadata_fun_raw <- data.frame(Sample = sample_names(ps_fun_raw)) |>
  mutate(Group = ifelse(grepl("^D", Sample), "Celiac", "Non_Celiac")) |>
  distinct() |>
  tibble::column_to_rownames("Sample")


sample_data(ps_fun_raw) <- sample_data(metadata_fun_raw)

# Cleaning

ps_fun_raw <- prune_samples(sample_sums(ps_fun_raw) > 0, ps_fun_raw)

ps_fun_raw <- prune_taxa(taxa_sums(ps_fun_raw) > 0, ps_fun_raw)


# Alpha diversity
p_alpha_fun <- plot_richness(ps_fun_raw,
                   color = "Group",
                   x = "Group",
                   measures = c("Chao1", "Simpson", "Shannon"))
p_alpha_fun$layers[[2]] <- NULL

plot_rich_fun <- p_alpha_fun + geom_boxplot(aes(fill = Group), alpha=.6, outlier.shape = NA, color = "black") + 
  scale_color_manual(values = c("#B2ABD2", "#FDE68A")) + 
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7, size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#B2ABD2", "#FDE68A")) + 
  theme_bw() +
  stat_summary(
    aes(group = Group),
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white") +
  labs(title = "Alpha Diversity (FungiGut)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 10, face = "bold"),
        axis.text.y  = element_text(size = 10, face = "bold"))
plot_rich_fun


# P values
alfas <- plot_rich_fun[["plot_env"]][["DF"]] %>% select(Chao1, Shannon, Simpson, Group)


compare_means(Chao1 ~ Group, data = alfas, method = "wilcox.test", paired = FALSE)
compare_means(Simpson ~ Group, data = alfas, method = "wilcox.test", paired = FALSE)
compare_means(Shannon ~ Group, data = alfas, method = "wilcox.test", paired = FALSE)



########################### MiCoP ###################################

# Load data
Celiac_raw_MI <- read_rds("Celiacos/MiCoP_Celiacos_0.99_raw.rds")
Meta_raw_MI <- read_rds("Metacardis/MiCoP_Metacardis_0.99_raw.rds")
View(otu_table(Celiac_raw_MI))
View(otu_table(Meta_raw_MI))


ps_mi_raw <- merge_phyloseq(Celiac_raw_MI, Meta_raw_MI)

# Build metadata 
metadata_mi_raw <- data.frame(Sample = sample_names(ps_mi_raw)) |>
  mutate(Group = ifelse(grepl("^D", Sample), "Celiac", "Non_Celiac")) |>
  distinct() |>
  tibble::column_to_rownames("Sample")


sample_data(ps_mi_raw) <- sample_data(metadata_mi_raw)

# Cleaning

ps_mi_raw <- prune_samples(sample_sums(ps_fun_raw) > 0, ps_fun_raw)

ps_mi_raw <- prune_taxa(taxa_sums(ps_fun_raw) > 0, ps_fun_raw)

# Alpha diversity
p_alpha_mi <- plot_richness(ps_mi_raw,
                             color = "Group",
                             x = "Group",
                             measures = c("Chao1", "Simpson", "Shannon"))
p_alpha_mi$layers[[2]] <- NULL

plot_rich_mi <- p_alpha_mi + geom_boxplot(aes(fill = Group), alpha=.6, outlier.shape = NA, color = "black") + 
  scale_color_manual(values = c("#B2ABD2", "#FDE68A")) + 
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7, size = 1.5, show.legend = FALSE) +
  scale_fill_manual(values = c("#B2ABD2", "#FDE68A")) + 
  theme_bw() +
  stat_summary(
    aes(group = Group),
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    color = "black",
    fill = "white") +
  labs(title = "Alpha Diversity (MiCoP)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 10, face = "bold"),
        axis.text.y  = element_text(size = 10, face = "bold"))
plot_rich_mi


# P values
alfas_mi <- plot_rich_fun[["plot_env"]][["DF"]] %>% select(Chao1, Shannon, Simpson, Group)


compare_means(Chao1 ~ Group, data = alfas_mi, method = "wilcox.test", paired = FALSE)
compare_means(Simpson ~ Group, data = alfas_mi, method = "wilcox.test", paired = FALSE)
compare_means(Shannon ~ Group, data = alfas_mi, method = "wilcox.test", paired = FALSE)

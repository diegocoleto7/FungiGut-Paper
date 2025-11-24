library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(rstatix)
library(ggpubr)


thr <- 0.001

map_taxa <- c(
  "Candida_auris" = "Candidozyma_auris",
  "Candida_duobushaemulonis" = "Candidozyma_duobushaemuli",
  "Candida_pseudohaemulonii" = "Candidozyma_pseudohaemuli")

bacteria <- c(
  "Bacteroides_fragilis",
  "Faecalibacterium_prausnitzii",
  "Escherichia_coli",
  "Ruminococcus_bromii",
  "Akkermansia_muciniphila")


# Metadata Set1 

metadata1 <- read_tsv("newmock/set1_microbiome.tsv")

metadata1 <- metadata1 %>%
  filter(!taxon %in% bacteria) %>%
  group_by(sample) %>%
  mutate(tax_abundance = tax_abundance / sum(tax_abundance)) %>%
  ungroup()

metadata1$tax_abundance <- metadata1$tax_abundance * 100

metadata1 <- metadata1 %>%
  mutate(taxon = recode(as.character(taxon), !!!map_taxa)) %>%
  mutate(Genus = sub("_.*", "", taxon)) %>%
  rename(Sample = sample, OTU = taxon)

metadata1_bacclean <- metadata1 %>%
  mutate(Sample = paste0(Sample, "_bacclean"))

metadata_set1_all <- bind_rows(metadata1, metadata1_bacclean)

# Set1 y Set1_filtered

ps1  <- read_rds("newmock/phyloseq_object_set1_newmock.rds")
ps1f <- read_rds("newmock/phyloseq_object_set1_f_newmock.rds")

ps_set1_all <- merge_phyloseq(ps1, ps1f)

df_ps_set1_all <- psmelt(ps_set1_all) %>%
  mutate(
    Species = sub(" ", "_", Species),
    OTU = sub(" ", "_", OTU),
    Abundance = as.numeric(Abundance)) %>%
  filter(Abundance >= thr) %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()


df_comp1 <- metadata_set1_all %>%
  full_join(df_ps_set1_all, by = c("Sample", "OTU")) %>%
  replace_na(list(tax_abundance = 0, Abundance = 0)) %>%
  mutate(
    diff  = tax_abundance - Abundance,
    abs_error  = abs(diff),
    truth  = tax_abundance > 0,
    pred = Abundance > 0)

pr_sample1 <- df_comp1 %>%
  group_by(Sample) %>%
  summarise(
    precision = ifelse(sum(pred)  == 0, NA_real_, sum(truth & pred) / sum(pred)),
    recall = ifelse(sum(truth) == 0, NA_real_, sum(truth & pred) / sum(truth)),
    mse = mean(diff^2),
    rmse = sqrt(mse)) %>%
  mutate(
    F1 = ifelse(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)),
    Group = if_else(grepl("_bacclean$", Sample), "Set1_filtered", "Set1")) %>%
  ungroup()



# Metadata Set2

metadata2 <- read_tsv("newmock/set2_mycobiome.tsv")

metadata2$tax_abundance <- metadata2$tax_abundance * 100

metadata2 <- metadata2 %>%
  mutate(
    num  = readr::parse_number(sample),
    pref = str_remove(sample, "\\d+$"),
    sample = if_else(num >= 1 & num <= 20,
                     paste0(pref, num + 20L),
                     sample)) %>%
  select(-num, -pref)

metadata2 <- metadata2 %>%
  mutate(taxon = recode(as.character(taxon), !!!map_taxa),
         Genus = sub("_.*", "", taxon)) %>%
  rename(Sample = sample, OTU = taxon)

# Set2

ps2 <- read_rds("newmock/phyloseq_object_set2_newmock.rds")

old_2 <- sample_names(ps2)
m <- stringr::str_match(old_2, "^(.*?)(\\d+)$")
pref <- m[, 2]
num  <- suppressWarnings(as.integer(m[, 3]))
to_shift <- !is.na(num) & num >= 1 & num <= 20
new_2 <- old_2
new_2[to_shift] <- paste0(pref[to_shift], num[to_shift] + 20L)
sample_names(ps2) <- new_2

df_ps_set2 <- psmelt(ps2) %>%
  mutate(
    Species = sub(" ", "_", Species),
    OTU  = sub(" ", "_", OTU),
    Abundance = as.numeric(Abundance)) %>%
  filter(Abundance >= thr) %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()



df_comp2 <- metadata2 %>%
  full_join(df_ps_set2, by = c("Sample", "OTU")) %>%
  replace_na(list(tax_abundance = 0, Abundance = 0)) %>%
  mutate(
    diff = tax_abundance - Abundance,
    abs_error  = abs(diff),
    truth  = tax_abundance > 0,
    pred = Abundance > 0)

pr_sample2 <- df_comp2 %>%
  group_by(Sample) %>%
  summarise(
    precision = ifelse(sum(pred)  == 0, NA_real_, sum(truth & pred) / sum(pred)),
    recall = ifelse(sum(truth) == 0, NA_real_, sum(truth & pred) / sum(truth)),
    mse = mean(diff^2),
    rmse = sqrt(mse)) %>%
  mutate(
    F1 = ifelse(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)),
    Group = "Set2") %>%
  ungroup()


# Join metrics

metrics_combined <- bind_rows(
  pr_sample1 %>% select(Sample, Group, F1, rmse),
  pr_sample2 %>% select(Sample, Group, F1, rmse))

metrics_long <- metrics_combined %>%
  pivot_longer(
    cols = c(F1, rmse),
    names_to  = "Metric",
    values_to = "Value") %>%
  mutate(
    Group  = factor(Group, levels = c("Set1", "Set1_filtered", "Set2")),
    Metric = factor(Metric, levels = c("F1", "rmse")))


my_comparisons <- list(
  c("Set1", "Set1_filtered"),
  c("Set2", "Set1_filtered"))

group_colors <- c(
  "Set1"          = "#E39A55",
  "Set1_filtered" = "#3FA79B",
  "Set2"          = "#E8CC4A")


# F1

df_F1 <- metrics_long %>%
  filter(Metric == "F1")

stat_F1 <- df_F1 %>%
  rstatix::wilcox_test(Value ~ Group, comparisons = my_comparisons) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "Group", comparisons = my_comparisons)

g_F1 <- ggplot(df_F1, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, alpha = 0.7, size = 1.8, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  ggpubr::stat_pvalue_manual(
    stat_F1,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(title = "F1-score", x = NULL, y = "Value") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x  = element_text(size = 10, face = "bold")
  )

### RMSE


df_RMSE <- metrics_long %>%
  filter(Metric == "rmse")

stat_RMSE <- df_RMSE %>%
  rstatix::wilcox_test(Value ~ Group, comparisons = my_comparisons) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  rstatix::add_significance("p.adj") %>%
  rstatix::add_xy_position(x = "Group", comparisons = my_comparisons)

g_RMSE <- ggplot(df_RMSE, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, alpha = 0.7, size = 1.8, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  ggpubr::stat_pvalue_manual(
    stat_RMSE,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE) +
  labs(title = "RMSE", x = NULL, y = "Value") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13),
    axis.text.x  = element_text(size = 10, face = "bold"))


# Joined plot

g_F1_l   <- g_F1 + labs(fill = "Mock Set", color = "Mock Set") + theme(legend.position = "bottom")
g_RMSE_l <- g_RMSE + labs(fill = "Mock Set", color = "Mock Set") + theme(legend.position = "bottom")

final_plot <- ggpubr::ggarrange(
  g_F1_l, g_RMSE_l,
  ncol = 2,
  align = "hv",
  common.legend = TRUE,
  legend = "bottom",
  labels = c("A", "B"))

ggsave(
  "Graph/RMSE-F1.png",
  plot = final_plot,
  dpi = 800,
  width  = 14,
  height = 16,
  units  = "cm")

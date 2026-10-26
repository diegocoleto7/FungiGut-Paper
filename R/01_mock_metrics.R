# Mock communities: F1-score and RMSE of FungiGut against the expected composition.
# Output: Figure 4, per-sample metrics, summary metrics and false positives.

source("R/00_setup.R")
library(stringr)
library(rstatix)

thr <- 0.001

map_taxa <- c(
  "Candida_auris"            = "Candidozyma_auris",
  "Candida_duobushaemulonis" = "Candidozyma_duobushaemuli",
  "Candida_pseudohaemulonii" = "Candidozyma_pseudohaemuli")

bacteria <- c(
  "Bacteroides_fragilis",
  "Faecalibacterium_prausnitzii",
  "Escherichia_coli",
  "Ruminococcus_bromii",
  "Akkermansia_muciniphila")

# Per-sample precision, recall, F1 and RMSE
sample_metrics <- function(expected, observed) {
  expected %>%
    full_join(observed, by = c("Sample", "OTU")) %>%
    replace_na(list(tax_abundance = 0, Abundance = 0)) %>%
    mutate(
      diff  = tax_abundance - Abundance,
      truth = tax_abundance > 0,
      pred  = Abundance > 0) %>%
    group_by(Sample) %>%
    summarise(
      precision = ifelse(sum(pred)  == 0, NA_real_, sum(truth & pred) / sum(pred)),
      recall    = ifelse(sum(truth) == 0, NA_real_, sum(truth & pred) / sum(truth)),
      rmse      = sqrt(mean(diff^2)),
      .groups = "drop") %>%
    mutate(F1 = ifelse(
      is.na(precision) | is.na(recall) | (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)))
}

# Shift sample1-20 to sample21-40 so Set2 names do not clash with Set1
shift_sample <- function(x) {
  m <- str_match(x, "^(.*?)(\\d+)$")
  num <- suppressWarnings(as.integer(m[, 3]))
  to_shift <- !is.na(num) & num >= 1 & num <= 20
  x[to_shift] <- paste0(m[to_shift, 2], num[to_shift] + 20L)
  x
}


########################### Set1 / Set1_filtered ###########################

# Expected: fungi only, rescaled to 100 %
expected1 <- read_tsv("data/mock_sets/set1_microbiome.tsv", show_col_types = FALSE) %>%
  filter(!taxon %in% bacteria) %>%
  group_by(sample) %>%
  mutate(tax_abundance = 100 * tax_abundance / sum(tax_abundance)) %>%
  ungroup() %>%
  mutate(taxon = recode(as.character(taxon), !!!map_taxa)) %>%
  rename(Sample = sample, OTU = taxon)

expected1_all <- bind_rows(
  expected1,
  expected1 %>% mutate(Sample = paste0(Sample, "_bacclean")))

ps1  <- read_rds("data/mock_sets/phyloseq_object_set1_newmock.rds")
ps1f <- read_rds("data/mock_sets/phyloseq_object_set1_f_newmock.rds")
observed1 <- rel_abundance(merge_phyloseq(ps1, ps1f), thr)

metrics1 <- sample_metrics(expected1_all, observed1) %>%
  mutate(Group = if_else(grepl("_bacclean$", Sample), "Set1_filtered", "Set1"))


################################## Set2 ####################################

expected2 <- read_tsv("data/mock_sets/set2_mycobiome.tsv", show_col_types = FALSE) %>%
  mutate(
    tax_abundance = 100 * tax_abundance,
    sample = shift_sample(sample),
    taxon  = recode(as.character(taxon), !!!map_taxa)) %>%
  rename(Sample = sample, OTU = taxon)

ps2 <- read_rds("data/mock_sets/phyloseq_object_set2_newmock.rds")
sample_names(ps2) <- shift_sample(sample_names(ps2))
observed2 <- rel_abundance(ps2, thr)

metrics2 <- sample_metrics(expected2, observed2) %>%
  mutate(Group = "Set2")


################################# Metrics ##################################

metrics <- bind_rows(metrics1, metrics2) %>%
  mutate(Group = factor(Group, levels = c("Set1", "Set1_filtered", "Set2"))) %>%
  select(Sample, Group, F1, rmse, precision, recall)

write_csv(metrics, "results/tables/mock_metrics_per_sample.csv")

metrics_summary <- metrics %>%
  group_by(Group) %>%
  summarise(
    F1        = mean(F1, na.rm = TRUE),
    RMSE      = mean(rmse, na.rm = TRUE),
    Precision = mean(precision, na.rm = TRUE),
    Recall    = mean(recall, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(across(c(F1, RMSE, Precision, Recall), ~ round(.x, 3)))

write_csv(metrics_summary, "results/tables/mock_metrics_summary.csv")
print(metrics_summary)


# False positives in Set1_filtered and Set2
expected_species <- unique(c(expected1$OTU, expected2$OTU))

false_positives <- bind_rows(
  observed1 %>% filter(grepl("_bacclean$", Sample)),
  observed2) %>%
  filter(!Species %in% expected_species) %>%
  group_by(Species) %>%
  summarise(
    mean_abund   = mean(Abundance, na.rm = TRUE),
    median_abund = median(Abundance, na.rm = TRUE),
    n_samples    = n_distinct(Sample),
    .groups = "drop") %>%
  arrange(desc(mean_abund))

write_csv(false_positives, "results/tables/mock_false_positives.csv")


################################# Figure 4 #################################

my_comparisons <- list(
  c("Set1", "Set1_filtered"),
  c("Set2", "Set1_filtered"))

metrics_long <- metrics %>%
  pivot_longer(cols = c(F1, rmse), names_to = "Metric", values_to = "Value")

plot_metric <- function(metric, title) {
  df <- metrics_long %>% filter(Metric == metric)

  stat <- df %>%
    wilcox_test(Value ~ Group, comparisons = my_comparisons) %>%
    adjust_pvalue(method = "holm") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Group", comparisons = my_comparisons)

  write_csv(
    stat %>% select(group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif),
    sprintf("results/tables/mock_wilcoxon_%s.csv", metric))

  ggplot(df, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = Group), width = 0.1, alpha = 0.7, size = 1.8, show.legend = FALSE) +
    scale_fill_manual(values = mock_colors) +
    scale_color_manual(values = mock_colors) +
    stat_pvalue_manual(stat, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
    labs(title = title, x = NULL, y = "Value", fill = "Mock Set", color = "Mock Set") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.y = element_text(face = "bold", size = 13),
      axis.text.x  = element_text(size = 10, face = "bold"))
}

set.seed(1)
figure4 <- ggarrange(
  plot_metric("F1", "F1-score"),
  plot_metric("rmse", "RMSE"),
  ncol = 2,
  align = "hv",
  common.legend = TRUE,
  legend = "bottom",
  labels = c("A", "B"))

ggsave(
  "results/figures/Figure4_mock_metrics.png",
  plot = figure4,
  width = 14, height = 16, units = "cm",
  dpi = 800,
  bg = "white")

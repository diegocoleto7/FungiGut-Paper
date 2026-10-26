# Alpha and beta diversity of the NRCD and MetaCardis cohorts, FungiGut vs MiCoP.
# Alpha diversity uses raw read counts; beta diversity uses the batch-corrected
# objects from 03_batch_correction.R.
# Output: Figure 5, alpha-diversity Wilcoxon tests and beta-diversity PERMANOVA.

source("R/00_setup.R")
library(vegan)


############################## Alpha diversity ##############################

alpha_measures <- c("Chao1", "Simpson", "Shannon")

load_raw <- function(nrcd_file, metacardis_file) {
  ps_nrcd <- read_rds(nrcd_file)
  ps_meta <- read_rds(metacardis_file)
  sample_names(ps_meta) <- gsub("\\.fastq$", "", sample_names(ps_meta))

  ps <- merge_phyloseq(ps_nrcd, ps_meta)
  sample_data(ps) <- sample_data(data.frame(Group = cohort_group(sample_names(ps)), row.names = sample_names(ps)))
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

plot_alpha <- function(ps, title) {
  p <- plot_richness(ps, color = "Group", x = "Group", measures = alpha_measures)
  p$layers[[2]] <- NULL

  p +
    geom_boxplot(aes(fill = Group), alpha = .6, outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Group), width = 0.2, alpha = 0.7, size = 1.5, show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    stat_summary(aes(group = Group), fun = mean, geom = "point", shape = 21, size = 3, color = "black", fill = "white") +
    labs(title = title) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 10, face = "bold"),
      axis.text.y  = element_text(size = 10, face = "bold"))
}

alpha_tests <- function(ps, tool) {
  alpha <- estimate_richness(ps, measures = alpha_measures)
  alpha$Group <- sample_data(ps)[rownames(alpha), ]$Group

  lapply(alpha_measures, function(m) {
    summ <- alpha %>%
      group_by(Group) %>%
      summarise(n = n(), mean = mean(.data[[m]]), median = median(.data[[m]]), .groups = "drop")
    test <- wilcox.test(alpha[[m]] ~ alpha$Group, exact = FALSE)
    tibble(
      Tool = tool, Measure = m,
      n_NRCD = summ$n[summ$Group == "NRCD"], n_Control = summ$n[summ$Group == "Control"],
      mean_NRCD = summ$mean[summ$Group == "NRCD"], mean_Control = summ$mean[summ$Group == "Control"],
      median_NRCD = summ$median[summ$Group == "NRCD"], median_Control = summ$median[summ$Group == "Control"],
      W = unname(test$statistic), p = test$p.value)
  }) %>%
    bind_rows()
}

ps_fun_raw <- load_raw("data/nrcd/phyloseq_celiacos_0.99_raw.rds", "data/metacardis/phyloseq_meta_0.99_raw.rds")
ps_mi_raw  <- load_raw("data/nrcd/MiCoP_Celiacos_0.99_raw.rds",   "data/metacardis/MiCoP_Metacardis_0.99_raw.rds")

alpha_results <- bind_rows(
  alpha_tests(ps_fun_raw, "FungiGut"),
  alpha_tests(ps_mi_raw,  "MiCoP"))

write_csv(alpha_results, "results/tables/alpha_diversity_wilcoxon.csv")
print(alpha_results)

set.seed(1)
plot_rich_fun <- plot_alpha(ps_fun_raw, "Alpha Diversity (FungiGut)")
plot_rich_mi  <- plot_alpha(ps_mi_raw,  "Alpha Diversity (MiCoP)")


############################## Beta diversity ###############################

ps_FUN <- read_rds("data/processed/batch_corrected_fungigut.rds")
ps_MI  <- read_rds("data/processed/batch_corrected_micop.rds")

beta_permanova <- function(ps, tool) {
  meta <- as(sample_data(ps), "data.frame")
  lapply(c(Bray_Curtis = "bray", Jaccard = "jaccard"), function(method) {
    set.seed(1)
    res <- adonis2(phyloseq::distance(ps, method = method) ~ Group, data = meta, permutations = 999)
    data.frame(R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1])
  }) %>%
    bind_rows(.id = "Distance") %>%
    mutate(Tool = tool, .before = 1)
}

beta_results <- bind_rows(
  beta_permanova(ps_FUN, "FungiGut"),
  beta_permanova(ps_MI,  "MiCoP"))

write_csv(beta_results, "results/tables/beta_diversity_permanova.csv")
print(beta_results)

plot_pcoa <- function(ps, title) {
  ord <- ordinate(ps, method = "MDS", distance = "bray")
  plotdf <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
  colnames(plotdf) <- c("Axis1", "Axis2")
  plotdf$Group <- sample_data(ps)[rownames(plotdf), ]$Group
  expl <- ord$values$Relative_eig

  ggplot(plotdf, aes(Axis1, Axis2)) +
    geom_point(aes(fill = Group), shape = 21, size = 3.5, alpha = 0.9, color = "grey20", stroke = 0.2) +
    scale_fill_manual(values = group_colors, name = NULL) +
    labs(
      title = title,
      x = paste0("PCoA1 (", round(expl[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(expl[2] * 100, 1), "%)")) +
    theme_bw() +
    theme(
      legend.position = "top",
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13))
}


################################# Figure 5 ##################################

figure5 <- ggarrange(
  plot_pcoa(ps_FUN, "PCoA Bray-Curtis (FungiGut)"), plot_rich_fun,
  plot_pcoa(ps_MI,  "PCoA Bray-Curtis (MiCoP)"),    plot_rich_mi,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "top",
  align = "hv")

ggsave(
  "results/figures/Figure5_diversity.png",
  plot = figure5,
  width = 10, height = 8.5, units = "in",
  dpi = 400,
  bg = "white")

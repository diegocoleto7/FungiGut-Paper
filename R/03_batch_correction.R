# Batch correction (MMUPHin) of the NRCD and MetaCardis cohorts, for FungiGut and MiCoP.
# Output: data/processed/batch_corrected_*.rds, PERMANOVA before/after and
# the supplementary batch-correction figure.

source("R/00_setup.R")
library(vegan)
library(MMUPHin)

# Merge both cohorts, add Group, drop empty samples/taxa and convert to relative abundance
load_cohorts <- function(nrcd_file, metacardis_file, drop_samples = character()) {
  ps <- merge_phyloseq(read_rds(nrcd_file), read_rds(metacardis_file))
  ps <- prune_samples(!sample_names(ps) %in% drop_samples, ps)
  sample_data(ps) <- sample_data(data.frame(Group = cohort_group(sample_names(ps)), row.names = sample_names(ps)))
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  microbiome::transform(ps, "compositional")
}

batch_correct <- function(ps, tool) {
  mat  <- as(otu_table(ps), "matrix")
  meta <- as(sample_data(ps), "data.frame")
  stopifnot(identical(colnames(mat), rownames(meta)))

  adj <- adjust_batch(
    feature_abd = mat,
    batch = "Group",
    data = meta,
    control = list(diagnostic_plot = NULL, verbose = FALSE))

  corrected <- phyloseq(
    tax_table(ps),
    otu_table(adj$feature_abd_adj, taxa_are_rows = TRUE),
    sample_data(ps))

  # PERMANOVA (Bray-Curtis) before and after correction
  permanova <- lapply(list(Unadjusted = mat, Adjusted = adj$feature_abd_adj), function(m) {
    set.seed(1)
    res <- adonis2(vegdist(t(m), method = "bray") ~ Group, data = meta, permutations = 999)
    data.frame(R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1])
  }) %>%
    bind_rows(.id = "Data") %>%
    mutate(Tool = tool, .before = 1)

  list(before = ps, after = corrected, permanova = permanova)
}

plot_pcoa_ellipse <- function(ps, title) {
  ord <- ordinate(ps, method = "MDS", distance = "bray")
  plotdf <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
  colnames(plotdf) <- c("Axis1", "Axis2")
  plotdf$Group <- sample_data(ps)[rownames(plotdf), ]$Group
  expl <- ord$values$Relative_eig

  ggplot(plotdf, aes(Axis1, Axis2)) +
    stat_ellipse(aes(color = Group), type = "norm", level = 0.8, linewidth = 0.7, show.legend = FALSE) +
    geom_point(aes(fill = Group), shape = 21, size = 3.5, alpha = 0.9, color = "grey", stroke = 0.05) +
    labs(
      title = title,
      x = paste0("PCoA1 (", round(expl[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(expl[2] * 100, 1), "%)")) +
    scale_fill_manual(values = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, guide = "none") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "top")
}


# FungiGut
fun <- batch_correct(
  load_cohorts("data/nrcd/phyloseq_object_celiacos_0.99.rds", "data/metacardis/phyloseq_object_meta_0.99.rds"),
  "FungiGut")

# MiCoP ("r1" is not a MetaCardis sample)
mi <- batch_correct(
  load_cohorts("data/nrcd/MiCoP_Celiacos_0.99.rds", "data/metacardis/MiCoP_Metacardis_0.99.rds", drop_samples = "r1"),
  "MiCoP")

saveRDS(fun$after, "data/processed/batch_corrected_fungigut.rds")
saveRDS(mi$after,  "data/processed/batch_corrected_micop.rds")

permanova <- bind_rows(fun$permanova, mi$permanova)
write_csv(permanova, "results/tables/batch_correction_permanova.csv")
print(permanova)


# Supplementary figure
panel_row <- function(res, tool) {
  ggarrange(
    plot_pcoa_ellipse(res$before, sprintf("PCoA Unadjusted (%s)", tool)),
    plot_pcoa_ellipse(res$after,  sprintf("PCoA Batch-adjusted (%s)", tool)),
    ncol = 2, nrow = 1,
    common.legend = TRUE,
    legend = "top")
}

figure_batch <- ggarrange(
  panel_row(fun, "FungiGut"),
  panel_row(mi, "MiCoP"),
  ncol = 1, nrow = 2,
  labels = c("A", "B"))

ggsave(
  "results/figures/FigureS_batch_correction.png",
  plot = figure_batch,
  width = 14, height = 10, units = "in",
  dpi = 600,
  bg = "white")

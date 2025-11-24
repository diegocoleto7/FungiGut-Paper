library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(ape) 
library(tibble)
library(vegan)
library(MMUPHin)
library(plyr)


# Load data
ps_FUN <- read_rds("Batch/phyloseq_object_corrected.rds")
ps_MI <- read_rds("Batch/MiCoP_corrected.rds")

# Unifrac
unifrac_FUN <- phyloseq::distance(ps_FUN, method = "unifrac", weighted = TRUE)
unifrac_MI <- phyloseq::distance(ps_MI, method = "unifrac", weighted = TRUE)

# PERMANOVA
permanova_FUN_u <- adonis2(unifrac_FUN ~ Group, data = as(sample_data(ps_FUN), "data.frame"))
permanova_MI_u <- adonis2(unifrac_MI ~ Group, data = as(sample_data(ps_MI), "data.frame"))

# Print results
print("PERMANOVA UniFrac FungiGut results")
print(permanova_FUN_u)
print("PERMANOVA UniFrac MiCoP results")
print(permanova_MI_u)

# Jaccard
jaccard_FUN <- phyloseq::distance(ps_FUN, method = "jaccard")
jaccard_MI <- phyloseq::distance(ps_MI, method = "jaccard")

# PERMANOVA
permanova_FUN_j <- adonis2(jaccard_FUN ~ Group, data = as(sample_data(ps_FUN), "data.frame"))
permanova_MI_j <- adonis2(jaccard_MI ~ Group, data = as(sample_data(ps_MI), "data.frame"))

# Print results
print("PERMANOVA UniFrac original results")
print(permanova_FUN_j)
print("PERMANOVA UniFrac corrected results")
print(permanova_MI_j)

# Bray Curtis
bray_FUN <- phyloseq::distance(ps_FUN, method = "bray", weighted = TRUE)
bray_MI <- phyloseq::distance(ps_MI, method = "bray", weighted = TRUE)

# PERMANOVA
permanova_FUN_b <- adonis2(bray_FUN ~ Group, data = as(sample_data(ps_FUN), "data.frame"))
permanova_MI_b <- adonis2(bray_MI ~ Group, data = as(sample_data(ps_MI), "data.frame"))

# Print results
print("PERMANOVA UniFrac original results")
print(permanova_FUN_b)
print("PERMANOVA UniFrac corrected results")
print(permanova_MI_b)


# Plot Beta diversity
ord_w_FUN <- ordinate(ps_FUN, method="MDS", distance="unifrac", weighted=TRUE)
ord_j_FUN <- ordinate(ps_FUN, method="MDS", distance="jaccard")
ord_b_FUN <- ordinate(ps_FUN, method="MDS", distance="bray")

dist_methods_FUN <- list(
  Unifrac     = ord_w_FUN,
  Jaccard     = ord_j_FUN,
  Bray_Curtis = ord_b_FUN)

meta_FUN <- as(sample_data(ps_FUN), "data.frame") |>
  tibble::rownames_to_column("Sample")

df_beta_FUN <- lapply(names(dist_methods_FUN), function(nm) {
  ord_obj <- dist_methods_FUN[[nm]]
  plotdf <- as.data.frame(ord_obj$vectors[, 1:2, drop = FALSE])
  colnames(plotdf) <- c("Axis1", "Axis2")
  plotdf$Sample <- rownames(plotdf)
  expl <- ord_obj$values$Relative_eig
  lab_facet <- sprintf(
    "%s\nPCoA1: %.1f%%, PCoA2: %.1f%%",
    nm, expl[1] * 100, expl[2] * 100)
  plotdf <- dplyr::left_join(plotdf, meta_FUN, by = "Sample")
  plotdf$distance <- lab_facet
  plotdf
}) |>
  dplyr::bind_rows()

p_beta_FUN <- ggplot(df_beta_FUN, aes(Axis1, Axis2, fill = Group)) +
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "grey20") +
  facet_wrap(~ distance, scales = "free") +
  scale_fill_manual(values = cols, name = NULL) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.text      = element_text(face = "bold", size = 10)) +
  labs(
    title = "Beta Diversity (FungiGut)",
    x = "PCoA1",
    y = "PCoA2")

p_beta_FUN



# Plot Beta diversity
ord_w_MI <- ordinate(ps_MI, method="MDS", distance="unifrac", weighted=TRUE)
ord_j_MI <- ordinate(ps_MI, method="MDS", distance="jaccard")
ord_b_MI <- ordinate(ps_MI, method="MDS", distance="bray")

dist_methods_MI <- list(
  Unifrac     = ord_w_MI,
  Jaccard     = ord_j_MI,
  Bray_Curtis = ord_b_MI)

meta_MI <- as(sample_data(ps_MI), "data.frame") |>
  tibble::rownames_to_column("Sample")

df_beta_MI <- lapply(names(dist_methods_MI), function(nm) {
  ord_obj <- dist_methods_MI[[nm]]
  plotdf <- as.data.frame(ord_obj$vectors[, 1:2, drop = FALSE])
  colnames(plotdf) <- c("Axis1", "Axis2")
  plotdf$Sample <- rownames(plotdf)
  expl <- ord_obj$values$Relative_eig
  lab_facet <- sprintf(
    "%s\nPCoA1: %.1f%%, PCoA2: %.1f%%",
    nm, expl[1] * 100, expl[2] * 100)
  plotdf <- dplyr::left_join(plotdf, meta_MI, by = "Sample")
  plotdf$distance <- lab_facet
  plotdf
}) |>
  dplyr::bind_rows()

p_beta_MI <- ggplot(df_beta_MI, aes(Axis1, Axis2, fill = Group)) +
  geom_point(shape = 21, size = 3, alpha = 0.8, color = "grey20") +
  facet_wrap(~ distance, scales = "free") +
  scale_fill_manual(values = cols, name = NULL) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.text      = element_text(face = "bold")) +
  labs(
    title = "Beta Diversity (MiCoP)",
    x = "PCoA1",
    y = "PCoA2")

p_beta_MI



## FungiGut – Unifrac

plotdf_uni_FUN <- as.data.frame(ord_w_FUN$vectors[, 1:2, drop = FALSE])
colnames(plotdf_uni_FUN) <- c("Axis1", "Axis2")
plotdf_uni_FUN$Sample <- rownames(plotdf_uni_FUN)

plotdf_uni_FUN <- dplyr::left_join(
  plotdf_uni_FUN,
  meta_FUN,
  by = "Sample")

expl_FUN <- ord_w_FUN$values$Relative_eig

p_unifrac_FUN <- ggplot(plotdf_uni_FUN, aes(Axis1, Axis2)) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey20",
    stroke = 0.2) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  labs(
    title = "PCoA Unifrac (FungiGut)",
    x = paste0("PCoA1 (", round(expl_FUN[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_FUN[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x    = element_text(face = "bold", size = 13),
    axis.title.y    = element_text(face = "bold", size = 13))

p_unifrac_FUN


## MiCoP – Unifrac

plotdf_uni_MI <- as.data.frame(ord_w_MI$vectors[, 1:2, drop = FALSE])
colnames(plotdf_uni_MI) <- c("Axis1", "Axis2")
plotdf_uni_MI$Sample <- rownames(plotdf_uni_MI)

plotdf_uni_MI <- dplyr::left_join(
  plotdf_uni_MI,
  meta_MI,
  by = "Sample")

expl_MI <- ord_w_MI$values$Relative_eig

p_unifrac_MI <- ggplot(plotdf_uni_MI, aes(Axis1, Axis2)) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey20",
    stroke = 0.2) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  labs(
    title = "PCoA Unifrac (MiCoP)",
    x = paste0("PCoA1 (", round(expl_MI[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_MI[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title.x    = element_text(face = "bold", size = 13),
    axis.title.y    = element_text(face = "bold", size = 13))

p_unifrac_MI



combo <- ggpubr::ggarrange(
  p_unifrac_FUN, plot_rich_fun,
  p_unifrac_MI,   plot_rich_mi,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = TRUE, legend = "top",
  align = "hv")
combo
 
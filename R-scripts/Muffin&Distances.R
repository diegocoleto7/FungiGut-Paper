library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(ape) 
library(vegan)
library(MMUPHin)


########################### FungiGut ################################

# Load data
FUN_celiac <- read_rds("Celiacos/phyloseq_object_celiacos_0.99.rds")
FUN_meta <- read_rds("Metacardis/phyloseq_object_meta_0.99.rds")

ps_FUN <- merge_phyloseq(FUN_celiac, FUN_meta)

# Build metadata 
metadata_FUN <- data.frame(Sample = sample_names(ps_FUN)) |>
  mutate(Group = ifelse(grepl("^D", Sample), "Celiac", "Non_Celiac")) |>
  distinct() |>
  tibble::column_to_rownames("Sample")


sample_data(ps_FUN) <- sample_data(metadata_FUN)

# Cleaning
ps_FUN <- prune_samples(sample_sums(ps_FUN) > 0, ps_FUN)
ps_FUN <- prune_taxa(taxa_sums(ps_FUN) > 0, ps_FUN)

ps_FUN <- microbiome::transform(ps_FUN,"compositional")

# As matrix
mat_FUN <- as(otu_table(ps_FUN), "matrix")
if (taxa_are_rows(ps_FUN)) {n_otus <- colSums(mat_FUN > 0)} else {n_otus <- rowSums(mat_FUN > 0)}

# As data.frame
meta_FUN <- as(sample_data(ps_FUN), "data.frame")
stopifnot(identical(colnames(mat_FUN), rownames(meta_FUN)))

#Adjusted batch
adj_FUN <- adjust_batch(
  feature_abd = mat_FUN,
  batch = "Group",
  data = meta_FUN
)

df_adj_FUN <- adj_FUN$feature_abd_adj

# Dist
dist_before_FUN <- vegdist(t(mat_FUN),method = "bray")
dist_after_FUN  <- vegdist(t(df_adj_FUN),method = "bray")

# Correction comparison
set.seed(1)
cat("PERMANOVA Before:\n")
print(adonis2(as.matrix(dist_before_FUN) ~ Group, data = meta_FUN, permutations = 999))
cat("\nPERMANOVA After:\n")
print(adonis2(as.matrix(dist_after_FUN)  ~ Group, data = meta_FUN, permutations = 999))

# Corrected objetc
corrected_FUN <- phyloseq(tax_table(ps_FUN),otu_table(df_adj_FUN,taxa_are_rows= TRUE), sample_data(ps_FUN))


tree_corrected <- rtree(ntaxa(corrected_FUN), rooted = TRUE, tip.label = taxa_names(corrected_FUN))
phy_tree(corrected_FUN) <- tree_corrected


#Plot differences
ord_before_FUN <- ordinate(ps_FUN, method="MDS", distance="bray")
ord_after_FUN <- ordinate(corrected_FUN, method="MDS", distance="bray")


#Plot before
plotdf_bf_FUN <- as.data.frame(ord_before_FUN$vectors[, 1:2, drop = FALSE])
colnames(plotdf_bf_FUN) <- c("Axis1","Axis2")
plotdf_bf_FUN$Sample <- rownames(plotdf_bf_FUN)
plotdf_bf_FUN <- cbind(plotdf_bf_FUN, meta_FUN[plotdf_bf_FUN$Sample, , drop = FALSE])

expl_bf_FUN <- ord_before_FUN$values$Relative_eig

cols <- c("Celiac" = "#B2ABD2", 
          "Non_Celiac"     = "#FDE68A")



p_before_FUN <- ggplot(plotdf_bf_FUN, aes(Axis1, Axis2)) +
  stat_ellipse(
    aes(color = Group),
    type = "norm",
    level = 0.8,
    linewidth = 0.7,
    show.legend = FALSE) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey",
    stroke = 0.05) +
  labs(
    title = "PCoA Original Data (FungiGut)",
    x = paste0("PCoA1 (", round(expl_bf_FUN[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_bf_FUN[2] * 100, 1), "%)")) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top")
p_before_FUN


#Plot corrected
plotdf_af_FUN <- as.data.frame(ord_after_FUN$vectors[, 1:2, drop = FALSE])
colnames(plotdf_af_FUN) <- c("Axis1","Axis2")
plotdf_af_FUN$Sample <- rownames(plotdf_af_FUN)
plotdf_af_FUN <- cbind(plotdf_af_FUN, meta_FUN[plotdf_af_FUN$Sample, , drop = FALSE])

expl_af_FUN <- ord_after_FUN$values$Relative_eig


p_after_FUN <- ggplot(plotdf_af_FUN, aes(Axis1, Axis2)) +
  stat_ellipse(
    aes(color = Group),
    type = "norm",
    level = 0.8,
    linewidth = 0.7,
    show.legend = FALSE) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey",
    stroke = 0.05) +
  labs(
    title = "PCoA Corrected Data (FungiGut)",
    x = paste0("PCoA1 (", round(expl_af_FUN[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_af_FUN[2] * 100, 1), "%)")) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top")
p_after_FUN




########################### MiCoP ###################################


# Load data
MI_celiac <- read_rds("Celiacos/MiCoP_Celiacos_0.99.rds")
MI_meta <- read_rds("Metacardis/MiCoP_Metacardis_0.99.rds")

ps_MI <- merge_phyloseq(MI_celiac, MI_meta)

ps_MI <- prune_samples(sample_names(ps_MI) != "r1", ps_MI)

# Build metadata 
metadata_MI <- data.frame(Sample = sample_names(ps_MI)) |>
  mutate(Group = ifelse(grepl("^D", Sample), "Celiac", "Non_Celiac")) |>
  distinct() |>
  tibble::column_to_rownames("Sample")


sample_data(ps_MI) <- sample_data(metadata_MI)

# Cleaning
ps_MI <- prune_samples(sample_sums(ps_MI) > 0, ps_MI)
ps_MI <- prune_taxa(taxa_sums(ps_MI) > 0, ps_MI)

ps_MI <- microbiome::transform(ps_MI,"compositional")

# As matrix
mat_MI <- as(otu_table(ps_MI), "matrix")
if (taxa_are_rows(ps_MI)) {n_otus <- colSums(mat_MI > 0)} else {n_otus <- rowSums(mat_MI > 0)}

# As data.frame
meta_MI <- as(sample_data(ps_MI), "data.frame")
stopifnot(identical(colnames(mat_MI), rownames(meta_MI)))

#Adjusted batch
adj_MI <- adjust_batch(
  feature_abd = mat_MI,
  batch = "Group",
  data = meta_MI)

df_adj_MI <- adj_MI$feature_abd_adj

# Dist
dist_before_MI <- vegdist(t(mat_MI),method = "bray")
dist_after_MI  <- vegdist(t(df_adj_MI),method = "bray")

# Correction comparison
set.seed(1)
cat("PERMANOVA Before:\n")
print(adonis2(as.matrix(dist_before_MI) ~ Group, data = meta_MI, permutations = 999))
cat("\nPERMANOVA After:\n")
print(adonis2(as.matrix(dist_after_MI)  ~ Group, data = meta_MI, permutations = 999))

# Corrected objetc
corrected_MI <- phyloseq(tax_table(ps_MI),otu_table(df_adj_MI,taxa_are_rows= TRUE), sample_data(ps_MI))

# Tree creation
# tree_original <- rtree(ntaxa(ps_clean), rooted = TRUE, tip.label = taxa_names(ps_clean))
# phy_tree(ps_clean) <- tree_original

tree_corrected_MI <- rtree(ntaxa(corrected_MI), rooted = TRUE, tip.label = taxa_names(corrected_MI))
phy_tree(corrected_MI) <- tree_corrected_MI


#Plot differences
ord_before_MI <- ordinate(ps_MI, method="MDS", distance="bray")
ord_after_MI <- ordinate(corrected_MI, method="MDS", distance="bray")


#Plot before
plotdf_bf_MI <- as.data.frame(ord_before_MI$vectors[, 1:2, drop = FALSE])
colnames(plotdf_bf_MI) <- c("Axis1","Axis2")
plotdf_bf_MI$Sample <- rownames(plotdf_bf_MI)
plotdf_bf_MI <- cbind(plotdf_bf_MI, meta_MI[plotdf_bf_MI$Sample, , drop = FALSE])

expl_bf_MI <- ord_before_MI$values$Relative_eig

p_before_MI <- ggplot(plotdf_bf_MI, aes(Axis1, Axis2)) +
  stat_ellipse(
    aes(color = Group),
    type = "norm",
    level = 0.8,
    linewidth = 0.7,
    show.legend = FALSE) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey",
    stroke = 0.05) +
  labs(
    title = "PCoA Original Data (MiCoP)",
    x = paste0("PCoA1 (", round(expl_bf_MI[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_bf_MI[2] * 100, 1), "%)")) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top")
p_before_MI



#Plot corrected
plotdf_af_MI <- as.data.frame(ord_after_MI$vectors[, 1:2, drop = FALSE])
colnames(plotdf_af_MI) <- c("Axis1","Axis2")
plotdf_af_MI$Sample <- rownames(plotdf_af_MI)
plotdf_af_MI <- cbind(plotdf_af_MI, meta_MI[plotdf_af_MI$Sample, , drop = FALSE])

expl_af_MI <- ord_after_MI$values$Relative_eig


p_after_MI <- ggplot(plotdf_af_MI, aes(Axis1, Axis2)) +
  stat_ellipse(
    aes(color = Group),
    type = "norm",
    level = 0.8,
    linewidth = 0.7,
    show.legend = FALSE) +
  geom_point(
    aes(fill = Group),
    shape = 21,
    size = 3.5,
    alpha = 0.9,
    color = "grey",
    stroke = 0.05) +
  labs(
    title = "PCoA Corrected Data (MiCoP)",
    x = paste0("PCoA1 (", round(expl_af_MI[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(expl_af_MI[2] * 100, 1), "%)")) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_color_manual(values = cols, guide = "none") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top")
p_after_MI




#Figure 


saveRDS(corrected_FUN, file = "Batch/phyloseq_object_corrected.rds")
saveRDS(corrected_MI, file = "Batch/MiCoP_corrected.rds")


library(ggpubr)

# Change legend
p_before_FUN_noleg <- p_before_FUN +  labs(title = "PCoA Unadjusted (FungiGut)")
p_after_FUN_noleg  <- p_after_FUN  + labs(title = "PCoA Batch-adjusted (FungiGut)")
p_before_MI_noleg  <- p_before_MI  + labs(title = "PCoA Unadjusted (MiCoP)")
p_after_MI_noleg   <- p_after_MI   + labs(title = "PCoA Batch-adjusted (MiCoP)")

# FungiGut (before / after)
panel_FUN <- ggarrange(
  p_before_FUN_noleg, p_after_FUN_noleg,
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "top")

# MiCoP (before / after)
panel_MI <- ggarrange(
  p_before_MI_noleg, p_after_MI_noleg,
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "top")

# Figure
fig_AB <- ggarrange(
  panel_FUN, panel_MI,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  heights = c(1, 1))

fig_AB


ggsave(
  "Batch/Batch_correction.png",
  plot = fig_AB,
  dpi = 800,)

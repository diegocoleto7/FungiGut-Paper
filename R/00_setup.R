# Shared settings for the FungiGut analysis.
# Every script sources this file and must be run from the repository root.

suppressPackageStartupMessages({
  library(readr)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Colours
group_colors <- c(
  "NRCD"    = "#B2ABD2",
  "Control" = "#FDE68A")

mock_colors <- c(
  "Set1"          = "#E39A55",
  "Set1_filtered" = "#3FA79B",
  "Set2"          = "#E8CC4A")

# NRCD samples start with "D"; the rest are MetaCardis controls
cohort_group <- function(samples) {
  factor(ifelse(grepl("^D", samples), "NRCD", "Control"), levels = c("NRCD", "Control"))
}

# Melt a phyloseq object, drop rows below `thr` and rescale each sample to 100 %
rel_abundance <- function(ps, thr) {
  psmelt(ps) %>%
    mutate(
      Species   = sub(" ", "_", Species),
      OTU       = sub(" ", "_", OTU),
      Abundance = as.numeric(Abundance)) %>%
    filter(Abundance >= thr) %>%
    group_by(Sample) %>%
    mutate(Abundance = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
    ungroup()
}

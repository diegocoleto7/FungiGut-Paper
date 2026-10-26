# Mean relative abundance and prevalence of each species per cohort, FungiGut vs MiCoP.
# Uses the batch-corrected objects from 03_batch_correction.R.
# Output: Supplementary Tables 10 (MetaCardis / Control) and 11 (NRCD).

source("R/00_setup.R")

thr <- 0.0001

tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

abundance_by_group <- function(ps) {
  df <- rel_abundance(ps, thr)

  taxonomy <- df %>%
    select(Species, all_of(tax_cols)) %>%
    distinct(Species, .keep_all = TRUE)

  prevalence <- df %>%
    group_by(Group, Species) %>%
    summarise(n_samples_present = n_distinct(Sample[Abundance > 0]), .groups = "drop")

  # Absent species count as 0 when averaging over all samples of the group
  df %>%
    select(Sample, Species, Abundance) %>%
    complete(Sample, Species, fill = list(Abundance = 0)) %>%
    left_join(distinct(df, Sample, Group), by = "Sample") %>%
    group_by(Group, Species) %>%
    summarise(Abundance_mean = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    filter(Abundance_mean > 0) %>%
    left_join(taxonomy, by = "Species") %>%
    left_join(prevalence, by = c("Group", "Species")) %>%
    select(Group, Species, all_of(tax_cols), Abundance_mean, n_samples_present) %>%
    arrange(Group, desc(Abundance_mean))
}

tools <- list(
  fungigut = read_rds("data/processed/batch_corrected_fungigut.rds"),
  micop    = read_rds("data/processed/batch_corrected_micop.rds"))

for (tool in names(tools)) {
  tab <- abundance_by_group(tools[[tool]])
  for (g in unique(as.character(tab$Group))) {
    out <- tab %>% filter(Group == g) %>% select(-Group)
    write_csv(out, sprintf("results/tables/%s_%s_abundance_mean.csv", tool, tolower(g)))
    cat(sprintf("%-9s %-8s %d species\n", tool, g, nrow(out)))
  }
}

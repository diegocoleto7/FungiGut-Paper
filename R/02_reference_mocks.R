# Reference mock communities (Zymo and Pool DNA) profiled with FungiGut and MiCoP.
# Output: Supplementary Tables 7 (Zymo) and 8 (Pool DNA).

source("R/00_setup.R")

thr <- 0.1

tax_cols <- c("OTU", "Sample", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

profile_table <- function(file) {
  rel_abundance(read_rds(file), thr) %>%
    select(any_of(tax_cols)) %>%
    arrange(desc(Abundance))
}

mocks <- list(
  zymo_fungigut     = "data/zymo/phyloseq_object_Zymo_0.99.rds",
  zymo_micop        = "data/zymo/MiCoP_Zymo_0.99.rds",
  pool_dna_fungigut = "data/pool_dna/phyloseq_object_44_0.99.rds",
  pool_dna_micop    = "data/pool_dna/MiCoP_44_0.99.rds")

for (nm in names(mocks)) {
  tab <- profile_table(mocks[[nm]])
  write_csv(tab, sprintf("results/tables/%s_abundance.csv", nm))
  cat(sprintf("%-18s %d taxa\n", nm, nrow(tab)))
}

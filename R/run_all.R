# Run the whole analysis from the repository root: Rscript R/run_all.R

scripts <- c(
  "R/01_mock_metrics.R",
  "R/02_reference_mocks.R",
  "R/03_batch_correction.R",
  "R/04_diversity.R",
  "R/05_cohort_abundance.R")

for (s in scripts) {
  message("== ", s)
  source(s, local = new.env())
}

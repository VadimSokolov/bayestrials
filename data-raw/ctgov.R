# Build ctgov_counts and ctgov_bayesian_records from the cached CSVs.
# Run from the book root: Rscript bayestrials/data-raw/ctgov.R
d <- "data/ctgov"
ctgov_counts <- read.csv(file.path(d, "counts_by_year.csv"), stringsAsFactors = FALSE)
ctgov_bayesian_records <- read.csv(file.path(d, "bayesian_records.csv"), stringsAsFactors = FALSE)
save(ctgov_counts, file = "bayestrials/data/ctgov_counts.rda", compress = "xz")
save(ctgov_bayesian_records, file = "bayestrials/data/ctgov_bayesian_records.rda", compress = "xz")

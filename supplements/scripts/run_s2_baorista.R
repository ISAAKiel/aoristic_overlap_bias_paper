#!/usr/bin/env Rscript

# Run the Supplement 2 Baorista reconstruction and cache its outputs.

source("R/baorista_helpers.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
set.seed(42)

input_file <- "outputs/s1_assigned_events.csv"
output_curve_file <- "outputs/s2_baorista_curve_binned.csv"
output_posterior_file <- "outputs/s2_baorista_posterior_binned.csv"
output_diagnostics_file <- "outputs/s2_baorista_diagnostics.csv"
output_fit_file <- "outputs/s2_baorista_fit.rds"

time_min <- -4100L
time_max <- -1800L
bin_width <- 50L
niter <- 50000L
nburnin <- 25000L
thin <- 50L
nchains <- 4L
parallel <- TRUE
sigma_prior <- "dexp(rate=1)"

if (!file.exists(input_file)) stop("Render Supplement 1 first.")
intervals <- readr::read_csv(input_file, show_col_types = FALSE)

cat("Running Supplement 2 Baorista model for ", nrow(intervals), " intervals.\n", sep = "")

result <- run_baorista_intervals(intervals = intervals, time_min = time_min,
                                 time_max = time_max, bin_width = bin_width,
                                 niter = niter, nburnin = nburnin, thin = thin,
                                 nchains = nchains, sigma_prior = sigma_prior,
                                 parallel = parallel, seeds = seq_len(nchains),
                                 keep_fit = TRUE)

readr::write_csv(result$curve, output_curve_file)
readr::write_csv(result$posterior_summary, output_posterior_file)
readr::write_csv(result$diagnostics, output_diagnostics_file)
saveRDS(result, output_fit_file)

cat("Written:\n", paste0("- ", c(
  output_curve_file, output_posterior_file, output_diagnostics_file, output_fit_file
), collapse = "\n"), "\n", sep = "")
print(result$diagnostics, n = Inf)

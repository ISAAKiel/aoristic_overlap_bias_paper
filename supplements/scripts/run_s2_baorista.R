#!/usr/bin/env Rscript

# Run the Supplement 2 Baorista reconstruction and cache its outputs.
#
# Prior sensitivity runs, e.g.:
#
#   BAORISTA_OUTPUT_TAG=rate_10 \
#   BAORISTA_SIGMA_PRIOR_RATE=10 \
#   Rscript scripts/run_s2_baorista.R

source("R/baorista_helpers.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
set.seed(42)

input_file <- "outputs/s1_assigned_events.csv"
output_tag <- Sys.getenv("BAORISTA_OUTPUT_TAG", "")
output_suffix <- if (nzchar(output_tag)) paste0("_", output_tag) else ""
output_prefix <- file.path("outputs", paste0("s2_baorista", output_suffix))

output_curve_file <- paste0(output_prefix, "_curve_binned.csv")
output_posterior_file <- paste0(output_prefix, "_posterior_binned.csv")
output_diagnostics_file <- paste0(output_prefix, "_diagnostics.csv")
output_fit_file <- paste0(output_prefix, "_fit.rds")
# keep the full fit only for the baseline run
save_fit <- !nzchar(output_tag)

time_min <- -4100L
time_max <- -1800L
bin_width <- 50L
niter <- 50000L
nburnin <- 25000L
thin <- 50L
nchains <- 4L
parallel <- TRUE
sigma_prior_rate <- as.numeric(Sys.getenv("BAORISTA_SIGMA_PRIOR_RATE", "1"))
sigma_prior <- sprintf("dexp(rate=%s)", format(sigma_prior_rate, scientific = FALSE))

if (!file.exists(input_file)) stop("Render Supplement 1 first.")
intervals <- readr::read_csv(input_file, show_col_types = FALSE)

cat("Running Supplement 2 Baorista model for ", nrow(intervals), " intervals, sigma prior ",
    sigma_prior, ".\n", sep = "")

result <- run_baorista_intervals(intervals = intervals, time_min = time_min,
                                 time_max = time_max, bin_width = bin_width,
                                 niter = niter, nburnin = nburnin, thin = thin,
                                 nchains = nchains, sigma_prior = sigma_prior,
                                 parallel = parallel, seeds = seq_len(nchains),
                                 keep_fit = save_fit)

readr::write_csv(result$curve, output_curve_file)
readr::write_csv(result$posterior_summary, output_posterior_file)
readr::write_csv(result$diagnostics, output_diagnostics_file)
if (save_fit) saveRDS(result, output_fit_file)

written <- c(output_curve_file, output_posterior_file, output_diagnostics_file)
if (save_fit) written <- c(written, output_fit_file)
cat("Written:\n", paste0("- ", written, collapse = "\n"), "\n", sep = "")
print(result$diagnostics, n = Inf)

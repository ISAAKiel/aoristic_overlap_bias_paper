#!/usr/bin/env Rscript

# Supplement 3: simulations with non-uniform signals
#
# Run with defaults:
#
#   Rscript scripts/run_s3_non_uniform_signal_sensitivity.R
#
# Example with custom settings:
#
#   N_RUNS=100 N_EVENTS=1000 N_MC=100 \
#   Rscript scripts/run_s3_non_uniform_signal_sensitivity.R
#
# Required input:
#
#   outputs/s1_synthetic_phase_catalogue.csv
#
# Main outputs:
#
#   outputs/s3_results.rds
#   outputs/s3_ground_truth_curves.csv
#   outputs/s3_example_curves_run_1.csv
#   outputs/s3_sensitivity_metrics.csv
#   outputs/s3_sensitivity_summary.csv
#   outputs/s3_best_method_by_run_rmsd.csv
#   outputs/s3_best_method_by_run_overlap_bias.csv
#   outputs/s3_run_metadata.csv

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)
library(readr)
library(rlang)

# Paths
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

helper_file <- "R/aoristic_helpers.R"
source(helper_file)

phase_file <- "outputs/s1_synthetic_phase_catalogue.csv"
if (!file.exists(phase_file)) stop("Render Supplement 1 first.")

# Settings
set.seed(42)

time_min <- -4100L
time_max <- -1800L

# Intervals are treated as half-open: [from, to).
year_grid <- seq(time_min, time_max - 1, by = 1)

bin_width <- 50L
n_events <- as.integer(Sys.getenv("N_EVENTS", "1000"))
n_runs <- as.integer(Sys.getenv("N_RUNS", "100"))
n_mc <- as.integer(Sys.getenv("N_MC", "100"))

base_seed <- 700000L

signal_patterns <- s3_signal_patterns

# Load phase catalogue
phase_catalogue <- read_csv(phase_file, show_col_types = FALSE) |>
  mutate(from = as.integer(from), to = as.integer(to), width = to - from)
# Ground-truth signal curves
truth_curves <- map_dfr(
  signal_patterns,
  function(pattern) {
    probability <- make_signal_probability(
      pattern = pattern, year_grid = year_grid
    )

    make_truth_curve_binned(probability = probability, year_grid = year_grid,
                            bin_width = bin_width, bin_start = time_min) |>
      mutate(pattern = pattern, .before = 1)
  }
)

write_csv(truth_curves, "outputs/s3_ground_truth_curves.csv")

# Run simulations

s3_results <- list()
counter <- 1

for (pattern in signal_patterns) {
  cat("Running pattern:", pattern, "\n")

  for (run_id in seq_len(n_runs)) {
    s3_results[[counter]] <- make_s3_simulation_run(
      pattern = pattern,
      run_id = run_id,
      phase_catalogue = phase_catalogue,
      year_grid = year_grid,
      n_events = n_events,
      bin_width = bin_width,
      bin_start = time_min,
      n_mc = n_mc,
      base_seed = base_seed
    )

    counter <- counter + 1
  }
}

saveRDS(s3_results, "outputs/s3_results.rds")

# Summaries
s3_metrics <- map_dfr(s3_results, "metrics") |>
  recode_s3_method_names()

s3_summary <- summarise_s3_metrics(s3_metrics)

s3_best_rmsd <- make_s3_best_method_by_run(
  metrics = s3_metrics, metric = "rmsd"
)

s3_best_overlap_bias <- make_s3_best_method_by_run(
  metrics = s3_metrics, metric = "overlap_attributable_rmsd"
)

example_results <- s3_results |>
  keep(~ .x$run == 1)

example_curves <- map_dfr(
  example_results,
  function(result) {
    result$curves_binned |>
      mutate(pattern = result$pattern, .before = 1)
  }
)

run_metadata <- tibble(
  setting = c("time_min", "time_max", "bin_width", "n_events", "n_runs", "n_mc",
              "base_seed", "signal_patterns", "phase_file", "helper_file"),
  value = c(as.character(time_min), as.character(time_max),
            as.character(bin_width), as.character(n_events),
            as.character(n_runs), as.character(n_mc), as.character(base_seed),
            paste(signal_patterns, collapse = ","), phase_file, helper_file)
)

# Write outputs
write_csv(s3_metrics, "outputs/s3_sensitivity_metrics.csv")

write_csv(s3_summary, "outputs/s3_sensitivity_summary.csv")

write_csv(s3_best_rmsd, "outputs/s3_best_method_by_run_rmsd.csv")

write_csv(s3_best_overlap_bias,
          "outputs/s3_best_method_by_run_overlap_bias.csv")

write_csv(example_curves, "outputs/s3_example_curves_run_1.csv")

write_csv(run_metadata, "outputs/s3_run_metadata.csv")

cat("Supplement 3 simulation finished.\n")
print(s3_summary, n = Inf)
print(s3_best_rmsd, n = Inf)
print(s3_best_overlap_bias, n = Inf)
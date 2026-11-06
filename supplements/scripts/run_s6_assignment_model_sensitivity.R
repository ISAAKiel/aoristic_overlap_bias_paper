#!/usr/bin/env Rscript

# Supplement 6: sensitivity to the phase-assignment model
#
# M0 = equal choice among compatible labels (as in S1-S3), M1-M3 = other
# assignment models. Also computes the iterated correction (EM, exploratory).
#
# Run from the supplements/ directory:
#
#   Rscript scripts/run_s6_assignment_model_sensitivity.R
#
# Custom settings, e.g.:
#
#   N_RUNS=100 N_EVENTS=1000 \
#   Rscript scripts/run_s6_assignment_model_sensitivity.R
#
# Required input:
#
#   outputs/s1_synthetic_phase_catalogue.csv
#
# Outputs:
#
#   outputs/s6_run_metadata.csv
#   outputs/s6_assignment_level_shares.csv
#   outputs/s6_expected_metrics.csv       (infinite-sample expectation)
#   outputs/s6_metrics.csv                (one row per model/pattern/run/method)
#   outputs/s6_summary.csv
#   outputs/s6_best_method_by_run_rmsd.csv
#   outputs/s6_example_curves_run_1.csv
#
# Runtime: a few minutes with the defaults (no Monte Carlo step; the Monte
# Carlo mean has the same expectation as the standard reconstruction).
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)

# Paths
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

helper_file <- "R/aoristic_helpers.R"
source(helper_file)

phase_file <- "outputs/s1_synthetic_phase_catalogue.csv"
if (!file.exists(phase_file)) stop("Render Supplement 1 first.")

# Settings
time_min <- -4100L
time_max <- -1800L

# Intervals are treated as half-open: [from, to).
year_grid <- seq(time_min, time_max - 1L, by = 1L)
n_years <- length(year_grid)

bin_width <- 50L
n_events <- as.integer(Sys.getenv("N_EVENTS", "1000"))
n_runs <- as.integer(Sys.getenv("N_RUNS", "100"))
base_seed <- 600000L
em_iters <- c(5L, 20L)

# Relative label weights per hierarchy level for model M3.
m3_level_weights <- c("umbrella" = 0.05, "broad" = 0.15, "cross-cutting" = 0.15,
                      "medium" = 0.30, "subphase" = 0.35)

assignment_models <- c("M0_uniform", "M1_prefer_fine", "M2_prefer_coarse",
                       "M3_level_weights")

signal_patterns <- s3_signal_patterns

if (n_years %% bin_width != 0) {
  stop("The time window must be divisible by bin_width.")
}

bin_index <- rep(seq_len(n_years / bin_width), each = bin_width)
bin_years <- year_grid[seq(1, n_years, by = bin_width)]

# Phase catalogue as a label x year incidence matrix
phase_catalogue <- read_csv(phase_file, show_col_types = FALSE) |>
  mutate(from = as.integer(from), to = as.integer(to), width = to - from)

n_labels <- nrow(phase_catalogue)

incidence <- t(vapply(
  seq_len(n_labels),
  function(j) {
    as.numeric(
      year_grid >= phase_catalogue$from[j] & year_grid < phase_catalogue$to[j]
    )
  },
  numeric(n_years)
))

if (any(colSums(incidence) == 0)) {
  stop("Some years in the time window are not covered by any phase label.")
}

catalogue_overlap <- colSums(incidence)
overlap_binned <- as.numeric(tapply(catalogue_overlap, bin_index, mean))

if (!all(phase_catalogue$level %in% names(m3_level_weights))) {
  stop("Model M3 needs a weight for every level in the phase catalogue.")
}

# Assignment models: unnormalised label weights w_j (constant in time)
#
# P(label j | t) = w_j * 1(t in J_j) / sum_k w_k * 1(t in J_k)
make_label_weights <- function(model) {
  w <- switch(
    model,
    # equal choice (as in S1-S3)
    "M0_uniform" = matrix(1, n_labels, n_years),

    # w_j ~ 1 / width
    "M1_prefer_fine" = matrix(1 / phase_catalogue$width, n_labels, n_years),

    # w_j ~ width
    "M2_prefer_coarse" = matrix(phase_catalogue$width, n_labels, n_years),

    # fixed weight per hierarchy level
    "M3_level_weights" = matrix(unname(m3_level_weights[phase_catalogue$level]),
                                n_labels, n_years),

    stop("Unknown assignment model: ", model)
  )

  w * incidence
}

make_assignment_probability <- function(model) {
  w <- make_label_weights(model)
  sweep(w, 2, colSums(w), "/")
}

# Reconstructions from label counts
#
# Same results as run_standard_aoristic() / run_phase_corrected_aoristic();
# lambda counted over the observed intervals.
distribute <- function(counts, weights) {
  # weights: label x year, zero outside each label's interval
  row_totals <- rowSums(weights)
  keep <- counts > 0 & row_totals > 0
  colSums(weights[keep, , drop = FALSE] *
            (counts[keep] / row_totals[keep]))
}

reconstruct_all <- function(counts, oracle_probability) {
  observed <- counts > 0
  inc_obs <- incidence * observed

  out <- list()

  out$standard <- distribute(counts, inc_obs)

  lambda_obs <- colSums(inc_obs)
  inv_lambda <- ifelse(lambda_obs > 0, 1 / pmax(lambda_obs, 1), 0)
  w_corr <- sweep(inc_obs, 2, inv_lambda, "*")

  out$phase_corrected <- distribute(counts, w_corr)

  # iterated correction (EM); iteration 1 = phase correction
  estimate <- rep(1 / n_years, n_years)

  for (k in seq_len(max(em_iters))) {
    estimate <- distribute(counts, sweep(w_corr, 2, estimate, "*"))
    estimate <- estimate / sum(estimate)

    if (k %in% em_iters) {
      out[[paste0("em_", k)]] <- estimate
    }
  }

  # reference: flat-prior posterior under the true assignment model
  out$oracle_flat_prior <- distribute(counts, oracle_probability)

  out
}

bin_normalised <- function(x) {
  normalise_curve(as.numeric(tapply(x, bin_index, sum)))
}

curve_metrics <- function(reconstructed, truth) {
  error <- reconstructed - truth

  tibble(
    rmsd = rmsd(reconstructed, truth),
    mae = mae(reconstructed, truth),
    cor_truth = safe_cor(reconstructed, truth),
    cor_overlap = safe_cor(reconstructed, overlap_binned),
    cor_overlap_error = safe_cor(error, overlap_binned),
    overlap_attributable_rmsd = overlap_aligned_rms(error, overlap_binned),
    overlap_relative_amplitude = overlap_relative_amplitude(reconstructed, overlap_binned)
  )
}

# Level shares per model (for describing the assignment models)
level_shares <- map_dfr(assignment_models, function(model) {
  prob <- make_assignment_probability(model)

  map_dfr(unique(phase_catalogue$level), function(lv) {
    rows <- phase_catalogue$level == lv
    share <- colSums(prob[rows, , drop = FALSE])

    tibble(model = model, level = lv, year = bin_years,
           share = as.numeric(tapply(share, bin_index, mean)))
  })
})

# Expected (infinite-sample) reconstructions: isolates systematic bias
expected_metrics <- map_dfr(assignment_models, function(model) {
  prob <- make_assignment_probability(model)

  map_dfr(signal_patterns, function(pattern) {
    f <- make_signal_probability(pattern, year_grid)
    truth <- bin_normalised(f)
    expected_counts <- as.numeric(prob %*% f)

    curves <- reconstruct_all(expected_counts, prob)

    imap_dfr(curves, function(curve, method) {
      curve_metrics(bin_normalised(curve), truth) |>
        mutate(model = model, pattern = pattern, method = method,
                      .before = 1)
    })
  })
})

write_csv(expected_metrics, "outputs/s6_expected_metrics.csv")

# Repeated finite-sample simulations
simulate_label_counts <- function(f, cum_prob, seed) {
  set.seed(seed)

  t_index <- sample.int(n_years, size = n_events, replace = TRUE, prob = f)
  # guard against cumulative sums slightly below one
  u <- runif(n_events) * (1 - 1e-9)

  cp <- cum_prob[, t_index, drop = FALSE]
  labels <- colSums(cp < matrix(u, nrow = n_labels, ncol = n_events,
                                byrow = TRUE)) + 1L
  tabulate(labels, nbins = n_labels)
}

metrics_list <- list()
example_list <- list()

for (mi in seq_along(assignment_models)) {
  model <- assignment_models[mi]
  prob <- make_assignment_probability(model)
  cum_prob <- apply(prob, 2, cumsum)

  for (pi in seq_along(signal_patterns)) {
    pattern <- signal_patterns[pi]
    f <- make_signal_probability(pattern, year_grid)
    truth <- bin_normalised(f)

    message("S6: ", model, " / ", pattern)

    for (run in seq_len(n_runs)) {
      counts <- simulate_label_counts(
        f = f,
        cum_prob = cum_prob,
        seed = base_seed + mi * 100000 + pi * 1000 + run
      )

      curves <- reconstruct_all(counts, prob)

      for (method in names(curves)) {
        binned <- bin_normalised(curves[[method]])

        metrics_list[[length(metrics_list) + 1]] <-
          curve_metrics(binned, truth) |>
          mutate(model = model, pattern = pattern, run = run,
                        method = method, .before = 1)

        if (run == 1) {
          example_list[[length(example_list) + 1]] <- tibble(
            model = model,
            pattern = pattern,
            method = method,
            year = bin_years,
            value = binned,
            true = truth,
            overlap = overlap_binned
          )
        }
      }
    }
  }
}

metrics <- bind_rows(metrics_list)
example_curves <- bind_rows(example_list)

summary_metrics <- metrics |>
  group_by(model, pattern, method) |>
  summarise(
    n_runs = n(),
    rmsd_median = median(rmsd),
    rmsd_q25 = quantile(rmsd, 0.25),
    rmsd_q75 = quantile(rmsd, 0.75),
    mae_median = median(mae),
    cor_truth_median = median(cor_truth, na.rm = TRUE),
    cor_overlap_error_median = median(cor_overlap_error, na.rm = TRUE),
    overlap_attributable_rmsd_median = median(overlap_attributable_rmsd, na.rm = TRUE),
    overlap_relative_amplitude_median = median(overlap_relative_amplitude, na.rm = TRUE),
    .groups = "drop"
  )

# best method per run: standard vs. phase-corrected only
best_by_run <- metrics |>
  filter(method %in% c("standard", "phase_corrected")) |>
  group_by(model, pattern, run) |>
  slice_min(rmsd, n = 1, with_ties = FALSE) |>
  ungroup() |>
  count(model, pattern, method, name = "best_n") |>
  group_by(model, pattern) |>
  mutate(best_share = best_n / sum(best_n)) |>
  ungroup()

metadata <- tibble(
  setting = c(
    "time_min", "time_max", "bin_width", "n_events", "n_runs", "base_seed",
    "em_iters", "m3_level_weights",
    "assignment_models", "signal_patterns", "phase_file", "helper_file",
    "r_version", "run_date"
  ),
  value = c(
    time_min, time_max, bin_width, n_events, n_runs, base_seed,
    paste(em_iters, collapse = ","),
    paste(names(m3_level_weights), m3_level_weights, sep = "=", collapse = ","),
    paste(assignment_models, collapse = ","),
    paste(signal_patterns, collapse = ","),
    phase_file, helper_file,
    R.version.string, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)

write_csv(metadata, "outputs/s6_run_metadata.csv")
write_csv(level_shares, "outputs/s6_assignment_level_shares.csv")
write_csv(metrics, "outputs/s6_metrics.csv")
write_csv(summary_metrics, "outputs/s6_summary.csv")
write_csv(best_by_run, "outputs/s6_best_method_by_run_rmsd.csv")
write_csv(example_curves, "outputs/s6_example_curves_run_1.csv")

message("Supplement 6 simulation finished.")

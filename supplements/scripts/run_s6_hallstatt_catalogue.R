#!/usr/bin/env Rscript

# Supplement 6: phase correction on the Hallstatt catalogue
# with label use as in the Hallstatt data (counterpart to the Danish run)
#
# Assignment models:
#   D0_uniform          equal choice among compatible labels
#   D1_calib_standard   label weights calibrated to the observed Hallstatt label
#                       frequencies, with the standard reconstruction as
#                       event distribution
#   D1_calib_corrected  the same, with the phase-corrected reconstruction
# (not possible under a uniform event process; see S6)
#
# Methods: standard, phase_corrected, oracle_flat_prior (true assignment
# probabilities)
#
# The signal patterns of S3 are defined on the Danish window (4100-1800 BCE).
# They are mapped linearly onto the Hallstatt window, so that each pattern
# has the same shape relative to the window.
#
# Run from the supplements/ directory:
#
#   Rscript scripts/run_s6_hallstatt_catalogue.R
#
# Custom settings, e.g.:
#
#   N_RUNS=100 N_EVENTS=2242 Rscript scripts/run_s6_hallstatt_catalogue.R
#
# Required input:
#
#   outputs/s5_hallstatt_interval_counts.csv   (from Supplement 5)
#
# Outputs:
#
#   outputs/s6_hallstatt_run_metadata.csv
#   outputs/s6_hallstatt_label_weights.csv
#   outputs/s6_hallstatt_expected_metrics.csv   (infinite-sample expectation)
#   outputs/s6_hallstatt_metrics.csv            (one row per model/pattern/run/method)
#   outputs/s6_hallstatt_summary.csv
#   outputs/s6_hallstatt_relative_rmsd.csv      (median RMSD relative to standard)
#   outputs/s6_hallstatt_label_use.csv
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)

# Paths
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

helper_file <- "R/aoristic_helpers.R"
source(helper_file)

interval_file <- "outputs/s5_hallstatt_interval_counts.csv"
if (!file.exists(interval_file)) stop("Render Supplement 5 first.")

# Settings
time_min <- -780L
time_max <- -450L

# Intervals are treated as half-open: [start, end).
year_grid <- seq(time_min, time_max - 1L, by = 1L)
n_years <- length(year_grid)

bin_width <- 10L
n_runs <- as.integer(Sys.getenv("N_RUNS", "100"))
base_seed <- 700000L
calibration_iterations <- 3000L

# The model order fixes the random seeds; keep it when adding models.
assignment_models <- c("D0_uniform", "D1_calib_standard", "D1_calib_corrected")
signal_patterns <- s3_signal_patterns

if (n_years %% bin_width != 0) {
  stop("The time window must be divisible by bin_width.")
}

bin_index <- rep(seq_len(n_years / bin_width), each = bin_width)
bin_years <- year_grid[seq(1, n_years, by = bin_width)]

# Years on the Danish window, used only to evaluate the signal patterns
signal_years <- -4100 + (year_grid - time_min) * 2300 / n_years

# Hallstatt catalogue as a label x year incidence matrix
catalogue <- read_csv(interval_file, show_col_types = FALSE) |>
  transmute(start = as.integer(start), end = as.integer(end),
            width = end - start, n_observed = as.integer(n)) |>
  filter(width > 0)

n_labels <- nrow(catalogue)

# By default, simulate as many events as there are Hallstatt records.
n_events <- as.integer(Sys.getenv("N_EVENTS", sum(catalogue$n_observed)))

incidence <- t(vapply(
  seq_len(n_labels),
  function(j) {
    as.numeric(year_grid >= catalogue$start[j] & year_grid < catalogue$end[j])
  },
  numeric(n_years)
))

if (any(colSums(incidence) == 0)) {
  stop("Some years in the time window are not covered by any Hallstatt interval.")
}

catalogue_overlap <- colSums(incidence)
overlap_binned <- as.numeric(tapply(catalogue_overlap, bin_index, mean))

# Label propensities
assignment_probability_from_weights <- function(w_label) {
  w <- incidence * w_label
  sweep(w, 2, colSums(w), "/")
}

# reference event distributions for the calibration
distribute_annual <- function(counts, weights) {
  row_totals <- rowSums(weights)
  keep <- counts > 0 & row_totals > 0
  colSums(weights[keep, , drop = FALSE] * (counts[keep] / row_totals[keep]))
}

f_ref_standard <- normalise_curve(distribute_annual(catalogue$n_observed, incidence))
f_ref_corrected <- normalise_curve(distribute_annual(
  catalogue$n_observed, sweep(incidence, 2, 1 / catalogue_overlap, "*")
))

# D1: fixed-point update of w until expected = observed label frequencies
calibrate_weights <- function(target, f_ref, iterations) {
  target_share <- target / sum(target)
  w <- rep(1, n_labels)

  for (i in seq_len(iterations)) {
    expected_share <- as.numeric(assignment_probability_from_weights(w) %*% f_ref)
    w <- w * target_share / pmax(expected_share, 1e-300)
    w <- w / sum(w)
  }

  expected_share <- as.numeric(assignment_probability_from_weights(w) %*% f_ref)

  list(weights = w,
       max_abs_share_error = max(abs(expected_share - target_share)))
}

calibration_standard <- calibrate_weights(
  catalogue$n_observed, f_ref_standard, calibration_iterations
)
calibration_corrected <- calibrate_weights(
  catalogue$n_observed, f_ref_corrected, calibration_iterations
)

if (!is.finite(calibration_standard$max_abs_share_error) ||
    !is.finite(calibration_corrected$max_abs_share_error)) {
  stop("Calibration of the label propensities failed.")
}

calibration <- list(
  max_abs_share_error = max(calibration_standard$max_abs_share_error,
                            calibration_corrected$max_abs_share_error)
)

label_weights <- list(D0_uniform = rep(1, n_labels),
                      D1_calib_standard = calibration_standard$weights,
                      D1_calib_corrected = calibration_corrected$weights)

assignment_probability <- map(label_weights, assignment_probability_from_weights)

message("D1 calibration: max. absolute deviation of expected from observed ",
        "label shares = ", signif(calibration$max_abs_share_error, 3))

# Reconstructions from label counts (as in the main S6 script)
distribute <- function(counts, weights) {
  row_totals <- rowSums(weights)
  keep <- counts > 0 & row_totals > 0
  colSums(weights[keep, , drop = FALSE] *
            (counts[keep] / row_totals[keep]))
}

# lambda over the observed intervals, as in aoristAAR
phase_corrected_weights <- function(counts) {
  inc_obs <- incidence * (counts > 0)
  lambda <- colSums(inc_obs)
  inv_lambda <- ifelse(lambda > 0, 1 / pmax(lambda, 1), 0)

  sweep(inc_obs, 2, inv_lambda, "*")
}

reconstruct_all <- function(counts, oracle_probability) {
  inc_obs <- incidence * (counts > 0)

  out <- list()
  out$standard <- distribute(counts, inc_obs)
  out$phase_corrected <- distribute(counts, phase_corrected_weights(counts))

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
    cor_truth = safe_cor(reconstructed, truth),
    cor_overlap = safe_cor(reconstructed, overlap_binned),
    overlap_attributable_rmsd = overlap_aligned_rms(error, overlap_binned),
    overlap_relative_amplitude = overlap_relative_amplitude(reconstructed, overlap_binned)
  )
}

# Expected (infinite-sample) reconstructions
#
# Expected counts are scaled to n_events, as in the finite-sample runs below.
expected_metrics <- map_dfr(assignment_models, function(model) {
  prob <- assignment_probability[[model]]

  map_dfr(signal_patterns, function(pattern) {
    f <- make_signal_probability(pattern, signal_years)
    truth <- bin_normalised(f)
    expected_counts <- n_events * as.numeric(prob %*% (f / sum(f)))

    curves <- reconstruct_all(expected_counts, prob)

    imap_dfr(curves, function(curve, method) {
      curve_metrics(bin_normalised(curve), truth) |>
        mutate(model = model, pattern = pattern, method = method,
                      .before = 1)
    })
  })
})

write_csv(expected_metrics, "outputs/s6_hallstatt_expected_metrics.csv")

# Repeated finite-sample simulations
simulate_label_counts <- function(f, cum_prob, seed) {
  set.seed(seed)

  t_index <- sample.int(n_years, size = n_events, replace = TRUE, prob = f)
  u <- runif(n_events) * (1 - 1e-9)

  cp <- cum_prob[, t_index, drop = FALSE]
  labels <- colSums(cp < matrix(u, nrow = n_labels, ncol = n_events,
                                byrow = TRUE)) + 1L
  tabulate(labels, nbins = n_labels)
}

metrics_list <- list()
label_use_list <- list()

for (mi in seq_along(assignment_models)) {
  model <- assignment_models[mi]
  prob <- assignment_probability[[model]]
  cum_prob <- apply(prob, 2, cumsum)

  for (pi in seq_along(signal_patterns)) {
    pattern <- signal_patterns[pi]
    f <- make_signal_probability(pattern, signal_years)
    truth <- bin_normalised(f)

    message("S6 Hallstatt catalogue: ", model, " / ", pattern)

    for (run in seq_len(n_runs)) {
      counts <- simulate_label_counts(
        f = f,
        cum_prob = cum_prob,
        seed = base_seed + mi * 100000 + pi * 1000 + run
      )

      label_use_list[[length(label_use_list) + 1]] <- tibble(
        model = model,
        pattern = pattern,
        run = run,
        n_labels_observed = sum(counts > 0),
        n_labels_used_once = sum(counts == 1),
        share_top2 = sum(sort(counts, decreasing = TRUE)[1:2]) / sum(counts)
      )

      curves <- reconstruct_all(counts, prob)

      for (method in names(curves)) {
        metrics_list[[length(metrics_list) + 1]] <-
          curve_metrics(bin_normalised(curves[[method]]), truth) |>
          mutate(model = model, pattern = pattern, run = run,
                        method = method, .before = 1)
      }
    }
  }
}

metrics <- bind_rows(metrics_list)
label_use <- bind_rows(label_use_list)

summary_metrics <- metrics |>
  group_by(model, pattern, method) |>
  summarise(
    n_runs = n(),
    rmsd_median = median(rmsd),
    rmsd_q25 = quantile(rmsd, 0.25),
    rmsd_q75 = quantile(rmsd, 0.75),
    cor_truth_median = median(cor_truth, na.rm = TRUE),
    overlap_attributable_rmsd_median = median(overlap_attributable_rmsd, na.rm = TRUE),
    overlap_relative_amplitude_median = median(overlap_relative_amplitude, na.rm = TRUE),
    .groups = "drop"
  )

# median RMSD relative to standard, share of runs better than standard
relative_rmsd <- metrics |>
  select(model, pattern, run, method, rmsd) |>
  left_join(
    metrics |>
      filter(method == "standard") |>
      select(model, pattern, run, rmsd_standard = rmsd),
    by = c("model", "pattern", "run")
  ) |>
  filter(method != "standard") |>
  group_by(model, pattern, method) |>
  summarise(relative_rmsd_median = median(rmsd / rmsd_standard),
            share_better_than_standard = mean(rmsd < rmsd_standard),
            .groups = "drop")

label_use_summary <- label_use |>
  group_by(model) |>
  summarise(n_labels_observed_median = median(n_labels_observed),
            n_labels_used_once_median = median(n_labels_used_once),
            share_top2_median = median(share_top2), .groups = "drop")

metadata <- tibble(
  setting = c(
    "time_min", "time_max", "bin_width", "n_events", "n_runs", "base_seed",
    "calibration_iterations", "calibration_max_abs_share_error",
    "n_labels", "observed_n_labels_used_once", "observed_share_top2",
    "assignment_models", "signal_patterns", "interval_file", "helper_file",
    "r_version", "run_date"
  ),
  value = c(
    time_min, time_max, bin_width, n_events, n_runs, base_seed,
    calibration_iterations,
    signif(calibration$max_abs_share_error, 3),
    n_labels, sum(catalogue$n_observed == 1),
    round(sum(sort(catalogue$n_observed, decreasing = TRUE)[1:2]) /
            sum(catalogue$n_observed), 3),
    paste(assignment_models, collapse = ","),
    paste(signal_patterns, collapse = ","),
    interval_file, helper_file,
    R.version.string, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)

label_weight_table <- catalogue |>
  mutate(w_D1_calib_standard = label_weights$D1_calib_standard,
         w_D1_calib_corrected = label_weights$D1_calib_corrected)

write_csv(metadata, "outputs/s6_hallstatt_run_metadata.csv")
write_csv(label_weight_table, "outputs/s6_hallstatt_label_weights.csv")
write_csv(metrics, "outputs/s6_hallstatt_metrics.csv")
write_csv(summary_metrics, "outputs/s6_hallstatt_summary.csv")
write_csv(relative_rmsd, "outputs/s6_hallstatt_relative_rmsd.csv")
write_csv(label_use_summary, "outputs/s6_hallstatt_label_use.csv")

message("Supplement 6 Hallstatt-catalogue simulation finished.")

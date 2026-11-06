#!/usr/bin/env Rscript

# Supplement 3: selected Baorista benchmark runs
# Fitting and summaries: R/baorista_helpers.R
helper_file <- "R/aoristic_helpers.R"
baorista_helper_file <- "R/baorista_helpers.R"
source(helper_file)
source(baorista_helper_file)

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
cache_dir <- "outputs/s3_baorista_cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

s3_results_file <- "outputs/s3_results.rds"
if (!file.exists(s3_results_file)) stop("Run scripts/run_s3_non_uniform_signal_sensitivity.R first.")

output_curve_file <- "outputs/s3_baorista_curves.csv"
output_metrics_file <- "outputs/s3_baorista_metrics.csv"
output_metric_posterior_file <- "outputs/s3_baorista_metric_posterior.csv"
output_metric_summary_file <- "outputs/s3_baorista_metric_posterior_summary.csv"
output_diagnostics_file <- "outputs/s3_baorista_diagnostics.csv"
output_metadata_file <- "outputs/s3_baorista_run_metadata.csv"

time_min <- -4100L
time_max <- -1800L
bin_width <- 50L
niter <- 50000L
nburnin <- 25000L
thin <- 50L
nchains <- 4L
parallel <- TRUE
sigma_prior <- "dexp(rate=1)"
refresh <- Sys.getenv("REFRESH_BAORISTA", "false") == "true"
save_fits <- FALSE

calculate_baorista_metric_posterior <- function(posterior_draws, years, curves_binned) {
  evaluation <- curves_binned |>
    dplyr::select(year, true, overlap) |>
    dplyr::filter(year %in% years) |>
    dplyr::arrange(match(year, years))

  if (nrow(evaluation) != ncol(posterior_draws)) {
    stop("S3 evaluation grid does not match Baorista posterior blocks.")
  }

  purrr::map_dfr(seq_len(nrow(posterior_draws)), function(draw_id) {
    input <- evaluation |>
      dplyr::mutate(baorista = as.numeric(posterior_draws[draw_id, ]))

    fit_metrics <- compute_curve_metrics(input, truth_col = "true") |>
      dplyr::filter(method == "baorista") |>
      dplyr::select(rmsd, mae, correlation)

    overlap_metrics <- compute_overlap_diagnostics(
      input, truth_col = "true", overlap_col = "overlap"
    ) |>
      dplyr::filter(method == "baorista") |>
      dplyr::select(overlap_curve_correlation, overlap_error_correlation,
                    overlap_abs_error_correlation,
                    overlap_squared_error_correlation,
                    overlap_attributable_rmsd)

    dplyr::bind_cols(tibble::tibble(draw = draw_id), fit_metrics, overlap_metrics)
  })
}

summarise_metric_posterior <- function(metric_posterior) {
  metric_posterior |>
    tidyr::pivot_longer(
      cols = c(rmsd, mae, correlation, overlap_attributable_rmsd),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::group_by(pattern, run, metric) |>
    dplyr::summarise(
      n_posterior_draws = sum(is.finite(value)),
      posterior_mean = ifelse(any(is.finite(value)), mean(value, na.rm = TRUE), NA_real_),
      posterior_median = ifelse(any(is.finite(value)), median(value, na.rm = TRUE), NA_real_),
      posterior_q025 = ifelse(any(is.finite(value)), quantile(value, 0.025, na.rm = TRUE, names = FALSE), NA_real_),
      posterior_q975 = ifelse(any(is.finite(value)), quantile(value, 0.975, na.rm = TRUE, names = FALSE), NA_real_),
      posterior_q10 = ifelse(any(is.finite(value)), quantile(value, 0.10, na.rm = TRUE, names = FALSE), NA_real_),
      posterior_q90 = ifelse(any(is.finite(value)), quantile(value, 0.90, na.rm = TRUE, names = FALSE), NA_real_),
      .groups = "drop"
    )
}

run_one <- function(result, index) {
  pattern <- result$pattern
  run_id <- result$run
  stem <- paste0("s3_baorista_", safe_slug(pattern), "_run_", run_id)
  paths <- list(
    curve = file.path(cache_dir, paste0(stem, "_curve.csv")),
    posterior = file.path(cache_dir, paste0(stem, "_posterior.csv")),
    diagnostics = file.path(cache_dir, paste0(stem, "_diagnostics.csv")),
    metrics = file.path(cache_dir, paste0(stem, "_metric_posterior.csv")),
    fit = file.path(cache_dir, paste0(stem, "_fit.rds"))
  )

  cache_ready <- all(file.exists(unlist(paths[c("curve", "posterior", "diagnostics", "metrics")])))

  # refit caches without overlap_attributable_rmsd
  if (cache_ready) {
    cached_cols <- names(readr::read_csv(paths$metrics, n_max = 1, show_col_types = FALSE))
    if (!"overlap_attributable_rmsd" %in% cached_cols) {
      cat("Cached metrics for ", pattern, " run ", run_id,
          " predate overlap_attributable_rmsd; refitting.\n", sep = "")
      cache_ready <- FALSE
    }
  }

  if (!refresh && cache_ready) {
    cat("Using cached Baorista result: ", pattern, " run ", run_id, "\n", sep = "")
    curve <- readr::read_csv(paths$curve, show_col_types = FALSE)
    diagnostics <- readr::read_csv(paths$diagnostics, show_col_types = FALSE)
    metric_posterior <- readr::read_csv(paths$metrics, show_col_types = FALSE)
  } else {
    cat("Running Baorista: ", pattern, " run ", run_id, "\n", sep = "")
    model <- run_baorista_intervals(
      intervals = result$assigned_events,
      time_min = time_min,
      time_max = time_max,
      bin_width = bin_width,
      niter = niter,
      nburnin = nburnin,
      thin = thin,
      nchains = nchains,
      sigma_prior = sigma_prior,
      parallel = parallel,
      seeds = 900000 + index * 1000 + seq_len(nchains),
      keep_fit = save_fits
    )

    curve <- model$curve
    diagnostics <- model$diagnostics
    metric_posterior <- calculate_baorista_metric_posterior(
      model$posterior_draws, model$years, result$curves_binned
    )

    readr::write_csv(curve, paths$curve)
    readr::write_csv(model$posterior_summary, paths$posterior)
    readr::write_csv(diagnostics, paths$diagnostics)
    readr::write_csv(metric_posterior, paths$metrics)
    if (save_fits) saveRDS(model, paths$fit)
  }

  curve <- curve |> dplyr::mutate(pattern = pattern, run = run_id, .before = 1)
  diagnostics <- diagnostics |> dplyr::mutate(pattern = pattern, run = run_id, .before = 1)
  metric_posterior <- metric_posterior |> dplyr::mutate(pattern = pattern, run = run_id, .before = 1)

  metric_input <- result$curves_binned |>
    dplyr::select(year, true, overlap) |>
    dplyr::left_join(curve |> dplyr::select(year, baorista), by = "year") |>
    dplyr::mutate(baorista = normalise_curve(tidyr::replace_na(baorista, 0)))

  fit_metrics <- compute_curve_metrics(metric_input, truth_col = "true") |>
    dplyr::filter(method == "baorista") |>
    dplyr::select(method, rmsd, mae, correlation)
  overlap_metrics <- compute_overlap_diagnostics(
    metric_input, truth_col = "true", overlap_col = "overlap"
  ) |>
    dplyr::filter(method == "baorista") |>
    dplyr::select(
      method, overlap_curve_correlation, overlap_error_correlation,
      overlap_abs_error_correlation, overlap_squared_error_correlation,
      overlap_attributable_rmsd
    )

  metrics <- fit_metrics |>
    dplyr::left_join(overlap_metrics, by = "method") |>
    dplyr::mutate(method = "Baorista", pattern = pattern, run = run_id, .before = 1)

  list(curves = curve, metrics = metrics, metric_posterior = metric_posterior,
       metric_summary = summarise_metric_posterior(metric_posterior),
       diagnostics = diagnostics)
}

s3_results <- readRDS(s3_results_file)
available_patterns <- unique(purrr::map_chr(s3_results, "pattern"))
selected_runs <- 1L
selected_patterns <- available_patterns

selected <- purrr::keep(
  s3_results,
  ~ .x$pattern %in% selected_patterns && .x$run %in% selected_runs
)
if (length(selected) == 0) stop("No S3 results match the selected patterns and runs.")

results <- purrr::imap(selected, run_one)
curves <- purrr::map_dfr(results, "curves")
metrics <- purrr::map_dfr(results, "metrics")
metric_posterior <- purrr::map_dfr(results, "metric_posterior")
metric_summary <- purrr::map_dfr(results, "metric_summary")
diagnostics <- purrr::map_dfr(results, "diagnostics")

metadata <- tibble::tibble(
  setting = c(
    "time_min", "time_max", "bin_width", "selected_patterns", "selected_run_ids",
    "n_selected_fits", "niter", "nburnin", "thin", "nchains", "parallel",
    "sigma_prior", "s3_results_file", "aoristic_helper_file",
    "baorista_helper_file", "cache_dir", "refresh_baorista", "save_fits"
  ),
  value = c(
    time_min, time_max, bin_width, paste(selected_patterns, collapse = ","),
    paste(selected_runs, collapse = ","), length(selected), niter, nburnin,
    thin, nchains, parallel, sigma_prior, s3_results_file, helper_file,
    baorista_helper_file, cache_dir, refresh, save_fits
  ) |> as.character()
)

readr::write_csv(curves, output_curve_file)
readr::write_csv(metrics, output_metrics_file)
readr::write_csv(metric_posterior, output_metric_posterior_file)
readr::write_csv(metric_summary, output_metric_summary_file)
readr::write_csv(diagnostics, output_diagnostics_file)
readr::write_csv(metadata, output_metadata_file)

run1_curves <- dplyr::filter(curves, run == 1)
run1_metrics <- dplyr::filter(metrics, run == 1)
run1_summary <- dplyr::filter(metric_summary, run == 1)
if (nrow(run1_curves) > 0) readr::write_csv(run1_curves, "outputs/s3_baorista_example_curves_run_1.csv")
if (nrow(run1_metrics) > 0) readr::write_csv(run1_metrics, "outputs/s3_baorista_example_metrics_run_1.csv")
if (nrow(run1_summary) > 0) readr::write_csv(run1_summary, "outputs/s3_baorista_example_metric_posterior_summary_run_1.csv")

cat("Supplement 3 Baorista benchmark finished.\n")
print(metrics, n = Inf)

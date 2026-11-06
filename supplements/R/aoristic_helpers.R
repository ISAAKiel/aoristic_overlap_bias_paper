# Shared helper functions for aoristic simulation supplements

# Tables as LaTeX floats in the PDF, so that they do not break across pages
kbl <- function(d, ...) {
  if (knitr::is_latex_output()) {
    tab <- knitr::kable(d, format = "latex", booktabs = TRUE, linesep = "", ...)
    structure(paste0("\\footnotesize\n", tab), format = "latex",
              class = "knitr_kable")
  } else {
    knitr::kable(d, ...)
  }
}
normalise_curve <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0

  total <- sum(x, na.rm = TRUE)

  if (is.na(total) || total <= 0) {
    return(rep(0, length(x)))
  }

  x / total
}

bin_year <- function(year, bin_width = 50, bin_start = min(year, na.rm = TRUE)) {
  floor((year - bin_start) / bin_width) * bin_width + bin_start
}

safe_cor <- function(x, y) {
  if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }

  stats::cor(x, y, use = "complete.obs")
}

rmsd <- function(observed, expected) {
  sqrt(mean((observed - expected)^2, na.rm = TRUE))
}

mae <- function(observed, expected) {
  mean(abs(observed - expected), na.rm = TRUE)
}

# Magnitude-based overlap diagnostics (population sd, as for RMSD)
sd_pop <- function(x) {
  x <- x[is.finite(x)]
  sqrt(mean((x - mean(x))^2))
}

# |cor(x, overlap)| * sd(x); applied to the error: overlap-attributable RMSD
overlap_aligned_rms <- function(x, overlap) {
  r <- safe_cor(x, overlap)
  if (is.na(r)) return(0)
  abs(r) * sd_pop(x)
}

# |beta| * sd(overlap) / mean(x): overlap-aligned relative amplitude
overlap_relative_amplitude <- function(x, overlap) {
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(m) || m == 0) return(NA_real_)
  overlap_aligned_rms(x, overlap) / m
}

# Density against overlap, one panel per method, common y axis
plot_overlap_scatter <- function(curves, method_cols, overlap_col = "overlap",
                                 method_labels = NULL, colours = NULL) {
  long <- curves |>
    dplyr::select(dplyr::all_of(c(overlap_col, method_cols))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(method_cols), names_to = "method",
                        values_to = "density") |>
    dplyr::rename(overlap = dplyr::all_of(overlap_col))

  if (!is.null(method_labels)) {
    long <- long |>
      dplyr::mutate(method = factor(
        dplyr::recode(method, !!!method_labels),
        levels = unname(method_labels[method_cols])
      ))
  }

  p <- ggplot2::ggplot(long, ggplot2::aes(overlap, density, colour = method)) +
    ggplot2::geom_hline(yintercept = 1 / nrow(curves), colour = "grey70",
                        linetype = "dashed") +
    ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    ggplot2::facet_wrap(~method, nrow = 1) +
    ggplot2::labs(
      x = "Overlap (number of intervals covering the bin)",
      y = "Reconstructed density per bin",
      caption = "Dashed line: density of a perfectly flat curve. Common y axis across panels."
    ) +
    plot_theme() +
    ggplot2::theme(legend.position = "none")

  if (!is.null(colours)) {
    p <- p + ggplot2::scale_colour_manual(values = colours)
  }

  p
}

plot_theme <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   plot.title.position = "plot", plot.caption.position = "plot")
}

# Phase catalogue helpers

count_phase_overlap <- function(phase_catalogue, year_grid) {
  tibble::tibble(year = year_grid) |>
    dplyr::mutate(
      overlap = purrr::map_int(
        year,
        function(y) {
          sum(phase_catalogue$from <= y & phase_catalogue$to > y)
        }
      )
    )
}

bin_overlap_curve <- function(overlap_curve, bin_width, bin_start) {
  overlap_curve |>
    dplyr::mutate(
      bin = bin_year(year, bin_width = bin_width, bin_start = bin_start)
    ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(overlap = mean(overlap), .groups = "drop") |>
    dplyr::rename(year = bin)
}

make_phase_catalogue_overlap_curve <- function(phase_catalogue, year_grid) {
  required_cols <- c("from", "to")

  if (!all(required_cols %in% names(phase_catalogue))) {
    stop("phase_catalogue must contain columns: ",
         paste(required_cols, collapse = ", "))
  }

  tibble::tibble(year = year_grid) |>
    dplyr::mutate(
      overlap = purrr::map_int(
        year,
        function(y) {
          sum(phase_catalogue$from <= y & phase_catalogue$to > y, na.rm = TRUE)
        }
      )
    )
}

# Simulation helpers


bin_curve <- function(data, value_col, bin_width, bin_start) {
  value_col <- rlang::ensym(value_col)

  data |>
    dplyr::mutate(
      bin = bin_year(year, bin_width = bin_width, bin_start = bin_start)
    ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(value = sum(!!value_col, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(year = bin) |>
    dplyr::mutate(value = normalise_curve(value))
}

assign_phase <- function(true_date, phase_catalogue) {
  candidates <- phase_catalogue |>
    dplyr::filter(from <= true_date, to > true_date)

  if (nrow(candidates) == 0) {
    stop("No compatible phase found for true_date = ", true_date,
         ". Check whether phase intervals are half-open [from, to).")
  }

  candidates[sample(seq_len(nrow(candidates)), size = 1), ]
}

assign_events_to_phases <- function(events, phase_catalogue) {
  purrr::map_dfr(
    seq_len(nrow(events)),
    function(i) {
      chosen_phase <- assign_phase(events$true_date[i], phase_catalogue)

      tibble::tibble(id = events$id[i], true_date = events$true_date[i],
                     phase = chosen_phase$name, level = chosen_phase$level,
                     start = chosen_phase$from, end = chosen_phase$to,
                     phase_width = chosen_phase$width)
    }
  )
}

# Aoristic reconstruction helpers
# Standard and phase-corrected reconstruction via aoristAAR
# (aorist(method = "weight") / aorist(method = "period_correction")).
# aoristAAR uses inclusive ends, so intervals are passed as [start, end - 1],
# after clipping to the grid.
# AORISTIC_BACKEND=internal switches to the internal implementation
# (see scripts/check_aoristAAR_equivalence.R).
aoristic_backend <- function() {
  backend <- tolower(trimws(Sys.getenv("AORISTIC_BACKEND", "aoristAAR")))
  if (!backend %in% c("aoristaar", "internal")) {
    stop("AORISTIC_BACKEND must be 'aoristAAR' or 'internal'.")
  }
  backend
}

check_aoristAAR <- function(min_version = "0.1.1") {
  if (!requireNamespace("aoristAAR", quietly = TRUE)) {
    stop("Package aoristAAR is not installed. Install it with\n",
         "  remotes::install_github(\"ISAAKiel/aoristAAR\")\n",
         "or set AORISTIC_BACKEND=internal.")
  }
  if (utils::packageVersion("aoristAAR") < min_version) {
    stop("aoristAAR >= ", min_version, " is required (installed: ",
         utils::packageVersion("aoristAAR"), "). Older versions overweight ",
         "short intervals in method = \"weight\".")
  }
  invisible(TRUE)
}

run_aoristAAR_on_grid <- function(intervals, year_grid, method) {
  check_aoristAAR()

  grid_min <- min(year_grid)
  grid_max <- max(year_grid)

  input <- data.frame(from = pmax(as.numeric(intervals$start), grid_min),
                      to = pmin(as.numeric(intervals$end) - 1, grid_max))
  input <- input[is.finite(input$from) & is.finite(input$to) & input$to >= input$from, ,
                 drop = FALSE]

  out <- numeric(length(year_grid))
  if (nrow(input) == 0) return(out)

  res <- aoristAAR::aorist(input, from = "from", to = "to", method = method)
  idx <- match(year_grid, res$date)
  ok <- !is.na(idx)
  out[ok] <- res$sum[idx[ok]]
  out[!is.finite(out)] <- 0
  out
}

run_standard_aoristic <- function(intervals, year_grid) {
  if (aoristic_backend() == "internal") {
    return(run_standard_aoristic_internal(intervals, year_grid))
  }

  tibble::tibble(
    year = year_grid,
    standard = normalise_curve(run_aoristAAR_on_grid(intervals, year_grid, "weight"))
  )
}

run_phase_corrected_aoristic <- function(intervals, year_grid) {
  if (aoristic_backend() == "internal") {
    return(run_phase_corrected_aoristic_internal(intervals, year_grid))
  }

  tibble::tibble(
    year = year_grid,
    phase_corrected = normalise_curve(
      run_aoristAAR_on_grid(intervals, year_grid, "period_correction")
    )
  )
}

# Internal reference implementations (half-open intervals).

run_standard_aoristic_internal <- function(intervals, year_grid) {
  out <- tibble::tibble(year = year_grid, standard = 0)

  for (i in seq_len(nrow(intervals))) {
    inside <- out$year >= intervals$start[i] & out$year < intervals$end[i]
    n_inside <- sum(inside)

    if (n_inside > 0) {
      out$standard[inside] <- out$standard[inside] + 1 / n_inside
    }
  }

  out |>
    dplyr::mutate(standard = normalise_curve(standard))
}




# Additional reconstruction helpers for Supplement 2 and 3
observed_interval_system <- function(intervals) {
  intervals |>
    dplyr::distinct(start, end) |>
    dplyr::arrange(start, end) |>
    dplyr::mutate(observed_interval_id = dplyr::row_number(),
                  width = end - start)
}

run_phase_corrected_aoristic_internal <- function(intervals, year_grid) {
  # weights 1/lambda(t), lambda = number of unique observed intervals at t;
  # each observation sums to one

  observed_intervals <- observed_interval_system(intervals)

  overlap <- tibble::tibble(year = year_grid, n_intervals = 0)

  for (i in seq_len(nrow(observed_intervals))) {
    inside <- overlap$year >= observed_intervals$start[i] &
      overlap$year < observed_intervals$end[i]

    overlap$n_intervals[inside] <- overlap$n_intervals[inside] + 1
  }

  out <- tibble::tibble(year = year_grid, phase_corrected = 0)

  for (i in seq_len(nrow(intervals))) {
    inside <- out$year >= intervals$start[i] & out$year < intervals$end[i]

    local_overlap <- overlap$n_intervals[inside]
    local_weights <- ifelse(local_overlap > 0, 1 / local_overlap, 0)
    local_total <- sum(local_weights, na.rm = TRUE)

    if (local_total > 0) {
      out$phase_corrected[inside] <- out$phase_corrected[inside] +
        local_weights / local_total
    }
  }

  out |>
    dplyr::mutate(phase_corrected = normalise_curve(phase_corrected))
}

make_reconstruction_binned <- function(curve, value_col, output_col, bin_width, bin_start) {
  value_col <- rlang::ensym(value_col)

  curve |>
    dplyr::select(year, value = !!value_col) |>
    bin_curve(value_col = value, bin_width = bin_width,
              bin_start = bin_start) |>
    dplyr::rename(!!output_col := value)
}

run_monte_carlo_aoristic <- function(intervals, year_grid, bin_width, bin_start,
                                     n_iter = 1000, seed = 42) {
  # one random year per observation and iteration; returns mean, median and
  # 90% interval per bin

  set.seed(seed)

  bin_levels <- sort(unique(bin_year(year_grid, bin_width = bin_width,
                                     bin_start = bin_start)))

  mc_matrix <- matrix(
    0,
    nrow = length(bin_levels),
    ncol = n_iter,
    dimnames = list(as.character(bin_levels), paste0("iter_", seq_len(n_iter)))
  )

  for (iter in seq_len(n_iter)) {
    sampled_years <- purrr::map_int(seq_len(nrow(intervals)), function(i) {
      possible_years <- year_grid[
        year_grid >= intervals$start[i] &
          year_grid < intervals$end[i]
      ]

      if (length(possible_years) == 0) {
        return(NA_integer_)
      }

      sample(possible_years, size = 1)
    })

    sampled_curve <- tibble::tibble(year = sampled_years) |>
      dplyr::filter(!is.na(year)) |>
      dplyr::count(year, name = "count") |>
      dplyr::right_join(tibble::tibble(year = year_grid), by = "year") |>
      dplyr::arrange(year) |>
      dplyr::mutate(count = tidyr::replace_na(count, 0),
                    count = normalise_curve(count))

    sampled_binned <- sampled_curve |>
      make_reconstruction_binned(value_col = count, output_col = "value",
                                 bin_width = bin_width,
                                 bin_start = bin_start) |>
      dplyr::right_join(tibble::tibble(year = bin_levels), by = "year") |>
      dplyr::arrange(year) |>
      dplyr::mutate(value = tidyr::replace_na(value, 0),
                    value = normalise_curve(value))

    mc_matrix[, iter] <- sampled_binned$value
  }

  tibble::tibble(
    year = bin_levels,
    mc_mean = rowMeans(mc_matrix, na.rm = TRUE),
    mc_median = apply(mc_matrix, 1, stats::median, na.rm = TRUE),
    mc_q05 = apply(mc_matrix, 1, stats::quantile, probs = 0.05, na.rm = TRUE),
    mc_q95 = apply(mc_matrix, 1, stats::quantile, probs = 0.95, na.rm = TRUE)
  ) |>
    dplyr::mutate(mc_mean = normalise_curve(mc_mean),
                  mc_median = normalise_curve(mc_median))
}


compute_curve_metrics <- function(curves, truth_col = "true") {
  method_cols <- setdiff(names(curves), c("year", truth_col, "overlap"))

  purrr::map_dfr(method_cols, function(method_col) {
    tibble::tibble(
      method = method_col,
      rmsd = rmsd(curves[[method_col]], curves[[truth_col]]),
      mae = mae(curves[[method_col]], curves[[truth_col]]),
      correlation = safe_cor(curves[[method_col]], curves[[truth_col]])
    )
  })
}

compute_overlap_diagnostics <- function(curves, truth_col = "true", overlap_col = "overlap") {
  method_cols <- setdiff(names(curves), c("year", truth_col, overlap_col))

  purrr::map_dfr(method_cols, function(method_col) {
    error <- curves[[method_col]] - curves[[truth_col]]

    tibble::tibble(
      method = method_col,
      overlap_curve_correlation = safe_cor(curves[[overlap_col]], curves[[method_col]]),
      overlap_error_correlation = safe_cor(curves[[overlap_col]], error),
      overlap_abs_error_correlation = safe_cor(curves[[overlap_col]], abs(error)),
      overlap_squared_error_correlation = safe_cor(curves[[overlap_col]], error^2),
      overlap_attributable_rmsd = overlap_aligned_rms(error, curves[[overlap_col]]),
      overlap_attributable_share = overlap_aligned_rms(error, curves[[overlap_col]]) /
        rmsd(curves[[method_col]], curves[[truth_col]]),
      overlap_relative_amplitude = overlap_relative_amplitude(
        curves[[method_col]], curves[[overlap_col]]
      )
    )
  })
}

# Signal-pattern helpers for Supplement 3
s3_signal_patterns <- c("uniform", "linear_increase", "linear_decline",
                        "exponential_increase", "exponential_decline",
                        "early_peak", "single_peak", "late_peak", "two_peaks",
                        "multi_peak", "plateau")

make_signal_probability <- function(pattern, year_grid) {
  yr <- year_grid

  probability <- switch(
    pattern,

    "uniform" = {
      rep(1, length(yr))
    },

    "linear_increase" = {
      seq(0.2, 1.8, length.out = length(yr))
    },

    "linear_decline" = {
      seq(1.8, 0.2, length.out = length(yr))
    },

    "exponential_increase" = {
      exp(seq(log(0.2), log(2.5), length.out = length(yr)))
    },

    "exponential_decline" = {
      rev(exp(seq(log(0.2), log(2.5), length.out = length(yr))))
    },

    "early_peak" = {
      stats::dnorm(yr, mean = -3650, sd = 220)
    },

    "single_peak" = {
      stats::dnorm(yr, mean = -3000, sd = 280)
    },

    "late_peak" = {
      stats::dnorm(yr, mean = -2200, sd = 220)
    },

    "two_peaks" = {
      0.55 * stats::dnorm(yr, mean = -3600, sd = 220) +
        0.45 * stats::dnorm(yr, mean = -2400, sd = 240)
    },

    "multi_peak" = {
      0.40 * stats::dnorm(yr, mean = -3700, sd = 170) +
        0.35 * stats::dnorm(yr, mean = -2950, sd = 230) +
        0.25 * stats::dnorm(yr, mean = -2200, sd = 160)
    },

    "plateau" = {
      left <- stats::plogis((yr - (-3700)) / 90)
      right <- stats::plogis(((-2300) - yr) / 90)
      left * right
    },

    stop("Unknown signal pattern: ", pattern)
  )

  normalise_curve(probability)
}

simulate_events_from_probability <- function(n, year_grid, probability, seed = NULL) {
  probability <- normalise_curve(probability)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  tibble::tibble(
    id = seq_len(n),
    true_date = sample(year_grid, size = n, replace = TRUE, prob = probability)
  )
}

make_truth_curve_binned <- function(probability, year_grid, bin_width, bin_start) {
  tibble::tibble(year = year_grid, true = normalise_curve(probability)) |>
    make_reconstruction_binned(value_col = true, output_col = "true",
                               bin_width = bin_width, bin_start = bin_start)
}

make_s3_simulation_run <- function(pattern, run_id, phase_catalogue, year_grid,
                                   n_events, bin_width, bin_start, n_mc = 100,
                                   base_seed = 100000) {
  pattern_id <- match(pattern, s3_signal_patterns)

  if (is.na(pattern_id)) {
    stop("Unknown signal pattern: ", pattern)
  }

  probability <- make_signal_probability(
    pattern = pattern, year_grid = year_grid
  )

  events <- simulate_events_from_probability(
    n = n_events,
    year_grid = year_grid,
    probability = probability,
    seed = base_seed + pattern_id * 1000 + run_id
  )

  assigned_events <- assign_events_to_phases(
    events = events, phase_catalogue = phase_catalogue
  )

  true_binned <- make_truth_curve_binned(
    probability = probability,
    year_grid = year_grid,
    bin_width = bin_width,
    bin_start = bin_start
  )

  overlap_curve <- make_phase_catalogue_overlap_curve(
    phase_catalogue = phase_catalogue, year_grid = year_grid
  )

  overlap_binned <- overlap_curve |>
    dplyr::mutate(
      bin = bin_year(year, bin_width = bin_width, bin_start = bin_start)
    ) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(overlap = mean(overlap, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(year = bin)

  standard_curve <- run_standard_aoristic(
    intervals = assigned_events, year_grid = year_grid
  )

  standard_binned <- standard_curve |>
    make_reconstruction_binned(value_col = standard, output_col = "standard",
                               bin_width = bin_width, bin_start = bin_start)

  phase_corrected_curve <- run_phase_corrected_aoristic(
    intervals = assigned_events, year_grid = year_grid
  )

  phase_corrected_binned <- phase_corrected_curve |>
    make_reconstruction_binned(value_col = phase_corrected,
                               output_col = "phase_corrected",
                               bin_width = bin_width, bin_start = bin_start)

  mc_binned <- run_monte_carlo_aoristic(
    intervals = assigned_events,
    year_grid = year_grid,
    bin_width = bin_width,
    bin_start = bin_start,
    n_iter = n_mc,
    seed = base_seed + pattern_id * 1000 + run_id + 500000
  )

  curves_binned <- true_binned |>
    dplyr::left_join(overlap_binned, by = "year") |>
    dplyr::left_join(standard_binned, by = "year") |>
    dplyr::left_join(phase_corrected_binned, by = "year") |>
    dplyr::left_join(
      mc_binned |>
        dplyr::select(year, mc_mean, mc_median, mc_q05, mc_q95),
      by = "year"
    ) |>
    dplyr::mutate(
      true = normalise_curve(tidyr::replace_na(true, 0)),
      standard = normalise_curve(tidyr::replace_na(standard, 0)),
      phase_corrected = normalise_curve(tidyr::replace_na(phase_corrected, 0)),
      mc_mean = normalise_curve(tidyr::replace_na(mc_mean, 0)),
      mc_median = normalise_curve(tidyr::replace_na(mc_median, 0)),
      overlap = tidyr::replace_na(overlap, 0)
    )

  metric_input <- curves_binned |>
    dplyr::select(year, true, overlap, standard, phase_corrected, mc_mean)

  fit_metrics <- compute_curve_metrics(curves = metric_input,
                                       truth_col = "true")

  overlap_diagnostics <- compute_overlap_diagnostics(
    curves = metric_input, truth_col = "true", overlap_col = "overlap"
  )

  metrics <- fit_metrics |>
    dplyr::left_join(overlap_diagnostics, by = "method") |>
    dplyr::mutate(pattern = pattern, run = run_id, .before = 1)

  list(pattern = pattern, run = run_id, probability = probability,
       events = events, assigned_events = assigned_events,
       curves_binned = curves_binned, metrics = metrics)
}

summarise_s3_metrics <- function(metrics) {
  metrics |>
    dplyr::group_by(pattern, method) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_rmsd = mean(rmsd, na.rm = TRUE),
      median_rmsd = stats::median(rmsd, na.rm = TRUE),
      sd_rmsd = stats::sd(rmsd, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      median_mae = stats::median(mae, na.rm = TRUE),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = stats::median(correlation, na.rm = TRUE),
      median_overlap_attributable_rmsd = stats::median(overlap_attributable_rmsd, na.rm = TRUE),
      median_overlap_relative_amplitude = stats::median(overlap_relative_amplitude, na.rm = TRUE),
      mean_overlap_attributable_rmsd = mean(overlap_attributable_rmsd, na.rm = TRUE),
      .groups = "drop"
    )
}

make_s3_best_method_by_run <- function(metrics, metric = "rmsd") {
  metric <- rlang::sym(metric)

  metrics |>
    dplyr::group_by(pattern, run) |>
    dplyr::filter(!!metric == min(!!metric, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::count(pattern, method, name = "best_n") |>
    dplyr::group_by(pattern) |>
    dplyr::mutate(best_share = best_n / sum(best_n)) |>
    dplyr::ungroup()
}

recode_s3_method_names <- function(data) {
  data |>
    dplyr::mutate(
      method = dplyr::recode(method, standard = "Standard aoristic",
                             phase_corrected = "Phase-corrected",
                             mc_mean = "Monte Carlo mean")
    )
}
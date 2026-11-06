# Shared helper functions for Baorista reconstructions
#
# Intervals are half-open (start <= year < end) and are converted to
# Baorista's descending, inclusive representation:
#   -4100 <= year < -3525  ->  4100 ... 3526
# Window -4100 to -1800 with 50-year blocks: timeRange = c(4100, 1801).
safe_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

make_baorista_time_range <- function(time_min, time_max, bin_width) {
  time_min <- as.integer(time_min)
  time_max <- as.integer(time_max)
  bin_width <- as.integer(bin_width)

  older_bp_like <- -time_min
  younger_bp_like_exclusive <- -time_max
  span_years <- older_bp_like - younger_bp_like_exclusive

  if (!is.finite(span_years) || span_years <= 0) {
    stop("Invalid BCE/CE time range after BP-like transformation.")
  }

  if (!is.finite(bin_width) || bin_width <= 0) {
    stop("bin_width must be a positive integer.")
  }

  if (span_years %% bin_width != 0) {
    stop("The requested time span is not divisible by bin_width.\n",
         "time_min = ", time_min, "\n", "time_max = ", time_max, "\n",
         "bin_width = ", bin_width, "\n", "span = ", span_years)
  }

  younger_bp_like_inclusive <- younger_bp_like_exclusive + 1L
  time_range <- c(older_bp_like, younger_bp_like_inclusive)
  n_blocks <- span_years / bin_width

  if (((time_range[1] - time_range[2] + 1L) %% bin_width) != 0) {
    stop("Internal error: Baorista timeRange is not divisible under its ",
         "inclusive block logic.\n",
         "timeRange = ", paste(time_range, collapse = " to "), "\n",
         "resolution = ", bin_width)
  }

  list(time_min = time_min, time_max = time_max, bin_width = bin_width,
       older_bp_like = older_bp_like,
       younger_bp_like_exclusive = younger_bp_like_exclusive,
       younger_bp_like_inclusive = younger_bp_like_inclusive,
       n_blocks = as.integer(n_blocks), time_range = time_range)
}

prepare_baorista_input <- function(intervals, time_range_info) {
  older_limit <- time_range_info$older_bp_like
  younger_limit <- time_range_info$younger_bp_like_inclusive

  out <- intervals |>
    dplyr::transmute(start = -as.numeric(start), end = -as.numeric(end) + 1) |>
    dplyr::mutate(start = pmin(start, older_limit),
                  end = pmax(end, younger_limit)) |>
    dplyr::filter(is.finite(start), is.finite(end), start >= end)

  if (nrow(out) == 0) {
    stop("No valid intervals remain after transformation for Baorista.")
  }

  if (nrow(out) != nrow(intervals)) {
    warning("Retained ", nrow(out), " of ", nrow(intervals),
            " intervals after Baorista transformation.")
  }

  out
}

extract_baorista_block_years <- function(
  prob_mat_object, time_range_info, bin_width = time_range_info$bin_width
) {
  n_model <- ncol(prob_mat_object$pmat)
  years <- seq(from = -time_range_info$older_bp_like, by = bin_width,
               length.out = n_model)

  if (n_model != time_range_info$n_blocks) {
    warning(
      "Baorista returned ", n_model, " time blocks; ",
      time_range_info$n_blocks, " were expected. Using the model-derived axis."
    )
  }

  years
}

extract_normalised_posterior_p <- function(fit, years = NULL) {
  posterior_p <- as.data.frame(fit$posterior.p)
  posterior_p[] <- lapply(posterior_p, as.numeric)

  if (!is.null(years) && ncol(posterior_p) != length(years)) {
    stop("Posterior matrix has ", ncol(posterior_p),
         " columns, but the year axis has ", length(years), ".")
  }

  row_totals <- rowSums(posterior_p, na.rm = TRUE)
  valid <- is.finite(row_totals) & row_totals > 0

  if (!all(valid)) {
    warning("Removing ", sum(!valid),
            " posterior draw(s) with non-finite or non-positive totals.")
    posterior_p <- posterior_p[valid, , drop = FALSE]
    row_totals <- row_totals[valid]
  }

  sweep(as.matrix(posterior_p), MARGIN = 1, STATS = row_totals, FUN = "/")
}

summarise_posterior_p <- function(posterior_p_norm, years) {
  if (ncol(posterior_p_norm) != length(years)) {
    stop("Normalised posterior matrix and year axis have incompatible sizes.")
  }

  purrr::map_dfr(seq_along(years), function(i) {
    samples <- posterior_p_norm[, i]

    hpdi_95 <- tryCatch(
      coda::HPDinterval(coda::mcmc(samples), prob = 0.95),
      error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1)
    )
    hpdi_80 <- tryCatch(
      coda::HPDinterval(coda::mcmc(samples), prob = 0.80),
      error = function(e) matrix(c(NA_real_, NA_real_), nrow = 1)
    )

    tibble::tibble(
      year = years[i],
      baorista_mean = mean(samples, na.rm = TRUE),
      baorista_median = stats::median(samples, na.rm = TRUE),
      baorista_q025 = stats::quantile(samples, 0.025, na.rm = TRUE, names = FALSE),
      baorista_q975 = stats::quantile(samples, 0.975, na.rm = TRUE, names = FALSE),
      baorista_q05 = stats::quantile(samples, 0.05, na.rm = TRUE, names = FALSE),
      baorista_q95 = stats::quantile(samples, 0.95, na.rm = TRUE, names = FALSE),
      baorista_q10 = stats::quantile(samples, 0.10, na.rm = TRUE, names = FALSE),
      baorista_q90 = stats::quantile(samples, 0.90, na.rm = TRUE, names = FALSE),
      baorista_hpdi95_lower = hpdi_95[1, 1],
      baorista_hpdi95_upper = hpdi_95[1, 2],
      baorista_hpdi80_lower = hpdi_80[1, 1],
      baorista_hpdi80_upper = hpdi_80[1, 2]
    )
  })
}

extract_baorista_diagnostics <- function(
  fit, time_range_info, niter, nburnin, thin, nchains, sigma_prior, parallel
) {
  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    fun(x)
  }

  rhat_values <- tryCatch(as.numeric(fit$rhat[[1]][, 1]),
                          error = function(e) NA_real_)
  ess_values <- tryCatch(as.numeric(fit$ess), error = function(e) NA_real_)
  posterior_p <- as.data.frame(fit$posterior.p)

  tibble::tibble(
    diagnostic = c("max_rhat", "median_rhat", "min_ess", "median_ess",
                   "n_posterior_samples", "n_time_blocks", "niter", "nburnin",
                   "thin", "nchains", "timeRange_start", "timeRange_end",
                   "resolution", "inclusive_span", "expected_n_blocks",
                   "sigma_prior", "parallel"),
    value = as.character(c(
      finite_summary(rhat_values, max),
      finite_summary(rhat_values, stats::median),
      finite_summary(ess_values, min),
      finite_summary(ess_values, stats::median),
      nrow(posterior_p),
      ncol(posterior_p),
      niter,
      nburnin,
      thin,
      nchains,
      time_range_info$time_range[1],
      time_range_info$time_range[2],
      time_range_info$bin_width,
      time_range_info$time_range[1] - time_range_info$time_range[2] + 1,
      time_range_info$n_blocks,
      sigma_prior,
      parallel
    ))
  )
}

run_baorista_intervals <- function(intervals, time_min, time_max, bin_width,
                                   niter, nburnin, thin, nchains,
                                   sigma_prior = "dexp(rate=1)",
                                   parallel = TRUE, seeds = seq_len(nchains),
                                   keep_fit = TRUE, rhat_max = 1.05,
                                   max_attempts = 3L) {
  # refit with doubled iterations, burn-in and thinning while max Rhat >
  # rhat_max, up to max_attempts fits

  time_range_info <- make_baorista_time_range(
    time_min = time_min, time_max = time_max, bin_width = bin_width
  )

  baorista_input <- prepare_baorista_input(
    intervals = intervals, time_range_info = time_range_info
  )

  if (length(seeds) != nchains) {
    stop("The number of seeds must equal nchains.")
  }

  prob_mat <- baorista::createProbMat(x = baorista_input,
                                      timeRange = time_range_info$time_range,
                                      resolution = bin_width)

  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    scale <- 2^(attempt - 1L)
    niter_k <- as.integer(niter * scale)
    nburnin_k <- as.integer(nburnin * scale)
    thin_k <- as.integer(thin * scale)

    fit <- baorista::icarfit(x = prob_mat, niter = niter_k, nburnin = nburnin_k,
                             thin = thin_k, nchains = nchains,
                             sigmaPrior = sigma_prior, parallel = parallel,
                             seeds = seeds + (attempt - 1L) * 1000L)

    max_rhat <- suppressWarnings(max(
      as.numeric(tryCatch(fit$rhat[[1]][, 1], error = function(e) NA_real_)),
      na.rm = TRUE
    ))
    converged <- is.finite(max_rhat) && max_rhat <= rhat_max

    message("Baorista attempt ", attempt, ": niter = ", niter_k,
            ", max R-hat = ", signif(max_rhat, 4),
            if (converged) " (criterion met)" else " (criterion not met)")

    if (converged || attempt >= max_attempts) break
  }

  if (!converged) {
    warning("Baorista did not reach max R-hat <= ", rhat_max, " after ",
            attempt, " attempt(s); last max R-hat = ", signif(max_rhat, 4), ".")
  }

  niter <- niter_k
  nburnin <- nburnin_k
  thin <- thin_k

  years <- extract_baorista_block_years(prob_mat_object = prob_mat,
                                        time_range_info = time_range_info,
                                        bin_width = bin_width)

  posterior_draws <- extract_normalised_posterior_p(fit, years)
  posterior_summary <- summarise_posterior_p(posterior_draws, years)
  diagnostics <- extract_baorista_diagnostics(
    fit = fit,
    time_range_info = time_range_info,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    sigma_prior = sigma_prior,
    parallel = parallel
  ) |>
    dplyr::bind_rows(
      tibble::tibble(
        diagnostic = c("rhat_criterion", "n_attempts", "converged"),
        value = as.character(c(rhat_max, attempt, converged))
      )
    )

  result <- list(
    curve = posterior_summary |>
      dplyr::transmute(year, baorista = baorista_median, baorista_mean,
                       baorista_q025, baorista_q975, baorista_q05, baorista_q95,
                       baorista_q10, baorista_q90, baorista_hpdi95_lower,
                       baorista_hpdi95_upper, baorista_hpdi80_lower,
                       baorista_hpdi80_upper),
    posterior_summary = posterior_summary,
    posterior_draws = posterior_draws,
    diagnostics = diagnostics,
    prob_mat = prob_mat,
    time_range_info = time_range_info,
    years = years,
    input = baorista_input,
    settings = list(time_min = time_min, time_max = time_max,
                    bin_width = bin_width, niter = niter, nburnin = nburnin,
                    thin = thin, nchains = nchains, sigma_prior = sigma_prior,
                    parallel = parallel, seeds = seeds)
  )

  if (keep_fit) result$fit <- fit
  result
}

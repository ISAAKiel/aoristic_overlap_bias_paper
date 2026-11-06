#!/usr/bin/env Rscript

# Check that the aoristAAR-based reconstructions are identical to the
# internal reference implementation in R/aoristic_helpers.R.
#
# Run from supplements/ (after S1, S4 and S5 have been rendered once, so
# that the prepared interval files exist):
#
#   Rscript scripts/check_aoristAAR_equivalence.R
#
# Writes outputs/aoristAAR_equivalence_check.csv and stops with an error if
# any difference exceeds the tolerance.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source("R/aoristic_helpers.R")
check_aoristAAR()

tolerance <- 1e-12

datasets <- list(
  "S1 synthetic (uniform)" = list(file = "outputs/s1_assigned_events.csv",
                                  grid = seq(-4100, -1801)),
  "S4 Danish Neolithic" = list(
    file = "outputs/s4_danish_neolithic_prepared_events.csv",
    grid = seq(-4100, -1801)
  ),
  "S5 Hallstatt" = list(file = "outputs/s5_hallstatt_prepared_events.csv",
                        grid = seq(-780, -451))
)

with_backend <- function(backend, expr) {
  old <- Sys.getenv("AORISTIC_BACKEND", NA)
  Sys.setenv(AORISTIC_BACKEND = backend)
  on.exit(if (is.na(old)) Sys.unsetenv("AORISTIC_BACKEND") else Sys.setenv(AORISTIC_BACKEND = old))
  force(expr)
}

results <- purrr::imap_dfr(datasets, function(d, name) {
  if (!file.exists(d$file)) {
    message("Skipping ", name, ": ", d$file, " not found.")
    return(NULL)
  }

  intervals <- read_csv(d$file, show_col_types = FALSE) |>
    select(start, end)

  std_pkg <- with_backend("aoristAAR", run_standard_aoristic(intervals, d$grid))$standard
  std_int <- with_backend("internal", run_standard_aoristic(intervals, d$grid))$standard
  pc_pkg <- with_backend("aoristAAR", run_phase_corrected_aoristic(intervals, d$grid))$phase_corrected
  pc_int <- with_backend("internal", run_phase_corrected_aoristic(intervals, d$grid))$phase_corrected

  tibble(
    dataset = name,
    n_intervals = nrow(intervals),
    method = c("Standard aoristic", "Phase-corrected"),
    max_abs_difference = c(max(abs(std_pkg - std_int)), max(abs(pc_pkg - pc_int))),
    max_relative_difference = c(max(abs(std_pkg - std_int)) / max(std_int),
                                max(abs(pc_pkg - pc_int)) / max(pc_int))
  )
})

results <- results |>
  mutate(aoristAAR_version = as.character(utils::packageVersion("aoristAAR")),
         equivalent = max_abs_difference <= tolerance)

write_csv(results, "outputs/aoristAAR_equivalence_check.csv")
print(results, width = Inf)

if (!all(results$equivalent)) {
  stop("aoristAAR and the internal implementation differ beyond the tolerance (", tolerance, ").")
}

message("aoristAAR and the internal implementation agree within ", tolerance, ".")

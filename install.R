# Install the R packages needed to reproduce the analyses.
# For an exactly pinned environment, prefer renv (see README).

cran_packages <- c(
  "dplyr", "tidyr", "ggplot2", "purrr", "tibble", "readr", "rlang",
  "knitr", "scales", "here", "coda", "rcarbon", "nimble", "baorista", "xronos"
)

missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) try(install.packages(missing))

# Fallback: xronos (R client for the XRONOS database) from GitHub.
if (!requireNamespace("xronos", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("xronos-ch/xronos.R")
}

# aoristAAR (>= 0.1.1) from GitHub.
if (!requireNamespace("aoristAAR", quietly = TRUE) ||
    utils::packageVersion("aoristAAR") < "0.1.1") {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ISAAKiel/aoristAAR@v0.1.1")
}

#!/usr/bin/env bash
set -euo pipefail

for rate in 1 10 50 100; do
  echo "Running Hallstatt Baorista fit with sigma prior rate ${rate}..."

  BAORISTA_OUTPUT_TAG="rate_${rate}" \
  BAORISTA_SIGMA_PRIOR_RATE="${rate}" \
  Rscript scripts/run_s5_baorista.R
done

echo "All Baorista sensitivity fits completed."

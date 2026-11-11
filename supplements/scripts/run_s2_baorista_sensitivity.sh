#!/usr/bin/env bash
set -euo pipefail

# Prior sensitivity of the Supplement 2 Baorista fit (uniform simulation).
# The baseline with rate 1 is computed by run_s2_baorista.R itself.

for rate in 10 50 100 1000; do
  echo "Running S2 Baorista fit with sigma prior rate ${rate}..."

  BAORISTA_OUTPUT_TAG="rate_${rate}" \
  BAORISTA_SIGMA_PRIOR_RATE="${rate}" \
  Rscript scripts/run_s2_baorista.R
done

echo "All S2 Baorista sensitivity fits completed."

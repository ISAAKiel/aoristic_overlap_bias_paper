#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSETS="$ROOT/website/assets"

rm -rf "$ASSETS"
mkdir -p "$ASSETS/figures" "$ASSETS/paper/supplements"

copy_required() {
  local source="$1"
  local target="$2"

  if [[ ! -f "$source" ]]; then
    printf 'Missing required website asset: %s\n' "$source" >&2
    exit 1
  fi

  mkdir -p "$(dirname "$target")"
  cp "$source" "$target"
}

copy_optional_file() {
  local source="$1"
  local target="$2"

  if [[ -f "$source" ]]; then
    mkdir -p "$(dirname "$target")"
    cp "$source" "$target"
  else
    printf 'Optional website asset not found: %s\n' "$source" >&2
  fi
}

copy_optional_dir() {
  local source="$1"
  local target="$2"

  if [[ -d "$source" ]]; then
    rm -rf "$target"
    mkdir -p "$(dirname "$target")"
    cp -R "$source" "$target"
  else
    printf 'Optional website support directory not found: %s\n' "$source" >&2
  fi
}

# -------------------------------------------------------------------------
# Landing-page figures
# -------------------------------------------------------------------------

copy_required \
  "$ROOT/figures/schematic_overlap_bias.png" \
  "$ASSETS/figures/schematic_overlap_bias.png"

copy_required \
  "$ROOT/figures/schematic_phase_correction.png" \
  "$ASSETS/figures/schematic_phase_correction.png"

copy_required \
  "$ROOT/figures/synthetic_phase_catalogue.png" \
  "$ASSETS/figures/synthetic_phase_catalogue.png"

copy_required \
  "$ROOT/supplements/figures/s1_overlap_structure.png" \
  "$ASSETS/figures/s1_overlap_structure.png"

copy_required \
  "$ROOT/supplements/figures/s1_uniform_ground_truth.png" \
  "$ASSETS/figures/s1_uniform_ground_truth.png"

copy_required \
  "$ROOT/supplements/figures/s1_ground_truth_vs_standard_aoristic.png" \
  "$ASSETS/figures/s1_ground_truth_vs_standard_aoristic.png"

copy_required \
  "$ROOT/supplements/figures/s2_reconstruction_methods_with_uncertainty.png" \
  "$ASSETS/figures/s2_reconstruction_methods_with_uncertainty.png"

copy_required \
  "$ROOT/supplements/figures/s3_ground_truth_signal_shapes.png" \
  "$ASSETS/figures/s3_ground_truth_signal_shapes.png"

copy_required \
  "$ROOT/supplements/figures/s3_rmsd_boxplot.png" \
  "$ASSETS/figures/s3_rmsd_boxplot.png"

copy_required \
  "$ROOT/supplements/figures/s4_danish_neolithic_methods_vs_overlap.png" \
  "$ASSETS/figures/s4_danish_neolithic_methods_vs_overlap.png"

copy_required \
  "$ROOT/supplements/figures/s5_hallstatt_methods_vs_overlap.png" \
  "$ASSETS/figures/s5_hallstatt_methods_vs_overlap.png"

# -------------------------------------------------------------------------
# Manuscript (HTML needs its *_files directory)
# -------------------------------------------------------------------------

copy_optional_file \
  "$ROOT/manuscript.html" \
  "$ASSETS/paper/manuscript.html"

copy_optional_file \
  "$ROOT/manuscript.pdf" \
  "$ASSETS/paper/manuscript.pdf"

copy_optional_dir \
  "$ROOT/manuscript_files" \
  "$ASSETS/paper/manuscript_files"

# The manuscript HTML refers to figures/ and supplements/figures/ by relative
# paths, so copy every image it uses next to it.
if [[ -f "$ROOT/manuscript.html" ]]; then
  grep -o '<img src="[^"]*"' "$ROOT/manuscript.html" \
    | sed -e 's/^<img src="//' -e 's/"$//' \
    | grep -v -e '^data:' -e '^https\?:' -e '^manuscript_files/' \
    | sort -u \
    | while read -r img; do
        copy_required "$ROOT/$img" "$ASSETS/paper/$img"
      done
fi

# -------------------------------------------------------------------------
# Supplements and their Quarto support directories
# -------------------------------------------------------------------------

for supplement in \
  s1_structural_phase_bias \
  s2_reconstruction_method_comparison \
  s3_non_uniform_signal_sensitivity \
  s4_empirical_1_danish_neolithic \
  s5_empirical_2_hallstatt \
  s6_assignment_model_sensitivity
do
  copy_optional_file \
    "$ROOT/supplements/${supplement}.html" \
    "$ASSETS/paper/supplements/${supplement}.html"

  copy_optional_file \
    "$ROOT/supplements/${supplement}.pdf" \
    "$ASSETS/paper/supplements/${supplement}.pdf"

  copy_optional_dir \
    "$ROOT/supplements/${supplement}_files" \
    "$ASSETS/paper/supplements/${supplement}_files"
done

printf 'Website assets prepared in %s\n' "$ASSETS"

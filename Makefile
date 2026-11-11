SHELL := /bin/bash

QUARTO ?= quarto
RSCRIPT ?= Rscript

MANUSCRIPT := manuscript.qmd
SUPPLEMENTS := \
	supplements/s1_structural_phase_bias.qmd \
	supplements/s2_reconstruction_method_comparison.qmd \
	supplements/s3_non_uniform_signal_sensitivity.qmd \
	supplements/s4_empirical_1_danish_neolithic.qmd \
	supplements/s5_empirical_2_hallstatt.qmd \
	supplements/s6_assignment_model_sensitivity.qmd

.PHONY: help site site-assets preview-site manuscript supplements paper all \
        clean-site clean-rendered clean data analysis \
        sensitivity s6 s6-danish s6-hallstatt check-aoristaar refresh-s2-baorista refresh-s2-baorista-sensitivity refresh-s3-baorista \
        refresh-s4-baorista refresh-s5-baorista refresh-expensive check-tools

# All analysis scripts expect to be run from supplements/.
SUPP := cd supplements &&

help:
	@printf '%s\n' \
	  'Available targets:' \
	  '  make site              Stage assets and render the website to docs/' \
	  '  make preview-site      Stage assets and preview the website locally' \
	  '  make manuscript        Render manuscript.qmd' \
	  '  make supplements       Render all supplementary QMD files' \
	  '  make paper             Render supplements, then the manuscript' \
	  '  make data              Check/rebuild the Danish settlement subset' \
	  '  make check-aoristaar   Compare aoristAAR with the internal implementation' \
	  '  make analysis          Rerun the S3 and S6 simulations (minutes)' \
	  '  make sensitivity       Rerun the S3 non-uniform simulation only' \
	  '  make s6                Rerun the S6 assignment-model simulation only' \
	  '  make s6-danish         Rerun the S6 simulation on the Danish catalogue (needs S4)' \
	  '  make s6-hallstatt      Rerun the S6 simulation on the Hallstatt catalogue (needs S5)' \
	  '  make refresh-s2-baorista-sensitivity  S2 Baorista with prior rates 10, 50, 100, 1000' \
	  '  make refresh-expensive Rerun all Baorista fits (S2–S5; hours)' \
	  '  make all               Render paper and website' \
	  '  make clean-site        Remove docs/ and staged website assets' \
	  '  make clean-rendered    Remove Quarto-generated HTML support directories'

check-tools:
	@command -v "$(QUARTO)" >/dev/null || { echo "Quarto not found: $(QUARTO)"; exit 1; }
	@command -v "$(RSCRIPT)" >/dev/null || { echo "Rscript not found: $(RSCRIPT)"; exit 1; }

site-assets:
	bash website/prepare_assets.sh

site: site-assets
	$(QUARTO) render website

preview-site: site-assets
	$(QUARTO) preview website

manuscript: check-tools
	$(QUARTO) render $(MANUSCRIPT)

supplements: check-tools
	@set -e; \
	for file in $(SUPPLEMENTS); do \
	  echo "Rendering $$file"; \
	  $(QUARTO) render "$$file"; \
	done

# The manuscript reads numbers from supplements/outputs, so supplements first.
paper: supplements manuscript

data:
	$(SUPP) $(RSCRIPT) scripts/prepare_s4_danish_settlements.R

# Requires rendered S1, S4 and S5 (prepared interval files).
check-aoristaar:
	$(SUPP) $(RSCRIPT) scripts/check_aoristAAR_equivalence.R

sensitivity:
	$(SUPP) $(RSCRIPT) scripts/run_s3_non_uniform_signal_sensitivity.R

s6:
	$(SUPP) $(RSCRIPT) scripts/run_s6_assignment_model_sensitivity.R

# Requires outputs/s4_danish_neolithic_interval_counts.csv: render S4 first.
s6-danish:
	$(SUPP) $(RSCRIPT) scripts/run_s6_danish_catalogue.R

# Requires outputs/s5_hallstatt_interval_counts.csv: render S5 first.
s6-hallstatt:
	$(SUPP) $(RSCRIPT) scripts/run_s6_hallstatt_catalogue.R

analysis: sensitivity s6

refresh-s2-baorista:
	$(SUPP) REFRESH_BAORISTA=true $(RSCRIPT) scripts/run_s2_baorista.R

# Prior sensitivity of the S2 fit (rates 10, 50, 100, 1000)
refresh-s2-baorista-sensitivity:
	$(SUPP) bash scripts/run_s2_baorista_sensitivity.sh

refresh-s3-baorista:
	$(SUPP) REFRESH_BAORISTA=true $(RSCRIPT) scripts/run_s3_baorista_examples.R

# Requires outputs/s4_danish_neolithic_prepared_events.csv: render S4 first.
refresh-s4-baorista:
	$(SUPP) bash scripts/run_s4_baorista_sensitivity.sh

# Requires outputs/s5_hallstatt_prepared_events.csv: render S5 first.
refresh-s5-baorista:
	$(SUPP) bash scripts/run_s5_baorista_sensitivity.sh

refresh-expensive: refresh-s2-baorista refresh-s2-baorista-sensitivity refresh-s3-baorista refresh-s4-baorista refresh-s5-baorista

all: paper site

clean-site:
	rm -rf docs website/assets

clean-rendered:
	rm -rf manuscript.html manuscript_files
	rm -rf supplements/*_files
	rm -f supplements/*.html supplements/*.pdf

clean: clean-site clean-rendered
	@echo "Source data, generated CSV outputs, and cached Bayesian fits were preserved."

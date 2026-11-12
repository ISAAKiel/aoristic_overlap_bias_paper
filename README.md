# Overlap bias in aoristic reconstructions

Code, data and supplements for

> Hinz, M. & Müller-Scheeßel, N. Systematic bias from overlapping temporal classifications requires correction in aoristic reconstructions. To be submitted to the *Journal of Archaeological Science*.

Archived at Zenodo: doi:… (TODO)

## What the paper is about

Aoristic analysis spreads each phase-dated observation evenly over its interval and adds everything up. That works badly when the phases overlap, as they do in most real chronologies (a site dated to "Early Neolithic" and another to "EN II" end up in the same sum). The standard aoristic curve then mostly reproduces how many labels cover each year, even if the underlying events are spread perfectly evenly. Monte Carlo aoristics and the Bayesian approach of Baorista don't fix this.

We derive a simple correction from the assumption that an event can receive any of the labels that cover its date with equal probability. It needs nothing beyond the intervals you already have. We test it in simulations and apply it to Danish Neolithic settlements and to Early Iron Age graves from southern Germany.

## Reading the paper

If you just want to read the paper and the supplements, you don't need to run anything. The PDFs are in the repository:

- `manuscript.pdf`
- `supplements/s1_*.pdf` to `supplements/s6_*.pdf`

S1–S3 are the simulations, S4 is the Danish Neolithic case study, S5 the Hallstatt burials and S6 looks at what happens if labels are not assigned the way the correction assumes.

The correction itself is implemented in the R package [aoristAAR](https://github.com/ISAAKiel/aoristAAR).

## Rebuilding the documents

You need R (we used 4.5.3), Quarto and a LaTeX installation. Get the package versions we used with

```
Rscript -e 'renv::restore()'
```

or, if you don't care about exact versions, `Rscript install.R`. Baorista needs nimble, and nimble needs a C++ compiler.

All results of the expensive steps are cached in `supplements/outputs/`, so

```
make paper
```

renders all six supplements and the manuscript in a few minutes. `make site` rebuilds the website in `docs/`.

## Running everything from scratch

This takes several hours, mostly because of the Baorista fits. Order matters, since some steps read the output of earlier ones:

```
make data                                                       # Danish subset from the full register export
quarto render supplements/s1_structural_phase_bias.qmd          # writes the synthetic catalogue
make analysis                                                   # simulations for S3 and S6
quarto render supplements/s4_empirical_1_danish_neolithic.qmd   # writes the event table for Baorista
make refresh-expensive                                          # Baorista for S2-S5 (hours)
make s6-danish s6-hallstatt                                     # S6 runs on the empirical catalogues
make paper
```

`make help` lists all targets. Scripts are meant to be run from `supplements/`, which the Makefile does for you.

Some things worth knowing:

- Seeds are fixed. `N_RUNS` and `N_EVENTS` change the size of the simulations, and each run writes its settings to `outputs/*_run_metadata.csv`.
- Baorista uses four chains and 50,000 iterations. If R-hat is above 1.05 the fit is repeated with twice the iterations, at most three times. Not all Danish fits converge, see S4.
- Intervals are half-open, `[start, end)`. Inclusive end dates in the source data are shifted accordingly.
- In the Danish data a few dates were entered freely (radiocarbon ranges and the like) instead of picked from the period list. These are mapped to the closest period. Set `HARMONISE_FREE_ENTRIES=false` to switch this off.
- `make check-aoristaar` checks that aoristAAR gives the same results as our own implementation in `supplements/R/aoristic_helpers.R`. `AORISTIC_BACKEND=internal` uses the latter.

## Data

- `data/anlaeg_all.csv`: full export of Fund og Fortidsminder (Slots- og Kulturstyrelsen), downloaded 5 April 2018, Latin-1. CC0 according to the provider, attribution requested. `data/anlaeg_settlements.csv` is the Neolithic settlement subset derived from it (`supplements/scripts/prepare_s4_danish_settlements.R`).
- `data/xronos_denmark_*`: Danish radiocarbon dates from XRONOS (CC BY). `REFRESH_XRONOS=true` downloads them again.
- `data/hallstatt.tab` and `data/hallstatt_phase_ranges.csv`: graves and phase ranges from Müller-Scheeßel 2013, also available via [Tübingen](http://nbn-resolving.de/urn:nbn:de:bsz:21-opus-70996) and [Zenodo](https://zenodo.org/records/15012694). CC BY 4.0.

## Licence

Code is MIT (`LICENSE`), text, figures and derived data are CC BY 4.0 (`LICENSE-CC-BY.md`). The third-party data keep their own licences, see above.

## Citing

Please cite the article. If you use code or data from here, also cite the Zenodo archive (see `CITATION.cff`, or "Cite this repository" on GitHub).

Martin Hinz's work on this was funded by the ERC project ESTER (101171036), Nils Müller-Scheeßel's by the Cluster of Excellence ROOTS (EXC 2150, 390870439).

Questions and bug reports: open an issue, or write to martin.hinz@ufg.uni-kiel.de.

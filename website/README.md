# Project website

Quarto source of the project website, published with GitHub Pages from `docs/`
at https://isaakiel.github.io/aoristic_overlap_bias_paper/. The page presents the main
argument and links to the manuscript, the supplements and the repository.

## Render

```bash
make site
```

This first runs `website/prepare_assets.sh`, which copies the selected figures,
manuscript, and supplement HTML files into `website/assets/`. Quarto then copies
those assets into `docs/`, producing stable links on GitHub Pages.

For a local preview:

```bash
make preview-site
```

## Rendering order

To ensure that the linked manuscript and supplements are current:

```bash
make all
```

This renders the manuscript and supplements first, stages the resulting files,
and then renders the website.

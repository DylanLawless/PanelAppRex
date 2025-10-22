#!/bin/bash
set -e

OUT=panelapprex2025lawless_appnote_source_v1.tar.gz

tar -czf "$OUT" -C latex \
  panelapprex2025lawless_appnote.tex \
  references.bib \
  head.tex \
  images/figure_1.pdf \
  images/validation_counts.pdf \
  images/plot_patch2_annotated_example.pdf \
  images/plot_patch_uniq_gene_moi_summary.pdf

echo "Created $OUT"


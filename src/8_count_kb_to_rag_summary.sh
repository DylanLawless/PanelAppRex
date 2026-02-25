#!/usr/bin/env bash

summary_dir=~/so/src/PanelAppRex/data
summary_file="$summary_dir/kb_to_rag_summary.tsv"

# compare_rag_sizes.sh
# Usage:
#   ./8_count_kb_to_rag_summary.sh [RAW_DIR] [RAG_DIR]
# Defaults:
#   RAW_DIR = uniprot_kb_up000005640_processed
#   RAG_DIR = uniprot_kb_up000005640_processed_rag

cd ~/mnt/atlas_data_big/data/panelapprex/ || {
  echo "Failed to cd into ~/mnt/atlas_data_big/data/panelapprex/" >&2
  exit 1
}

RAW_DIR="${1:-uniprot_kb_up000005640_processed}"
RAG_DIR="${2:-uniprot_kb_up000005640_processed_rag}"

if [ ! -d "$RAW_DIR" ]; then
  echo "Raw directory not found: $RAW_DIR" >&2
  exit 1
fi

if [ ! -d "$RAG_DIR" ]; then
  echo "RAG directory not found: $RAG_DIR" >&2
  exit 1
fi

echo "Raw input (UniProt parsed TSV) directory: $RAW_DIR"
echo "RAG summaries (panel info boxes) directory: $RAG_DIR"
echo

# Raw side: single gene level TSV
raw_pattern="$RAW_DIR/uniprot_gene_level.tsv"
if [ ! -f "$raw_pattern" ]; then
  echo "Raw file not found: $raw_pattern" >&2
  exit 1
fi

echo "Raw file pattern:"
echo "  $raw_pattern"
echo "Raw file:"
ls "$raw_pattern" | sed 's/^/  /'
echo

raw_total_line=$(wc "$raw_pattern")
read -r raw_lines raw_words raw_bytes _ <<< "$raw_total_line"

echo "Raw UniProt derived text (gene level TSV):"
echo "  lines: $raw_lines"
echo "  words: $raw_words"
echo "  bytes: $raw_bytes"
echo

# RAG side: full panel summaries
rag_pattern="$RAG_DIR"/*full*txt
if ! ls $rag_pattern >/dev/null 2>&1; then
  echo "No *full*txt files in $RAG_DIR" >&2
  exit 1
fi

echo "RAG file pattern:"
echo "  $rag_pattern"
echo "Example RAG files (first 3):"
ls $rag_pattern | head -n 3 | sed 's/^/  /'
echo

rag_total_line=$(wc $rag_pattern | tail -n 1)
read -r rag_lines rag_words rag_bytes _ <<< "$rag_total_line"

echo "AI RAG panel summaries (full text files):"
echo "  lines: $rag_lines"
echo "  words: $rag_words"
echo "  bytes: $rag_bytes"
echo

# Ratios
compression_ratio_bytes=$(awk -v a="$raw_bytes" -v b="$rag_bytes" 'BEGIN { if (a == 0) {print "NA"} else {printf "%.6f", b / a} }')
reduction_factor_bytes=$(awk -v a="$raw_bytes" -v b="$rag_bytes" 'BEGIN { if (b == 0) {print "NA"} else {printf "%.2f", a / b} }')

compression_ratio_words=$(awk -v a="$raw_words" -v b="$rag_words" 'BEGIN { if (a == 0) {print "NA"} else {printf "%.6f", b / a} }')
reduction_factor_words=$(awk -v a="$raw_words" -v b="$rag_words" 'BEGIN { if (b == 0) {print "NA"} else {printf "%.2f", a / b} }')

echo "Compression summary (RAG vs raw):"
echo "  bytes: RAG / raw = $compression_ratio_bytes   (raw is ${reduction_factor_bytes}x larger)"
echo "  words: RAG / raw = $compression_ratio_words   (raw is ${reduction_factor_words}x larger)"
echo

# Save summary table

if [ ! -f "$summary_file" ]; then
  echo -e "timestamp\traw_dir\trag_dir\traw_file\trag_pattern\traw_lines\traw_words\traw_bytes\trag_lines\trag_words\trag_bytes\tcompression_ratio_bytes\tcompression_ratio_words" > "$summary_file"
fi

timestamp=$(date -Iseconds)
raw_file_name=$(basename "$raw_pattern")

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$timestamp" \
  "$RAW_DIR" \
  "$RAG_DIR" \
  "$raw_file_name" \
  "$rag_pattern" \
  "$raw_lines" \
  "$raw_words" \
  "$raw_bytes" \
  "$rag_lines" \
  "$rag_words" \
  "$rag_bytes" \
  "$compression_ratio_bytes" \
  "$compression_ratio_words" >> "$summary_file"

echo "Summary appended to: $summary_file"

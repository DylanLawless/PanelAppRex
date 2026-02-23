#!/usr/bin/env bash
set -euo pipefail

# UniProtKB reference proteome for human (UP000005640, taxon 9606)
# Source:
#   https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Eukaryota/UP000005640/
#
# This script downloads:
#   - RELEASE.metalink                   meta information and mirror URLs for this proteome release
#   - UP000005640_9606.dat.gz            UniProtKB flat file with full reviewed annotations
#   - UP000005640_9606.xml.gz            UniProtKB XML with full reviewed annotations
#   - UP000005640_9606.fasta.gz          protein sequences
#   - UP000005640_9606_DNA.fasta.gz      coding DNA sequences
#   - UP000005640_9606.gene2acc.gz       mapping between genes and accessions
#   - UP000005640_9606.idmapping.gz      cross reference mappings
#   - UP000005640_9606.score.gz          proteome quality scores
#   - UP000005640_9606_additional.*      additional entries associated with the proteome
#
# Target directory (adjust if needed):
#   $HOME/mnt/atlas_data_big/data/panelapprex/uniprot_up000005640

BASE_DIR="$HOME/mnt/atlas_data_big/data/panelapprex/uniprot_up000005640"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR"

BASE_URL="https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/reference_proteomes/Eukaryota/UP000005640"

FILES=(
  "RELEASE.metalink"
  "UP000005640_9606.dat.gz"
  "UP000005640_9606.fasta.gz"
  "UP000005640_9606.gene2acc.gz"
  "UP000005640_9606.idmapping.gz"
  "UP000005640_9606.score.gz"
  "UP000005640_9606.xml.gz"
  "UP000005640_9606_DNA.fasta.gz"
  "UP000005640_9606_additional.dat.gz"
  "UP000005640_9606_additional.fasta.gz"
  "UP000005640_9606_additional.fasta_NORMAL.gz"
  "UP000005640_9606_additional.score.gz"
  "UP000005640_9606_additional.xml.gz"
)

echo "Downloading UniProt human reference proteome files to:"
echo "  $BASE_DIR"
echo

for f in "${FILES[@]}"; do
  echo "→ $f"
  curl -fSL "${BASE_URL}/${f}" -o "${f}"
done

echo
echo "All downloads completed."

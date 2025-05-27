#!/bin/bash

# Exit on errors
set -e

# ---------- Parse arguments ----------
if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <counts_matrix.csv> <metadata.csv> <design_formula> <extra_info.txt>"
  exit 1
fi

COUNTS_MATRIX="$1"
METADATA="$2"
DESIGN="$3"
EXTRA_INFO_FILE="$4"

# Optional args with defaults
SPECIES="${5:-human}"
TXI_STATUS="${6:-TRUE}"
ANNOTATION="${7:-biomart}"
TREATMENT="${8:-untreated}"
DO_PATHWAY="${9:-TRUE}"

# Read extra info from file (supports multiline)
EXTRA_INFO=$(<"$EXTRA_INFO_FILE")

# R script
R_SCRIPT="run_pipeline.R"

# Output folder
OUTPUT_DIR="./results"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Run R script
echo "Running R analysis..."
Rscript "../$R_SCRIPT" "$COUNTS_MATRIX" "$METADATA" "$DESIGN" "$SPECIES" "$TXI_STATUS" "$ANNOTATION" "$TREATMENT" "$DO_PATHWAY" "$EXTRA_INFO"

# Read the title from R output
TITLE=$(<title.txt)
SAFE_TITLE=$(echo "$TITLE" | tr ' /:' '_' | tr -cd '[:alnum:]_-')
PDF_OUTPUT="${SAFE_TITLE}.pdf"

# Convert all PNGs to PDF
echo "Combining PNGs into a PDF..."
convert *.png "$PDF_OUTPUT"

echo "✅ Analysis complete. PDF saved as: $OUTPUT_DIR/$PDF_OUTPUT"

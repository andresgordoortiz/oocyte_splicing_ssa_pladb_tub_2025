#!/usr/bin/env bash
set -euo pipefail

# Usage check
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <chromosome> <start> <end>"
  echo "Example: $0 16 32677894 32680807"
  exit 1
fi

CHROM=$1
START=$2
END=$3

INPUT_GTF="gencode.vM10.annotation.gtf"
OUTPUT_TSV="features_chr${CHROM}_${START}_${END}.exons_introns.lengths.tsv"

awk -F $'\t' \
  -v chr="chr${CHROM}" -v lo="$START" -v hi="$END" \
  'BEGIN {
     # Print header once
     print "chrom\tfeature\tstart\tend\tlen\tattribute"
   }
   $1 == chr        &&
   $4 >= lo         &&
   $5 <= hi         &&
   ($3 == "exon" || $3 == "intron") {
     len = $5 - $4;
     print $1 "\t" $3 "\t" $4 "\t" $5 "\t" len "\t" $9;
   }
  ' "$INPUT_GTF" \
| column -t -s $'\t' \
> "$OUTPUT_TSV"

echo "Wrote: $OUTPUT_TSV"

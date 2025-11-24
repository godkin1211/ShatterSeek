#!/bin/bash
#
# CORRECTED script for analyzing chromothripsis from TSV output
# This version uses CORRECT fragment joins test logic: p > 0.05 (not < 0.05)
#
# Usage: ./analyze_chromothripsis_CORRECTED.sh chromothripsis.tsv
#

if [ $# -eq 0 ]; then
  echo "Usage: $0 <chromothripsis.tsv>"
  exit 1
fi

INPUT_FILE=$1

if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: File $INPUT_FILE not found"
  exit 1
fi

echo "================================================================================"
echo "  CHROMOTHRIPSIS CLASSIFICATION - CORRECTED VERSION"
echo "================================================================================"
echo ""
echo "Input file: $INPUT_FILE"
echo ""
echo "Classification criteria (CORRECT LOGIC):"
echo "  High confidence Type 1:"
echo "    - Intrachromosomal SVs ≥ 6"
echo "    - CN oscillations (2-state) ≥ 7"
echo "    - Fragment joins p-value > 0.05 (CORRECTED: random joining)"
echo "    - At least one clustering test p < 0.05"
echo ""
echo "  Low confidence:"
echo "    - Intrachromosomal SVs ≥ 6"
echo "    - CN oscillations (2-state) 4-6"
echo "    - Fragment joins p-value > 0.05 (CORRECTED: random joining)"
echo "    - At least one clustering test p < 0.05"
echo ""
echo "================================================================================"
echo ""

# Header
echo "Results:"
echo ""
printf "%-8s %-8s %-8s %-8s %-12s %-12s %-12s %-20s\n" \
  "Chrom" "IntraSVs" "TotalSVs" "CN_Osc" "Frag_Joins" "Chr_Enrich" "Exp_Chr" "Classification"
echo "--------------------------------------------------------------------------------"

# Read file line by line (skip header)
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r chrom clusterSize clusterSize_TRA number_TRA \
  number_DEL number_h2hINV number_t2tINV number_DUP number_INV number_translocation \
  pval_fragment_joins inter_pval_fragment_joins \
  chr_breakpoint_enrichment pval_exp_chr pval_exp_cluster \
  max_number_oscillating_CN_segments_2_states max_number_oscillating_CN_segments_3_states; do

  # Skip if no cluster
  if [[ "$clusterSize" -eq 0 ]]; then
    continue
  fi

  # Parse numeric values (handle scientific notation)
  numIntra=$clusterSize
  numTotal=$clusterSize_TRA
  cnOsc=$max_number_oscillating_CN_segments_2_states

  # Convert p-values to comparable format
  pval_frag=$(echo "$pval_fragment_joins" | awk '{printf "%.5f", $1}')
  pval_chr_enrich=$(echo "$chr_breakpoint_enrichment" | awk '{printf "%.2e", $1}')
  pval_exp=$(echo "$pval_exp_chr" | awk '{printf "%.2e", $1}')

  classification="Not chromothripsis"

  # CORRECTED High confidence Type 1
  # Fragment joins p > 0.05 (CORRECTED from < to >)
  if [[ "$numIntra" -ge 6 ]] && [[ "$cnOsc" -ge 7 ]] &&
    awk -v p="$pval_fragment_joins" 'BEGIN {exit !(p > 0.05)}' &&
    (awk -v p="$chr_breakpoint_enrichment" 'BEGIN {exit !(p < 0.05)}' ||
      awk -v p="$pval_exp_chr" 'BEGIN {exit !(p < 0.05)}'); then
    classification="High confidence"
  # CORRECTED Low confidence
  # Fragment joins p > 0.05 (CORRECTED from < to >)
  elif [[ "$numIntra" -ge 6 ]] && [ "$cnOsc" -ge 4 ] && [ "$cnOsc" -le 6 ] &&
    awk -v p="$pval_fragment_joins" 'BEGIN {exit !(p > 0.05)}' &&
    (awk -v p="$chr_breakpoint_enrichment" 'BEGIN {exit !(p < 0.05)}' ||
      awk -v p="$pval_exp_chr" 'BEGIN {exit !(p < 0.05)}'); then
    classification="Low confidence"
  fi

  # Only print if classified as chromothripsis
  if [[ "$classification" != "Not chromothripsis" ]]; then
    printf "%-8s %-8s %-8s %-8s %-12s %-12s %-12s %-20s\n" \
      "$chrom" "$numIntra" "$numTotal" "$cnOsc" \
      "$pval_frag" "$pval_chr_enrich" "$pval_exp" "$classification"
  fi
done

echo ""
echo "================================================================================"
echo ""
echo "Summary:"
echo ""

# Count results
high_count=$(tail -n +2 "$INPUT_FILE" | awk -F'\t' '
    $2 > 0 && $16 >= 7 && $2 >= 6 && $11 > 0.05 && ($13 < 0.05 || $14 < 0.05)
' | wc -l)

low_count=$(tail -n +2 "$INPUT_FILE" | awk -F'\t' '
    $2 > 0 && $16 >= 4 && $16 <= 6 && $2 >= 6 && $11 > 0.05 && ($13 < 0.05 || $14 < 0.05)
' | wc -l)

echo "  High confidence chromothripsis: $high_count"
echo "  Low confidence chromothripsis:  $low_count"
echo ""
echo "================================================================================"

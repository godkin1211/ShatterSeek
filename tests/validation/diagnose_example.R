#!/usr/bin/env Rscript
#
# Diagnostic analysis of example data to understand classification changes
#

cat("\n")
cat("========================================================================\n")
cat(" DIAGNOSTIC ANALYSIS: Why are detections missing?\n")
cat("========================================================================\n\n")

# Load package
suppressPackageStartupMessages({
    library(ShatterSeek)
})

# Load example data
data(DO17373)

cat("Step 1: Preparing data...\n")

# Prepare SV data
SV_data <- SVs(
    chrom1 = as.character(SV_DO17373$chrom1),
    pos1 = as.numeric(SV_DO17373$start1),
    chrom2 = as.character(SV_DO17373$chrom2),
    pos2 = as.numeric(SV_DO17373$end2),
    SVtype = as.character(SV_DO17373$svclass),
    strand1 = as.character(SV_DO17373$strand1),
    strand2 = as.character(SV_DO17373$strand2)
)

# Prepare CNV data
CN_data <- CNVsegs(
    chrom = as.character(SCNA_DO17373$chromosome),
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

cat("  ✓ Data prepared\n\n")

# Run detection
cat("Step 2: Running chromothripsis detection...\n")

chromoth_result <- detect_chromothripsis(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg19"
)

cat("  ✓ Detection complete\n\n")

# Classify events
cat("Step 3: Classifying chromothripsis events...\n")

classification <- classify_chromothripsis(chromoth_result)

cat("  ✓ Classification complete\n\n")

# Show quick summary
cat("========================================================================\n")
cat(" QUICK SUMMARY\n")
cat("========================================================================\n\n")

cat(sprintf("Total chromosomes analyzed: %d\n", nrow(classification)))
cat(sprintf("  - High confidence: %d\n",
           sum(classification$classification == "High confidence")))
cat(sprintf("  - Low confidence:  %d\n",
           sum(classification$classification == "Low confidence")))
cat(sprintf("  - Not chromothripsis: %d\n\n",
           sum(classification$classification == "Not chromothripsis")))

# Now run diagnostic analysis
cat("========================================================================\n")
cat(" DETAILED DIAGNOSTIC ANALYSIS\n")
cat("========================================================================\n\n")

cat("Running diagnostic on all chromosomes with clusters...\n\n")

# Run diagnostic
diagnostic <- diagnose_chromothripsis_classification(chromoth_result,
                                                     show_all_chromosomes = FALSE)

# Print detailed report
print_diagnostic_report(diagnostic)

# Additional analysis: Show what would happen with different thresholds
cat("\n")
cat("========================================================================\n")
cat(" THRESHOLD SENSITIVITY ANALYSIS\n")
cat("========================================================================\n\n")

cat("Chromosomes that ALMOST meet criteria:\n\n")

chromSummary <- chromoth_result@chromSummary
chromSummary <- chromSummary[chromSummary$clusterSize > 0, ]

for (i in 1:nrow(chromSummary)) {
    row <- chromSummary[i, ]

    # Get CN oscillations
    cn_osc <- row$max_number_oscillating_CN_segments_2_states
    if (is.na(cn_osc)) {
        cn_osc <- row$max_number_oscillating_CN_segments_3_states
    }

    # Check fragment joins
    frag_joins_ok <- !is.na(row$pval_fragment_joins) && row$pval_fragment_joins > 0.05

    # Check clustering
    clustering_ok <- (!is.na(row$pval_exp_cluster) && row$pval_exp_cluster < 0.05) ||
                     (!is.na(row$chr_breakpoint_enrichment) && row$chr_breakpoint_enrichment < 0.05)

    # Identify "almost" cases
    if (row$clusterSize >= 4 && row$clusterSize < 6) {
        cat(sprintf("Chr %s: %d intrachromosomal SVs (need 6, have %d - CLOSE!)\n",
                   row$chrom, row$clusterSize, row$clusterSize))
    }

    if (!is.na(cn_osc) && cn_osc >= 3 && cn_osc < 4) {
        cat(sprintf("Chr %s: %d CN oscillations (need 4, have %d - CLOSE!)\n",
                   row$chrom, cn_osc, cn_osc))
    }

    if (!is.na(row$pval_fragment_joins) && row$pval_fragment_joins > 0.03 &&
        row$pval_fragment_joins <= 0.05) {
        cat(sprintf("Chr %s: Fragment joins p=%.3f (need >0.05, BORDERLINE!)\n",
                   row$chrom, row$pval_fragment_joins))
    }
}

cat("\n")
cat("========================================================================\n\n")

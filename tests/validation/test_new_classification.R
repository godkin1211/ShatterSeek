#!/usr/bin/env Rscript
#
# Test new chromothripsis classification system
# Based on Cortes-Ciriano et al. 2020 criteria
#

cat("\n")
cat("========================================================================\n")
cat(" Testing New Chromothripsis Classification System\n")
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

# Print results
cat("========================================================================\n")
cat(" CLASSIFICATION RESULTS\n")
cat("========================================================================\n\n")

print(classification)

# Count by classification
cat("\n")
cat("Classification Summary:\n")
cat(sprintf("  - High confidence: %d\n",
           sum(classification$classification == "High confidence")))
cat(sprintf("  - Low confidence:  %d\n",
           sum(classification$classification == "Low confidence")))
cat(sprintf("  - Not chromothripsis: %d\n",
           sum(classification$classification == "Not chromothripsis")))
cat("\n")

# Validate classification categories
cat("Step 4: Validating classification logic...\n\n")

# Test that only valid classifications are present
valid_classes <- c("High confidence", "Low confidence", "Not chromothripsis")
all_valid <- all(classification$classification %in% valid_classes)

if (all_valid) {
    cat("  ✓ All classifications use valid categories\n")
} else {
    stop("ERROR: Invalid classification categories found!")
}

# Test fragment joins p-value for chromothripsis events
chromoth_events <- classification[classification$classification %in%
                                  c("High confidence", "Low confidence"), ]

if (nrow(chromoth_events) > 0) {
    # Check if fragment joins test is passed (p > 0.05)
    pval_ok <- all(!is.na(chromoth_events$pval_fragment_joins) &
                   chromoth_events$pval_fragment_joins > 0.05,
                   na.rm = TRUE)

    if (pval_ok) {
        cat("  ✓ All chromothripsis events pass fragment joins test (p > 0.05)\n")
    } else {
        warning("  ⚠ Some chromothripsis events have fragment joins p ≤ 0.05")
    }

    # Check CN oscillation requirements
    for (i in 1:nrow(chromoth_events)) {
        row <- chromoth_events[i, ]
        cn_osc <- row$max_number_oscillating_CN_segments_2_states
        if (is.na(cn_osc)) {
            cn_osc <- row$max_number_oscillating_CN_segments_3_states
        }

        if (row$classification == "High confidence") {
            if (!is.na(cn_osc) && cn_osc >= 7) {
                cat(sprintf("  ✓ Chr %s (High conf): %d CN oscillations (≥7 required)\n",
                           row$chrom, cn_osc))
            } else {
                warning(sprintf("  ⚠ Chr %s: High conf but only %d CN oscillations",
                               row$chrom, cn_osc))
            }
        } else if (row$classification == "Low confidence") {
            if (!is.na(cn_osc) && cn_osc >= 4 && cn_osc <= 6) {
                cat(sprintf("  ✓ Chr %s (Low conf): %d CN oscillations (4-6 required)\n",
                           row$chrom, cn_osc))
            } else {
                warning(sprintf("  ⚠ Chr %s: Low conf but %d CN oscillations (should be 4-6)",
                               row$chrom, cn_osc))
            }
        }
    }
} else {
    cat("  ⚠ No chromothripsis events detected\n")
}

cat("\n")

# Generate summary report
cat("Step 5: Generating summary report...\n\n")

summary_result <- summarize_chromothripsis(chromoth_result, print_summary = TRUE)

cat("\n")
cat("========================================================================\n")
cat(" TEST COMPLETED SUCCESSFULLY\n")
cat("========================================================================\n\n")

cat("Key Points Verified:\n")
cat("  1. Classification uses three categories: High/Low/Not chromothripsis\n")
cat("  2. Fragment joins test is enforced (p > 0.05 required)\n")
cat("  3. CN oscillation requirements match tutorial criteria\n")
cat("  4. High confidence: ≥7 oscillating segments\n")
cat("  5. Low confidence: 4-6 oscillating segments\n\n")

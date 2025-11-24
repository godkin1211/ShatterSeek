#!/usr/bin/env Rscript
#
# Test script to verify the classify_mixed_mechanisms bug fix
#

library(ShatterSeek)

cat("\n")
cat("================================================================================\n")
cat("  Testing classify_mixed_mechanisms Bug Fix\n")
cat("================================================================================\n\n")

# Load example data
cat("Loading example data (DO17373)...\n")
data(DO17373)

SV_data <- SVs(
    chrom1 = as.character(SV_DO17373$chrom1),
    pos1 = as.numeric(SV_DO17373$start1),
    chrom2 = as.character(SV_DO17373$chrom2),
    pos2 = as.numeric(SV_DO17373$end2),
    SVtype = as.character(SV_DO17373$svclass),
    strand1 = as.character(SV_DO17373$strand1),
    strand2 = as.character(SV_DO17373$strand2)
)

CN_data <- CNVsegs(
    chrom = as.character(SCNA_DO17373$chromosome),
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

cat("  ✓ Data loaded successfully\n\n")

# Run comprehensive analysis
cat("Running comprehensive chromoanagenesis detection...\n")
results <- detect_chromoanagenesis(SV_data, CN_data, genome = "hg19")
cat("  ✓ Detection complete\n\n")

# Test the fixed classify_mixed_mechanisms function
cat("Testing classify_mixed_mechanisms (this previously failed)...\n")
tryCatch({
    mixed_class <- classify_mixed_mechanisms(results)
    cat("  ✓ classify_mixed_mechanisms completed successfully!\n\n")

    # Print results
    cat("Results:\n")
    cat(sprintf("  Classification: %s\n", mixed_class$sample_classification$classification))
    cat(sprintf("  Category: %s\n", mixed_class$sample_classification$category))
    cat(sprintf("  Number of mechanisms: %d\n", mixed_class$sample_classification$n_mechanisms))
    cat(sprintf("  Mechanisms present: %s\n",
               paste(mixed_class$sample_classification$mechanisms_present, collapse = ", ")))
    cat(sprintf("  Complexity score: %.3f\n", mixed_class$complexity$overall_score))
    cat(sprintf("  Complexity level: %s\n", mixed_class$complexity$complexity_level))
    cat("\n")

    cat("✓ BUG FIX VERIFIED - classify_mixed_mechanisms now works correctly!\n\n")

}, error = function(e) {
    cat("✗ ERROR: classify_mixed_mechanisms still failing!\n")
    cat(sprintf("Error message: %s\n", e$message))
    cat("\nPlease report this error.\n\n")
    quit(status = 1)
})

cat("================================================================================\n")
cat("  Test Complete\n")
cat("================================================================================\n\n")

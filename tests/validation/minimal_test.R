#!/usr/bin/env Rscript
# Minimal validation test - Quick proof that core functions work
# Run time: ~10 seconds

cat("\n")
cat("================================================\n")
cat("ShatterSeek Extended Edition - Minimal Test\n")
cat("================================================\n\n")

# Load package
cat("Loading ShatterSeek... ")
suppressPackageStartupMessages({
    library(ShatterSeek)
})
cat("✓\n")

# Test 1: Load data
cat("Loading built-in dataset... ")
data(DO17373, package = "ShatterSeek")

sv_data <- SVs(
    chrom1 = as.character(SV_DO17373$chrom1),
    pos1 = as.numeric(SV_DO17373$start1),
    chrom2 = as.character(SV_DO17373$chrom2),
    pos2 = as.numeric(SV_DO17373$end2),
    SVtype = as.character(SV_DO17373$svclass),
    strand1 = as.character(SV_DO17373$strand1),
    strand2 = as.character(SV_DO17373$strand2)
)

cnv_data <- CNVsegs(
    chrom = as.character(SCNA_DO17373$chromosome),
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

cat(sprintf("✓ (%d SVs, %d CNVs)\n", length(sv_data), length(cnv_data)))

# Test 2: Chromothripsis detection
cat("Running chromothripsis detection... ")
result <- detect_chromothripsis(
    SV.sample = sv_data,
    seg.sample = cnv_data,
    genome = "hg19"
)

if (!is.null(result) && nrow(result@chromSummary) > 0) {
    cat(sprintf("✓ (detected on %d chr)\n", nrow(result@chromSummary)))
} else {
    cat("✗ FAILED\n")
    quit(status = 1)
}

# Test 3: Extended edition
cat("Running extended chromoanagenesis detection... ")
results <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    genome = "hg19"
)

if (!is.null(results) && !is.null(results$chromothripsis)) {
    cat("✓\n")
} else {
    cat("✗ FAILED\n")
    quit(status = 1)
}

# Summary
cat("\n")
cat("================================================\n")
cat("✓ MINIMAL TEST PASSED\n")
cat("================================================\n")
cat("\n")
cat("ShatterSeek Extended Edition is working!\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Run full validation: Rscript phase1_validation.R\n")
cat("  2. Try with your data: see example_usage.R\n")
cat("  3. Generate plots: see README.md\n")
cat("\n")

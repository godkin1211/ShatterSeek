#!/usr/bin/env Rscript
# Quick validation check for ShatterSeek Extended Edition
# Fast sanity check during development

suppressPackageStartupMessages({
    library(ShatterSeek)
})

cat("\n=== Quick Validation Check ===\n\n")

success <- TRUE

# Quick test 1: Load built-in data
cat("1. Loading built-in dataset... ")
tryCatch({
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
    cat("✓\n")
}, error = function(e) {
    cat("✗\n")
    cat("  Error:", as.character(e), "\n")
    success <<- FALSE
})

# Quick test 2: Chromothripsis detection
cat("2. Testing chromothripsis detection... ")
tryCatch({
    result <- detect_chromothripsis(
        SV.sample = sv_data,
        seg.sample = cnv_data,
        genome = "hg19"
    )
    stopifnot(!is.null(result))
    stopifnot(nrow(result@chromSummary) > 0)
    cat("✓\n")
}, error = function(e) {
    cat("✗\n")
    cat("  Error:", as.character(e), "\n")
    success <<- FALSE
})

# Quick test 3: Extended edition detection
cat("3. Testing comprehensive detection... ")
tryCatch({
    results <- detect_chromoanagenesis(
        SV.sample = sv_data,
        CNV.sample = cnv_data,
        genome = "hg19"
    )
    stopifnot(!is.null(results))
    cat("✓\n")
}, error = function(e) {
    cat("✗\n")
    cat("  Error:", as.character(e), "\n")
    success <<- FALSE
})

# Quick test 4: BND parsing
cat("4. Testing BND parsing... ")
tryCatch({
    # Test regex pattern
    test_alts <- c(
        "[chr13:49291490[TT",
        "T]chr1:245088964]",
        "]chr5:15886744]T"
    )
    pattern <- "(\\[|\\])([^:]+):([0-9]+)(\\[|\\])"

    for (alt in test_alts) {
        match <- regexec(pattern, alt)
        stopifnot(match[[1]][1] != -1)
    }
    cat("✓\n")
}, error = function(e) {
    cat("✗\n")
    cat("  Error:", as.character(e), "\n")
    success <<- FALSE
})

# Summary
cat("\n")
if (success) {
    cat("✓ Quick check PASSED - All core functions working\n")
    cat("  Run 'Rscript phase1_validation.R' for comprehensive validation\n\n")
    quit(status = 0)
} else {
    cat("✗ Quick check FAILED - See errors above\n\n")
    quit(status = 1)
}

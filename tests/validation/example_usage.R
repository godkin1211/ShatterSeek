#!/usr/bin/env Rscript
# Example usage of ShatterSeek Extended Edition validation tools

library(ShatterSeek)

# Load synthetic data generator functions
source("synthetic_data_generator.R")

cat("===========================================\n")
cat("ShatterSeek Extended Edition\n")
cat("Validation Examples\n")
cat("===========================================\n\n")

# =============================================================================
# Example 1: Test with built-in dataset
# =============================================================================

cat("Example 1: Built-in Dataset (DO17373)\n")
cat("---------------------------------------\n")

# Load data
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

cat(sprintf("Loaded %d SVs and %d CNV segments\n", length(sv_data), length(cnv_data)))

# Detect chromoanagenesis
cat("\nRunning chromoanagenesis detection...\n")
results <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    genome = "hg19",
    detect_chromothripsis = TRUE,
    detect_chromoplexy = TRUE,
    detect_chromosynthesis = TRUE
)

# Summary
cat("\nResults:\n")
if (!is.null(results$chromothripsis)) {
    cat(sprintf("  Chromothripsis: Detected on %d chromosomes\n",
                nrow(results$chromothripsis$classification)))
}
if (!is.null(results$chromoplexy)) {
    cat(sprintf("  Chromoplexy: %d chains detected\n",
                length(results$chromoplexy$chains)))
}
if (!is.null(results$chromosynthesis)) {
    cat(sprintf("  Chromosynthesis: %d regions detected\n",
                nrow(results$chromosynthesis$summary)))
}

cat("\n")

# =============================================================================
# Example 2: Synthetic chromothripsis
# =============================================================================

cat("Example 2: Synthetic Chromothripsis\n")
cat("------------------------------------\n")

# Generate synthetic data
synthetic_ct <- generate_chromothripsis(
    chromosome = "1",
    region_start = 50e6,
    region_end = 100e6,
    num_breakpoints = 20,
    seed = 12345
)

cat(sprintf("Generated synthetic chromothripsis:\n"))
cat(sprintf("  Chromosome: %s\n", synthetic_ct$metadata$chromosome))
cat(sprintf("  Region: %d - %d (%.1f Mb)\n",
            synthetic_ct$metadata$region_start,
            synthetic_ct$metadata$region_end,
            (synthetic_ct$metadata$region_end - synthetic_ct$metadata$region_start) / 1e6))
cat(sprintf("  Breakpoints: %d\n", synthetic_ct$metadata$num_breakpoints))
cat(sprintf("  SVs: %d\n", nrow(synthetic_ct$SV)))
cat(sprintf("  CNV segments: %d\n", nrow(synthetic_ct$CNV)))

# Convert to ShatterSeek objects
sv_obj <- SVs(
    chrom1 = synthetic_ct$SV$chrom1,
    pos1 = synthetic_ct$SV$pos1,
    chrom2 = synthetic_ct$SV$chrom2,
    pos2 = synthetic_ct$SV$pos2,
    SVtype = synthetic_ct$SV$SVtype,
    strand1 = synthetic_ct$SV$strand1,
    strand2 = synthetic_ct$SV$strand2
)

cnv_obj <- CNVsegs(
    chrom = synthetic_ct$CNV$chrom,
    start = synthetic_ct$CNV$start,
    end = synthetic_ct$CNV$end,
    total_cn = synthetic_ct$CNV$total_cn
)

# Detect
cat("\nDetecting chromothripsis in synthetic data...\n")
ct_result <- detect_chromothripsis(
    SV.sample = sv_obj,
    seg.sample = cnv_obj,
    genome = "hg19"
)

if (!is.null(ct_result) && nrow(ct_result@chromSummary) > 0) {
    cat("✓ Successfully detected chromothripsis in synthetic data\n")
    print(ct_result@chromSummary)
} else {
    cat("✗ Failed to detect chromothripsis (may need parameter tuning)\n")
}

cat("\n")

# =============================================================================
# Example 3: Synthetic chromoplexy
# =============================================================================

cat("Example 3: Synthetic Chromoplexy\n")
cat("---------------------------------\n")

# Generate synthetic chromoplexy
synthetic_cp <- generate_chromoplexy(
    chromosomes = c("1", "2", "3", "4"),
    num_translocations = 6,
    seed = 54321
)

cat(sprintf("Generated synthetic chromoplexy:\n"))
cat(sprintf("  Chromosomes: %s\n", paste(synthetic_cp$metadata$chromosomes, collapse = ", ")))
cat(sprintf("  Translocations: %d\n", synthetic_cp$metadata$num_translocations))
cat(sprintf("  SVs: %d\n", nrow(synthetic_cp$SV)))

# Show translocation chain
cat("\nTranslocation chain:\n")
for (i in 1:min(5, nrow(synthetic_cp$SV))) {
    cat(sprintf("  %d. %s:%d -> %s:%d\n",
                i,
                synthetic_cp$SV$chrom1[i], synthetic_cp$SV$pos1[i],
                synthetic_cp$SV$chrom2[i], synthetic_cp$SV$pos2[i]))
}

cat("\n")

# =============================================================================
# Example 4: Mixed chromoanagenesis
# =============================================================================

cat("Example 4: Mixed Chromoanagenesis\n")
cat("----------------------------------\n")

# Generate mixed event
synthetic_mixed <- generate_mixed_chromoanagenesis(seed = 99999)

cat(sprintf("Generated mixed chromoanagenesis:\n"))
cat(sprintf("  Components: %s\n", paste(synthetic_mixed$metadata$components, collapse = " + ")))
cat(sprintf("  Total SVs: %d\n", nrow(synthetic_mixed$SV)))
cat(sprintf("  Total CNV segments: %d\n", nrow(synthetic_mixed$CNV)))

# Convert and detect
sv_obj <- SVs(
    chrom1 = synthetic_mixed$SV$chrom1,
    pos1 = synthetic_mixed$SV$pos1,
    chrom2 = synthetic_mixed$SV$chrom2,
    pos2 = synthetic_mixed$SV$pos2,
    SVtype = synthetic_mixed$SV$SVtype,
    strand1 = synthetic_mixed$SV$strand1,
    strand2 = synthetic_mixed$SV$strand2
)

cnv_obj <- CNVsegs(
    chrom = synthetic_mixed$CNV$chrom,
    start = synthetic_mixed$CNV$start,
    end = synthetic_mixed$CNV$end,
    total_cn = synthetic_mixed$CNV$total_cn
)

cat("\nDetecting multiple mechanisms...\n")
mixed_results <- detect_chromoanagenesis(
    SV.sample = sv_obj,
    CNV.sample = cnv_obj,
    genome = "hg19"
)

cat("\nDetection results:\n")
detected <- c()
if (!is.null(mixed_results$chromothripsis)) {
    detected <- c(detected, "chromothripsis")
}
if (!is.null(mixed_results$chromoplexy)) {
    detected <- c(detected, "chromoplexy")
}
if (!is.null(mixed_results$chromosynthesis)) {
    detected <- c(detected, "chromosynthesis")
}

if (length(detected) > 0) {
    cat(sprintf("  Detected: %s\n", paste(detected, collapse = ", ")))
} else {
    cat("  No mechanisms detected\n")
}

cat("\n")
cat("===========================================\n")
cat("Examples completed\n")
cat("===========================================\n")

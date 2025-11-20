#!/usr/bin/env Rscript
# Phase 1 Validation for ShatterSeek Extended Edition
# Tests: Built-in data, synthetic data, BND parsing

library(ShatterSeek)

cat("===========================================\n")
cat("ShatterSeek Extended Edition - Phase 1 Validation\n")
cat("===========================================\n\n")

# Initialize results storage
validation_results <- list(
    timestamp = Sys.time(),
    version = packageVersion("ShatterSeek"),
    tests_passed = 0,
    tests_failed = 0,
    tests_total = 0,
    details = list()
)

#' Test wrapper with error handling
run_test <- function(test_name, test_func) {
    validation_results$tests_total <<- validation_results$tests_total + 1

    cat(sprintf("\n[Test %d] %s\n", validation_results$tests_total, test_name))
    cat(strrep("-", 60), "\n")

    result <- tryCatch({
        test_func()
        validation_results$tests_passed <<- validation_results$tests_passed + 1
        list(status = "PASS", error = NULL)
    }, error = function(e) {
        validation_results$tests_failed <<- validation_results$tests_failed + 1
        list(status = "FAIL", error = as.character(e))
    })

    validation_results$details[[test_name]] <<- result

    if (result$status == "PASS") {
        cat("✓ PASSED\n")
    } else {
        cat("✗ FAILED\n")
        cat("  Error:", result$error, "\n")
    }

    return(result$status == "PASS")
}

# =============================================================================
# TEST 1: Built-in Dataset (DO17373)
# =============================================================================

test_builtin_data <- function() {
    cat("Loading built-in dataset DO17373...\n")
    data(DO17373, package = "ShatterSeek")

    # Prepare SV data
    sv_data <- SVs(
        chrom1 = as.character(SV_DO17373$chrom1),
        pos1 = as.numeric(SV_DO17373$start1),
        chrom2 = as.character(SV_DO17373$chrom2),
        pos2 = as.numeric(SV_DO17373$end2),
        SVtype = as.character(SV_DO17373$svclass),
        strand1 = as.character(SV_DO17373$strand1),
        strand2 = as.character(SV_DO17373$strand2)
    )

    # Prepare CNV data
    cnv_data <- CNVsegs(
        chrom = as.character(SCNA_DO17373$chromosome),
        start = SCNA_DO17373$start,
        end = SCNA_DO17373$end,
        total_cn = SCNA_DO17373$total_cn
    )

    cat(sprintf("  SV records: %d\n", length(sv_data)))
    cat(sprintf("  CNV segments: %d\n", length(cnv_data)))

    # Test original chromothripsis detection
    cat("\n  Testing chromothripsis detection...\n")
    ct_result <- detect_chromothripsis(
        SV.sample = sv_data,
        seg.sample = cnv_data,
        genome = "hg19"
    )

    stopifnot(!is.null(ct_result))
    stopifnot(length(ct_result@chromSummary) > 0)
    cat(sprintf("  ✓ Detected chromothripsis on %d chromosomes\n",
                nrow(ct_result@chromSummary)))

    # Test extended edition comprehensive detection
    cat("\n  Testing comprehensive chromoanagenesis detection...\n")
    results <- detect_chromoanagenesis(
        SV.sample = sv_data,
        CNV.sample = cnv_data,
        genome = "hg19",
        detect_chromothripsis = TRUE,
        detect_chromoplexy = TRUE,
        detect_chromosynthesis = TRUE
    )

    stopifnot(!is.null(results))
    stopifnot(!is.null(results$chromothripsis))

    cat("  ✓ Comprehensive detection successful\n")

    # Test classification
    cat("\n  Testing classification...\n")
    classification <- classify_chromothripsis(
        chromoth_output = ct_result,
        min_cluster_size = 6,
        min_cn_oscillations = 4
    )

    stopifnot(nrow(classification) > 0)
    cat(sprintf("  ✓ Classification complete: %d chromosomes classified\n",
                nrow(classification)))

    # Verify consistency between original and extended
    cat("\n  Verifying consistency between original and extended...\n")
    original_chroms <- ct_result@chromSummary$chrom
    extended_chroms <- results$chromothripsis$detection_output@chromSummary$chrom

    stopifnot(all(original_chroms %in% extended_chroms))
    cat("  ✓ Original and extended results are consistent\n")

    cat("\n✓ Built-in data test PASSED\n")
    return(TRUE)
}

# =============================================================================
# TEST 2: Synthetic Data Generation and Detection
# =============================================================================

#' Generate synthetic chromothripsis data
generate_synthetic_chromothripsis <- function(
    chromosome = "1",
    region_start = 50e6,
    region_end = 100e6,
    num_breakpoints = 20,
    seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)

    # Generate clustered breakpoints
    breakpoints <- sort(sample(region_start:region_end, num_breakpoints))

    # Create SVs
    sv_list <- list()
    for (i in seq(1, length(breakpoints)-1, 2)) {
        if (i+1 > length(breakpoints)) break

        sv_type <- sample(c("DEL", "h2hINV", "t2tINV"), 1, prob = c(0.5, 0.25, 0.25))

        # Assign correct strands based on SV type
        if (sv_type == "DEL") {
            strand1 <- "+"
            strand2 <- "-"
        } else if (sv_type == "h2hINV") {
            strand1 <- "+"
            strand2 <- "+"
        } else {  # t2tINV
            strand1 <- "-"
            strand2 <- "-"
        }

        sv_list[[length(sv_list) + 1]] <- data.frame(
            chrom1 = chromosome,
            pos1 = breakpoints[i],
            chrom2 = chromosome,
            pos2 = breakpoints[i+1],
            SVtype = sv_type,
            strand1 = strand1,
            strand2 = strand2,
            stringsAsFactors = FALSE
        )
    }

    sv_data <- do.call(rbind, sv_list)

    # Generate oscillating CN segments
    num_segments <- num_breakpoints + 1
    cn_values <- rep(c(1, 3, 1, 4, 2, 3, 1, 2), length.out = num_segments)

    cnv_list <- list()
    seg_starts <- c(region_start, breakpoints)
    seg_ends <- c(breakpoints, region_end)

    for (i in 1:num_segments) {
        cnv_list[[i]] <- data.frame(
            chrom = chromosome,
            start = seg_starts[i],
            end = seg_ends[i],
            total_cn = cn_values[i],
            stringsAsFactors = FALSE
        )
    }

    cnv_data <- do.call(rbind, cnv_list)

    return(list(SV = sv_data, CNV = cnv_data))
}

test_synthetic_data <- function() {
    cat("Generating synthetic chromothripsis data...\n")

    n_tests <- 10
    detected <- 0

    for (i in 1:n_tests) {
        cat(sprintf("  Test case %d/%d: ", i, n_tests))

        # Generate synthetic data
        synthetic <- generate_synthetic_chromothripsis(
            chromosome = sample(c("1", "2", "3"), 1),
            num_breakpoints = sample(15:25, 1),
            seed = 12345 + i
        )

        # Convert to ShatterSeek objects
        sv_obj <- SVs(
            chrom1 = synthetic$SV$chrom1,
            pos1 = synthetic$SV$pos1,
            chrom2 = synthetic$SV$chrom2,
            pos2 = synthetic$SV$pos2,
            SVtype = synthetic$SV$SVtype,
            strand1 = synthetic$SV$strand1,
            strand2 = synthetic$SV$strand2
        )

        cnv_obj <- CNVsegs(
            chrom = synthetic$CNV$chrom,
            start = synthetic$CNV$start,
            end = synthetic$CNV$end,
            total_cn = synthetic$CNV$total_cn
        )

        # Detect chromothripsis
        result <- tryCatch({
            detect_chromothripsis(
                SV.sample = sv_obj,
                seg.sample = cnv_obj,
                genome = "hg19"
            )
        }, error = function(e) {
            cat("ERROR\n")
            return(NULL)
        })

        if (!is.null(result) && nrow(result@chromSummary) > 0) {
            detected <- detected + 1
            cat("DETECTED\n")
        } else {
            cat("NOT DETECTED\n")
        }
    }

    sensitivity <- detected / n_tests
    cat(sprintf("\n  Sensitivity on synthetic data: %.1f%% (%d/%d)\n",
                sensitivity * 100, detected, n_tests))

    # Require at least 70% detection rate
    stopifnot(sensitivity >= 0.7)

    cat("\n✓ Synthetic data test PASSED\n")
    return(TRUE)
}

# =============================================================================
# TEST 3: BND/TRA Parsing
# =============================================================================

test_bnd_parsing <- function() {
    cat("Testing BND ALT field parsing...\n")

    # Create test VCF content with BND records
    vcf_content <- c(
        "##fileformat=VCFv4.2",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##ALT=<ID=BND,Description=\"Breakend\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t100000\tbnd1\tN\t[chr13:50000000[A\t.\tPASS\tSVTYPE=BND",
        "chr13\t50000000\tbnd2\tN\t]chr1:100000]T\t.\tPASS\tSVTYPE=BND",
        "chr2\t200000\tbnd3\tN\tT]chr5:30000000]\t.\tPASS\tSVTYPE=BND",
        "chr5\t30000000\tbnd4\tN\tG[chr2:200000[\t.\tPASS\tSVTYPE=BND"
    )

    # Write to temporary file
    tmp_vcf <- tempfile(fileext = ".vcf")
    writeLines(vcf_content, tmp_vcf)

    cat(sprintf("  Created test VCF: %s\n", tmp_vcf))

    # Test parsing
    cat("  Parsing VCF...\n")
    sv_data <- tryCatch({
        read_sv_vcf(tmp_vcf, caller = "generic", min_sv_size = 0, include_tra = TRUE)
    }, error = function(e) {
        cat("  Error during parsing:", as.character(e), "\n")
        NULL
    })

    # Clean up
    unlink(tmp_vcf)

    stopifnot(!is.null(sv_data))
    cat(sprintf("  ✓ Parsed %d SV records\n", length(sv_data@chrom1)))

    # Check for TRA records
    tra_idx <- which(sv_data@SVtype == "TRA")
    cat(sprintf("  ✓ Found %d TRA records\n", length(tra_idx)))

    stopifnot(length(tra_idx) > 0)

    # Verify mate coordinates
    cat("\n  Verifying mate coordinates...\n")
    for (i in tra_idx) {
        cat(sprintf("    TRA %d: %s:%d -> %s:%d\n",
                    which(tra_idx == i),
                    sv_data@chrom1[i], sv_data@pos1[i],
                    sv_data@chrom2[i], sv_data@pos2[i]))

        # Should be interchromosomal
        stopifnot(sv_data@chrom1[i] != sv_data@chrom2[i])
    }

    cat("  ✓ All TRA records are interchromosomal\n")

    # Test different BND formats
    cat("\n  Testing various BND bracket formats...\n")

    test_formats <- list(
        list(alt = "[chr13:49291490[TT", expected_chr = "13", expected_pos = 49291490),
        list(alt = "T]chr1:245088964]", expected_chr = "1", expected_pos = 245088964),
        list(alt = "]chr5:15886744]T", expected_chr = "5", expected_pos = 15886744),
        list(alt = "G[chr2:200000[", expected_chr = "2", expected_pos = 200000)
    )

    for (fmt in test_formats) {
        # Test regex pattern
        pattern <- "(\\[|\\])([^:]+):([0-9]+)(\\[|\\])"
        match <- regexec(pattern, fmt$alt)

        if (match[[1]][1] != -1) {
            matches <- regmatches(fmt$alt, match)[[1]]
            mate_chr <- gsub("^chr", "", matches[3], ignore.case = TRUE)
            mate_pos <- as.numeric(matches[4])

            cat(sprintf("    '%s' -> chr=%s, pos=%d ", fmt$alt, mate_chr, mate_pos))

            stopifnot(mate_chr == fmt$expected_chr)
            stopifnot(mate_pos == fmt$expected_pos)

            cat("✓\n")
        } else {
            stop(sprintf("Failed to parse: %s", fmt$alt))
        }
    }

    cat("\n✓ BND parsing test PASSED\n")
    return(TRUE)
}

# =============================================================================
# Run All Tests
# =============================================================================

cat("\nStarting Phase 1 Validation Tests...\n")
cat("====================================\n")

# Test 1: Built-in data
run_test("Built-in Dataset (DO17373)", test_builtin_data)

# Test 2: Synthetic data
run_test("Synthetic Data Generation & Detection", test_synthetic_data)

# Test 3: BND parsing
run_test("BND/TRA ALT Field Parsing", test_bnd_parsing)

# =============================================================================
# Summary Report
# =============================================================================

cat("\n")
cat("===========================================\n")
cat("VALIDATION SUMMARY\n")
cat("===========================================\n")
cat(sprintf("Total tests: %d\n", validation_results$tests_total))
cat(sprintf("Passed: %d\n", validation_results$tests_passed))
cat(sprintf("Failed: %d\n", validation_results$tests_failed))
cat(sprintf("Success rate: %.1f%%\n",
            100 * validation_results$tests_passed / validation_results$tests_total))
cat("===========================================\n")

if (validation_results$tests_failed > 0) {
    cat("\nFailed tests:\n")
    for (test_name in names(validation_results$details)) {
        if (validation_results$details[[test_name]]$status == "FAIL") {
            cat(sprintf("  - %s\n", test_name))
            cat(sprintf("    Error: %s\n", validation_results$details[[test_name]]$error))
        }
    }
}

# Save results
results_file <- "phase1_validation_results.rds"
saveRDS(validation_results, results_file)
cat(sprintf("\nResults saved to: %s\n", results_file))

# Exit with appropriate code
if (validation_results$tests_failed == 0) {
    cat("\n✓ ALL TESTS PASSED\n\n")
    quit(status = 0)
} else {
    cat("\n✗ SOME TESTS FAILED\n\n")
    quit(status = 1)
}

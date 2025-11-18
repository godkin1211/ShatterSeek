#' Comprehensive chromoanagenesis detection
#'
#' Performs integrated analysis to detect chromothripsis, chromoplexy, and
#' chromosynthesis from the same dataset, providing a complete view of complex
#' chromosomal rearrangements.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample An instance of class CNVsegs or data frame with CNV data
#' @param genome Reference genome ("hg19" or "hg38", default: "hg19")
#' @param detect_chromothripsis Detect chromothripsis events (default: TRUE)
#' @param detect_chromoplexy Detect chromoplexy events (default: TRUE)
#' @param detect_chromosynthesis Detect chromosynthesis events (default: TRUE)
#' @param min_chromothripsis_size Minimum SV cluster size for chromothripsis (default: 1)
#' @param min_chromoplexy_chromosomes Minimum chromosomes for chromoplexy (default: 3)
#' @param min_chromosynthesis_tandem_dups Minimum tandem dups for chromosynthesis (default: 3)
#' @param verbose Print progress messages (default: TRUE)
#' @return A list containing results for all three analyses
#' @details
#' This function provides a comprehensive analysis of chromoanagenesis events
#' by detecting chromothripsis, chromoplexy, and chromosynthesis in the same sample.
#'
#' Chromothripsis characteristics:
#' - Localized to one or few chromosomes
#' - Many clustered breakpoints
#' - Oscillating copy number patterns
#' - Random fragment joins
#'
#' Chromoplexy characteristics:
#' - Involves multiple chromosomes (3-8)
#' - Chained translocations
#' - Minimal copy number changes
#' - May form cycles
#'
#' Chromosynthesis characteristics:
#' - Replication-based mechanism (FoSTeS/MMBIR)
#' - Gradual copy number increases
#' - Tandem duplications and inversions
#' - Localized to specific regions
#'
#' The function returns integrated results allowing comparison and
#' classification of all complex rearrangement patterns.
#'
#' @examples
#' \dontrun{
#' library(ShatterSeek)
#' data(DO17373)
#'
#' # Prepare data
#' SV_data <- SVs(chrom1=SV_DO17373$chrom1, ...)
#' CN_data <- CNVsegs(chrom=SCNA_DO17373$chromosome, ...)
#'
#' # Run comprehensive analysis
#' results <- detect_chromoanagenesis(SV_data, CN_data)
#'
#' # View results
#' print(results)
#' summary(results)
#'
#' # Access specific results
#' results$chromothripsis
#' results$chromoplexy
#' }
#'
#' @export
detect_chromoanagenesis <- function(SV.sample,
                                   CNV.sample,
                                   genome = "hg19",
                                   detect_chromothripsis = TRUE,
                                   detect_chromoplexy = TRUE,
                                   detect_chromosynthesis = TRUE,
                                   min_chromothripsis_size = 1,
                                   min_chromoplexy_chromosomes = 3,
                                   min_chromosynthesis_tandem_dups = 3,
                                   verbose = TRUE) {

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("     COMPREHENSIVE CHROMOANAGENESIS ANALYSIS\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    results <- list()

    # 1. Data quality check
    if (verbose) cat("Step 1: Checking data quality...\n")
    quality <- check_data_quality(SV.sample, CNV.sample, verbose = FALSE)

    if (quality$has_issues && verbose) {
        cat("  WARNING: Data quality issues detected.\n")
        cat(sprintf("  Number of warnings: %d\n", length(quality$warnings)))
    } else if (verbose) {
        cat("  Data quality: OK\n")
    }

    results$quality_check <- quality

    # 2. Detect chromothripsis
    if (detect_chromothripsis) {
        if (verbose) {
            cat("\nStep 2: Detecting chromothripsis...\n")
        }

        chromothripsis_result <- shatterseek(
            SV.sample = SV.sample,
            seg.sample = CNV.sample,
            min.Size = min_chromothripsis_size,
            genome = genome
        )

        # Calculate scores and classifications
        if (verbose) cat("  Calculating confidence scores...\n")
        chromothripsis_classification <- classify_chromothripsis(chromothripsis_result)

        results$chromothripsis <- list(
            shatterseek_output = chromothripsis_result,
            classification = chromothripsis_classification,
            n_likely = sum(chromothripsis_classification$classification == "Likely chromothripsis"),
            n_possible = sum(chromothripsis_classification$classification == "Possible chromothripsis")
        )

        if (verbose) {
            cat(sprintf("  Found %d likely and %d possible chromothripsis events.\n",
                       results$chromothripsis$n_likely,
                       results$chromothripsis$n_possible))
        }
    }

    # 3. Detect chromoplexy
    if (detect_chromoplexy) {
        if (verbose) {
            cat("\nStep 3: Detecting chromoplexy...\n")
        }

        chromoplexy_result <- detect_chromoplexy(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_chromosomes = min_chromoplexy_chromosomes
        )

        results$chromoplexy <- chromoplexy_result

        if (verbose) {
            cat(sprintf("  Found %d likely and %d possible chromoplexy events.\n",
                       chromoplexy_result$likely_chromoplexy,
                       chromoplexy_result$possible_chromoplexy))
        }
    }

    # 4. Detect chromosynthesis
    if (detect_chromosynthesis) {
        if (verbose) {
            cat("\nStep 4: Detecting chromosynthesis...\n")
        }

        chromosynthesis_result <- detect_chromosynthesis(
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            min_tandem_dups = min_chromosynthesis_tandem_dups
        )

        results$chromosynthesis <- chromosynthesis_result

        if (verbose) {
            cat(sprintf("  Found %d likely and %d possible chromosynthesis events.\n",
                       chromosynthesis_result$likely_chromosynthesis,
                       chromosynthesis_result$possible_chromosynthesis))
        }
    }

    # 5. Integrated classification
    if (verbose) cat("\nStep 5: Integrated classification...\n")

    integrated_summary <- create_integrated_summary(results)
    results$integrated_summary <- integrated_summary

    if (verbose) {
        cat("\n")
        cat(rep("=", 70), "\n", sep = "")
        cat("ANALYSIS COMPLETE\n")
        cat(rep("=", 70), "\n\n", sep = "")
    }

    class(results) <- c("chromoanagenesis", "list")
    return(results)
}


#' Create integrated summary of chromoanagenesis results
#'
#' @param results Results list from detect_chromoanagenesis
#' @return Integrated summary data frame
#' @keywords internal
create_integrated_summary <- function(results) {

    summary <- list()

    # Chromothripsis summary
    if (!is.null(results$chromothripsis)) {
        summary$chromothripsis_likely <- results$chromothripsis$n_likely
        summary$chromothripsis_possible <- results$chromothripsis$n_possible

        if (results$chromothripsis$n_likely > 0) {
            likely_chroms <- results$chromothripsis$classification[
                results$chromothripsis$classification$classification == "Likely chromothripsis",
                "chrom"
            ]
            summary$chromothripsis_chromosomes <- paste(likely_chroms, collapse = ", ")
        } else {
            summary$chromothripsis_chromosomes <- "None"
        }
    } else {
        summary$chromothripsis_likely <- NA
        summary$chromothripsis_possible <- NA
        summary$chromothripsis_chromosomes <- "Not analyzed"
    }

    # Chromoplexy summary
    if (!is.null(results$chromoplexy)) {
        summary$chromoplexy_likely <- results$chromoplexy$likely_chromoplexy
        summary$chromoplexy_possible <- results$chromoplexy$possible_chromoplexy
        summary$chromoplexy_chains <- results$chromoplexy$total_chains

        if (results$chromoplexy$likely_chromoplexy > 0) {
            likely_chains <- results$chromoplexy$summary[
                results$chromoplexy$summary$classification == "Likely chromoplexy",
            ]
            chr_list <- unique(unlist(strsplit(likely_chains$chromosomes_involved, ",")))
            summary$chromoplexy_chromosomes <- paste(chr_list, collapse = ", ")
        } else {
            summary$chromoplexy_chromosomes <- "None"
        }
    } else {
        summary$chromoplexy_likely <- NA
        summary$chromoplexy_possible <- NA
        summary$chromoplexy_chains <- NA
        summary$chromoplexy_chromosomes <- "Not analyzed"
    }

    # Chromosynthesis summary
    if (!is.null(results$chromosynthesis)) {
        summary$chromosynthesis_likely <- results$chromosynthesis$likely_chromosynthesis
        summary$chromosynthesis_possible <- results$chromosynthesis$possible_chromosynthesis
        summary$chromosynthesis_regions <- results$chromosynthesis$total_regions

        if (results$chromosynthesis$likely_chromosynthesis > 0) {
            likely_regions <- results$chromosynthesis$summary[
                results$chromosynthesis$summary$classification == "Likely chromosynthesis",
            ]
            chr_list <- unique(likely_regions$chrom)
            summary$chromosynthesis_chromosomes <- paste(chr_list, collapse = ", ")
        } else {
            summary$chromosynthesis_chromosomes <- "None"
        }
    } else {
        summary$chromosynthesis_likely <- NA
        summary$chromosynthesis_possible <- NA
        summary$chromosynthesis_regions <- NA
        summary$chromosynthesis_chromosomes <- "Not analyzed"
    }

    # Overall assessment
    has_chromothripsis <- !is.null(results$chromothripsis) &&
                         results$chromothripsis$n_likely > 0
    has_chromoplexy <- !is.null(results$chromoplexy) &&
                      results$chromoplexy$likely_chromoplexy > 0
    has_chromosynthesis <- !is.null(results$chromosynthesis) &&
                          results$chromosynthesis$likely_chromosynthesis > 0

    # Create detailed classification
    mechanisms <- c()
    if (has_chromothripsis) mechanisms <- c(mechanisms, "chromothripsis")
    if (has_chromoplexy) mechanisms <- c(mechanisms, "chromoplexy")
    if (has_chromosynthesis) mechanisms <- c(mechanisms, "chromosynthesis")

    if (length(mechanisms) == 0) {
        summary$overall_classification <- "No chromoanagenesis detected"
    } else if (length(mechanisms) == 1) {
        summary$overall_classification <- paste(tools::toTitleCase(mechanisms[1]), "only")
    } else {
        summary$overall_classification <- paste("Multiple mechanisms:",
                                               paste(mechanisms, collapse = ", "))
    }

    return(as.data.frame(summary, stringsAsFactors = FALSE))
}


#' Print method for chromoanagenesis results
#'
#' @param x Chromoanagenesis result object
#' @param ... Additional arguments
#' @export
print.chromoanagenesis <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOANAGENESIS ANALYSIS SUMMARY\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat("Overall Classification:", x$integrated_summary$overall_classification, "\n\n")

    # Chromothripsis results
    if (!is.null(x$chromothripsis)) {
        cat("CHROMOTHRIPSIS:\n")
        cat(sprintf("  - Likely events:   %d\n", x$chromothripsis$n_likely))
        cat(sprintf("  - Possible events: %d\n", x$chromothripsis$n_possible))
        cat(sprintf("  - Chromosomes:     %s\n",
                   x$integrated_summary$chromothripsis_chromosomes))
        cat("\n")
    }

    # Chromoplexy results
    if (!is.null(x$chromoplexy)) {
        cat("CHROMOPLEXY:\n")
        cat(sprintf("  - Likely events:   %d\n", x$chromoplexy$likely_chromoplexy))
        cat(sprintf("  - Possible events: %d\n", x$chromoplexy$possible_chromoplexy))
        cat(sprintf("  - Total chains:    %d\n", x$chromoplexy$total_chains))
        cat(sprintf("  - Chromosomes:     %s\n",
                   x$integrated_summary$chromoplexy_chromosomes))
        cat("\n")
    }

    # Data quality
    if (!is.null(x$quality_check) && x$quality_check$has_issues) {
        cat("DATA QUALITY WARNINGS:\n")
        for (w in x$quality_check$warnings) {
            cat(sprintf("  - %s\n", w))
        }
        cat("\n")
    }

    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for chromoanagenesis results
#'
#' @param object Chromoanagenesis result object
#' @param ... Additional arguments
#' @export
summary.chromoanagenesis <- function(object, ...) {

    cat("\nIntegrated Summary Table:\n")
    print(object$integrated_summary)

    cat("\n")

    if (!is.null(object$chromothripsis)) {
        cat("Chromothripsis Details:\n")
        print(object$chromothripsis$classification[, c("chrom", "classification",
                                                       "confidence_score",
                                                       "clusterSize")])
        cat("\n")
    }

    if (!is.null(object$chromoplexy) && object$chromoplexy$total_chains > 0) {
        cat("Chromoplexy Details:\n")
        print(object$chromoplexy$summary[, c("chain_id", "n_chromosomes",
                                            "n_translocations", "classification")])
        cat("\n")
    }

    invisible(object)
}

#' Detect chromosynthesis events from SV and CNV data
#'
#' Chromosynthesis is characterized by replication-based mechanisms
#' (Fork Stalling and Template Switching - FoSTeS) resulting in
#' gradual copy number increases, tandem duplications, and complex
#' local rearrangements.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample An instance of class CNVsegs or data frame with CNV data
#' @param min_tandem_dups Minimum number of tandem duplications (default: 3)
#' @param min_cn_segments Minimum CN segments for gradient analysis (default: 5)
#' @param gradient_threshold Minimum gradient correlation (default: 0.5)
#' @param max_region_size Maximum region size in bp (default: 10e6)
#' @return A list containing chromosynthesis detection results
#' @details
#' Chromosynthesis differs from chromothripsis and chromoplexy:
#' - Mechanism: Replication-dependent (FoSTeS/MMBIR)
#' - CN pattern: Gradual increases, not oscillations
#' - SV pattern: Tandem duplications, inverted repeats
#' - Localization: Usually 1-2 chromosomes
#' - Directionality: Serial template switching
#'
#' Detection criteria:
#' 1. Gradual CN increase (positive correlation with position)
#' 2. Enrichment of tandem duplications
#' 3. Complex local rearrangements
#' 4. Minimal inter-chromosomal events
#'
#' @references
#' Lee et al. (2015) Nat Genet. Complex chromosomal rearrangements by
#' single catastrophic pathways in cancer genomes.
#'
#' @export
detect_chromosynthesis <- function(SV.sample,
                                  CNV.sample,
                                  min_tandem_dups = 3,
                                  min_cn_segments = 5,
                                  gradient_threshold = 0.5,
                                  max_region_size = 10e6) {

    # Convert to data frames if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (is.null(CNV.sample)) {
        stop("CNV data is required for chromosynthesis detection.")
    }

    if (is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    if (nrow(CNV.sample) < min_cn_segments) {
        warning(sprintf("Insufficient CNV segments (%d). Need at least %d for chromosynthesis detection.",
                       nrow(CNV.sample), min_cn_segments))
        return(create_empty_chromosynthesis_result())
    }

    cat("Detecting chromosynthesis patterns...\n")

    # Step 1: Identify regions with CN gradients
    cat("Step 1: Analyzing copy number gradients...\n")
    gradient_regions <- detect_cn_gradients(
        CNV.sample,
        min_segments = min_cn_segments,
        gradient_threshold = gradient_threshold,
        max_region_size = max_region_size
    )

    if (length(gradient_regions) == 0) {
        cat("No significant CN gradient regions found.\n")
        return(create_empty_chromosynthesis_result())
    }

    cat(sprintf("Found %d region(s) with CN gradients.\n", length(gradient_regions)))

    # Step 2: Detect tandem duplications
    cat("Step 2: Detecting tandem duplications...\n")
    tandem_dup_analysis <- analyze_tandem_duplications(
        SV.sample,
        gradient_regions
    )

    # Step 3: Evaluate each region
    cat("Step 3: Evaluating chromosynthesis criteria...\n")
    region_results <- list()

    for (i in 1:length(gradient_regions)) {
        region_results[[i]] <- evaluate_chromosynthesis_region(
            region = gradient_regions[[i]],
            SV.sample = SV.sample,
            CNV.sample = CNV.sample,
            tandem_dups = tandem_dup_analysis[[i]],
            min_tandem_dups = min_tandem_dups
        )
    }

    # Step 4: Create summary
    summary_df <- do.call(rbind, lapply(region_results, function(x) x$summary))

    # Classify events
    summary_df$classification <- classify_chromosynthesis_event(summary_df)

    result <- list(
        regions = gradient_regions,
        region_details = region_results,
        summary = summary_df,
        total_regions = length(gradient_regions),
        likely_chromosynthesis = sum(summary_df$classification == "Likely chromosynthesis"),
        possible_chromosynthesis = sum(summary_df$classification == "Possible chromosynthesis")
    )

    class(result) <- c("chromosynthesis", "list")

    cat(sprintf("\nDetection complete: %d likely, %d possible chromosynthesis events.\n",
               result$likely_chromosynthesis,
               result$possible_chromosynthesis))

    return(result)
}


#' Detect copy number gradients in chromosomes
#'
#' Identifies regions with gradual CN increases characteristic of chromosynthesis.
#'
#' @param CNV.sample CNV data frame
#' @param min_segments Minimum segments for analysis
#' @param gradient_threshold Minimum correlation coefficient
#' @param max_region_size Maximum region size
#' @return List of gradient regions
#' @keywords internal
detect_cn_gradients <- function(CNV.sample,
                                min_segments = 5,
                                gradient_threshold = 0.5,
                                max_region_size = 10e6) {

    gradient_regions <- list()
    region_id <- 0

    # Analyze each chromosome
    chroms <- unique(CNV.sample$chrom)

    for (chr in chroms) {
        chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]

        if (nrow(chr_cnv) < min_segments) next

        # Sort by position
        chr_cnv <- chr_cnv[order(chr_cnv$start), ]

        # Use sliding window to detect gradient regions
        window_results <- find_gradient_windows(
            chr_cnv,
            min_segments = min_segments,
            gradient_threshold = gradient_threshold,
            max_size = max_region_size
        )

        if (length(window_results) > 0) {
            for (window in window_results) {
                region_id <- region_id + 1
                window$region_id <- region_id
                gradient_regions[[region_id]] <- window
            }
        }
    }

    return(gradient_regions)
}


#' Find windows with CN gradients using sliding window approach
#'
#' @param chr_cnv Chromosome CNV data
#' @param min_segments Minimum segments
#' @param gradient_threshold Threshold for correlation
#' @param max_size Maximum window size
#' @return List of gradient windows
#' @keywords internal
find_gradient_windows <- function(chr_cnv,
                                  min_segments = 5,
                                  gradient_threshold = 0.5,
                                  max_size = 10e6) {

    windows <- list()
    n_segs <- nrow(chr_cnv)

    # Try different window sizes
    for (win_size in min_segments:min(n_segs, min_segments + 10)) {

        for (start_idx in 1:(n_segs - win_size + 1)) {
            end_idx <- start_idx + win_size - 1

            window_cnv <- chr_cnv[start_idx:end_idx, ]

            # Check window size constraint
            region_start <- window_cnv$start[1]
            region_end <- window_cnv$end[nrow(window_cnv)]
            region_size <- region_end - region_start

            if (region_size > max_size) next

            # Calculate gradient correlation
            segment_order <- 1:nrow(window_cnv)
            cn_values <- window_cnv$total_cn

            # Spearman correlation (robust to outliers)
            if (length(unique(cn_values)) > 1) {
                gradient_cor <- cor(segment_order, cn_values, method = "spearman")

                # Look for positive gradients (CN increases)
                if (!is.na(gradient_cor) && gradient_cor >= gradient_threshold) {

                    # Calculate additional metrics
                    cn_range <- max(cn_values) - min(cn_values)
                    cn_mean <- mean(cn_values)
                    cn_trend <- lm(cn_values ~ segment_order)
                    slope <- coef(cn_trend)[2]

                    # Only keep if there's actual increase
                    if (slope > 0.1) {
                        windows[[length(windows) + 1]] <- list(
                            chrom = window_cnv$chrom[1],
                            start = region_start,
                            end = region_end,
                            size = region_size,
                            n_segments = nrow(window_cnv),
                            gradient_correlation = gradient_cor,
                            cn_range = cn_range,
                            cn_mean = cn_mean,
                            slope = slope,
                            segments = window_cnv
                        )
                    }
                }
            }
        }
    }

    # Remove overlapping windows, keep those with highest correlation
    if (length(windows) > 0) {
        windows <- remove_overlapping_windows(windows)
    }

    return(windows)
}


#' Remove overlapping gradient windows
#'
#' @param windows List of windows
#' @return Non-overlapping windows
#' @keywords internal
remove_overlapping_windows <- function(windows) {

    if (length(windows) <= 1) return(windows)

    # Sort by correlation (descending)
    cors <- sapply(windows, function(x) x$gradient_correlation)
    windows <- windows[order(cors, decreasing = TRUE)]

    kept_windows <- list()
    kept_windows[[1]] <- windows[[1]]

    for (i in 2:length(windows)) {
        window <- windows[[i]]

        # Check overlap with kept windows
        overlaps <- FALSE
        for (kept in kept_windows) {
            if (window$chrom == kept$chrom &&
                window$start < kept$end &&
                window$end > kept$start) {
                overlaps <- TRUE
                break
            }
        }

        if (!overlaps) {
            kept_windows[[length(kept_windows) + 1]] <- window
        }
    }

    return(kept_windows)
}


#' Analyze tandem duplications in gradient regions
#'
#' @param SV.sample SV data
#' @param gradient_regions Gradient regions
#' @return List of tandem duplication analyses
#' @keywords internal
analyze_tandem_duplications <- function(SV.sample, gradient_regions) {

    tandem_dup_list <- list()

    for (i in 1:length(gradient_regions)) {
        region <- gradient_regions[[i]]

        # Extract SVs in this region
        region_svs <- SV.sample[
            SV.sample$chrom1 == region$chrom &
            SV.sample$chrom2 == region$chrom &
            SV.sample$pos1 >= region$start &
            SV.sample$pos2 <= region$end,
        ]

        # Identify tandem duplications
        # Tandem DUP: same chromosome, -/+ strands, close proximity
        tandem_dups <- region_svs[
            region_svs$SVtype == "DUP" &
            abs(region_svs$pos2 - region_svs$pos1) < 1e6,  # Within 1 Mb
        ]

        # Check for inverted repeats (inversions)
        inversions <- region_svs[
            region_svs$SVtype %in% c("h2hINV", "t2tINV"),
        ]

        tandem_dup_list[[i]] <- list(
            n_tandem_dups = nrow(tandem_dups),
            n_inversions = nrow(inversions),
            n_total_svs = nrow(region_svs),
            tandem_dups = tandem_dups,
            inversions = inversions,
            all_svs = region_svs
        )
    }

    return(tandem_dup_list)
}


#' Evaluate a chromosynthesis region
#'
#' @param region Gradient region
#' @param SV.sample SV data
#' @param CNV.sample CNV data
#' @param tandem_dups Tandem duplication analysis
#' @param min_tandem_dups Minimum tandem duplications
#' @return Evaluation results
#' @keywords internal
evaluate_chromosynthesis_region <- function(region,
                                           SV.sample,
                                           CNV.sample,
                                           tandem_dups,
                                           min_tandem_dups = 3) {

    # Calculate tandem duplication enrichment
    tandem_dup_enrichment <- if (tandem_dups$n_total_svs > 0) {
        tandem_dups$n_tandem_dups / tandem_dups$n_total_svs
    } else {
        0
    }

    # Calculate complexity score
    complexity_score <- calculate_chromosynthesis_complexity(
        region = region,
        tandem_dups = tandem_dups
    )

    # Check for inter-chromosomal events (should be minimal)
    inter_chr_svs <- SV.sample[
        SV.sample$chrom1 != SV.sample$chrom2 &
        ((SV.sample$chrom1 == region$chrom &
          SV.sample$pos1 >= region$start &
          SV.sample$pos1 <= region$end) |
         (SV.sample$chrom2 == region$chrom &
          SV.sample$pos2 >= region$start &
          SV.sample$pos2 <= region$end)),
    ]

    n_inter_chr <- nrow(inter_chr_svs)

    # Chromosynthesis should have minimal inter-chromosomal events
    inter_chr_penalty <- min(n_inter_chr / 5, 1.0)  # Penalty increases with more inter-chr

    # Create summary
    summary <- data.frame(
        region_id = if (!is.null(region$region_id)) region$region_id else 1,
        chrom = region$chrom,
        start = region$start,
        end = region$end,
        size = region$size,
        n_cn_segments = region$n_segments,
        gradient_correlation = region$gradient_correlation,
        cn_slope = region$slope,
        cn_range = region$cn_range,
        n_tandem_dups = tandem_dups$n_tandem_dups,
        n_inversions = tandem_dups$n_inversions,
        tandem_dup_enrichment = tandem_dup_enrichment,
        n_inter_chr_svs = n_inter_chr,
        complexity_score = complexity_score,
        stringsAsFactors = FALSE
    )

    return(list(
        summary = summary,
        region = region,
        tandem_dups = tandem_dups
    ))
}


#' Calculate chromosynthesis complexity score
#'
#' @param region Gradient region
#' @param tandem_dups Tandem duplication data
#' @return Complexity score (0-1)
#' @keywords internal
calculate_chromosynthesis_complexity <- function(region, tandem_dups) {

    # Complexity based on:
    # 1. Number of CN segments (more = more complex)
    # 2. Gradient strength
    # 3. Tandem duplication presence
    # 4. Inversion presence

    segment_score <- min(region$n_segments / 10, 1.0)
    gradient_score <- region$gradient_correlation
    tandem_score <- min(tandem_dups$n_tandem_dups / 5, 1.0)
    inversion_score <- min(tandem_dups$n_inversions / 3, 0.5)

    complexity <- (
        segment_score * 0.3 +
        gradient_score * 0.3 +
        tandem_score * 0.3 +
        inversion_score * 0.1
    )

    return(complexity)
}


#' Classify chromosynthesis event
#'
#' @param summary_df Summary data frame
#' @return Classification vector
#' @keywords internal
classify_chromosynthesis_event <- function(summary_df) {

    classifications <- character(nrow(summary_df))

    for (i in 1:nrow(summary_df)) {
        row <- summary_df[i, ]

        # Criteria for chromosynthesis:
        # 1. Strong CN gradient (correlation ≥ 0.6)
        # 2. Multiple tandem duplications (≥ 3)
        # 3. Complex local rearrangements (score ≥ 0.4)
        # 4. Minimal inter-chromosomal events (< 2)

        meets_gradient <- row$gradient_correlation >= 0.6
        meets_tandem <- row$n_tandem_dups >= 3
        meets_complexity <- row$complexity_score >= 0.4
        meets_locality <- row$n_inter_chr_svs < 2

        criteria_met <- sum(c(meets_gradient, meets_tandem,
                             meets_complexity, meets_locality))

        if (criteria_met >= 4) {
            classifications[i] <- "Likely chromosynthesis"
        } else if (criteria_met >= 3) {
            classifications[i] <- "Possible chromosynthesis"
        } else if (criteria_met >= 2) {
            classifications[i] <- "Unlikely chromosynthesis"
        } else {
            classifications[i] <- "Not chromosynthesis"
        }
    }

    return(classifications)
}


#' Create empty chromosynthesis result
#'
#' @return Empty result structure
#' @keywords internal
create_empty_chromosynthesis_result <- function() {
    result <- list(
        regions = list(),
        region_details = list(),
        summary = data.frame(
            region_id = integer(0),
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            size = numeric(0),
            n_cn_segments = integer(0),
            gradient_correlation = numeric(0),
            cn_slope = numeric(0),
            cn_range = numeric(0),
            n_tandem_dups = integer(0),
            n_inversions = integer(0),
            tandem_dup_enrichment = numeric(0),
            n_inter_chr_svs = integer(0),
            complexity_score = numeric(0),
            classification = character(0)
        ),
        total_regions = 0,
        likely_chromosynthesis = 0,
        possible_chromosynthesis = 0
    )

    class(result) <- c("chromosynthesis", "list")
    return(result)
}


#' Print method for chromosynthesis results
#'
#' @param x Chromosynthesis result object
#' @param ... Additional arguments
#' @export
print.chromosynthesis <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOSYNTHESIS DETECTION RESULTS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Total regions detected: %d\n", x$total_regions))
    cat(sprintf("  - Likely chromosynthesis:   %d\n", x$likely_chromosynthesis))
    cat(sprintf("  - Possible chromosynthesis: %d\n", x$possible_chromosynthesis))
    cat("\n")

    if (x$total_regions > 0) {
        cat("Region Summary:\n")
        cat(rep("-", 70), "\n", sep = "")
        print(x$summary[, c("region_id", "chrom", "size", "gradient_correlation",
                           "n_tandem_dups", "classification")])
    } else {
        cat("No chromosynthesis regions detected.\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for chromosynthesis results
#'
#' @param object Chromosynthesis result object
#' @param ... Additional arguments
#' @export
summary.chromosynthesis <- function(object, ...) {
    print(object)
    invisible(object)
}

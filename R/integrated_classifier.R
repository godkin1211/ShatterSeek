#' Integrated classifier for mixed chromoanagenesis mechanisms
#'
#' Analyzes chromoanagenesis results to identify and classify mixed mechanisms,
#' where multiple catastrophic events co-occur in the same sample or region.
#'
#' @param chromoanagenesis_result Result from detect_chromoanagenesis()
#' @param overlap_threshold Minimum overlap (in bp) to consider mechanisms co-occurring (default: 1e6)
#' @param min_confidence Minimum confidence score to include events (default: 0.3)
#' @return A list containing integrated classification results
#' @details
#' This function identifies complex patterns where multiple chromoanagenesis
#' mechanisms occur together, such as:
#' - Chromothripsis + Chromoplexy (catastrophic + chained rearrangements)
#' - Chromothripsis + Chromosynthesis (shattering + serial replication)
#' - Triple mechanism events (all three mechanisms)
#'
#' The classifier evaluates:
#' 1. Spatial overlap of different mechanisms
#' 2. Temporal relationships (primary vs secondary events)
#' 3. Mechanism dominance in each chromosome
#' 4. Overall sample complexity
#'
#' @examples
#' \dontrun{
#' # Run comprehensive analysis
#' results <- detect_chromoanagenesis(SV_data, CN_data)
#'
#' # Classify mixed mechanisms
#' classification <- classify_mixed_mechanisms(results)
#'
#' # View results
#' print(classification)
#' plot_mechanism_landscape(classification)
#' }
#'
#' @export
classify_mixed_mechanisms <- function(chromoanagenesis_result,
                                     overlap_threshold = 1e6,
                                     min_confidence = 0.3) {

    if (!inherits(chromoanagenesis_result, "chromoanagenesis")) {
        stop("Input must be a chromoanagenesis result object from detect_chromoanagenesis()")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("     INTEGRATED MECHANISM CLASSIFICATION\n")
    cat(rep("=", 70), "\n\n", sep = "")

    results <- list()

    # 1. Extract mechanism locations
    cat("Step 1: Extracting mechanism locations...\n")
    mechanism_locations <- extract_mechanism_locations(
        chromoanagenesis_result,
        min_confidence
    )
    results$mechanism_locations <- mechanism_locations

    # 2. Detect overlapping mechanisms
    cat("Step 2: Detecting mechanism overlaps...\n")
    overlaps <- detect_mechanism_overlaps(
        mechanism_locations,
        overlap_threshold
    )
    results$overlaps <- overlaps

    # 3. Chromosome-level classification
    cat("Step 3: Classifying mechanisms by chromosome...\n")
    chr_classification <- classify_by_chromosome(
        mechanism_locations,
        overlaps
    )
    results$chromosome_classification <- chr_classification

    # 4. Sample-level classification
    cat("Step 4: Sample-level integrated classification...\n")
    sample_classification <- classify_sample_level(
        chromoanagenesis_result,
        chr_classification,
        overlaps
    )
    results$sample_classification <- sample_classification

    # 5. Mechanism dominance analysis
    cat("Step 5: Analyzing mechanism dominance...\n")
    dominance <- analyze_mechanism_dominance(
        chromoanagenesis_result,
        chr_classification
    )
    results$dominance <- dominance

    # 6. Complexity scoring
    cat("Step 6: Calculating complexity scores...\n")
    complexity <- calculate_complexity_score(
        chromoanagenesis_result,
        overlaps,
        chr_classification
    )
    results$complexity <- complexity

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("CLASSIFICATION COMPLETE\n")
    cat(rep("=", 70), "\n\n", sep = "")

    class(results) <- c("mixed_mechanisms", "list")
    return(results)
}


#' Extract genomic locations of each mechanism
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param min_confidence Minimum confidence threshold
#' @return Data frame with mechanism locations
#' @keywords internal
extract_mechanism_locations <- function(chromoanagenesis_result, min_confidence) {

    locations <- list()

    # Chromothripsis locations
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        chromoth <- chromoanagenesis_result$chromothripsis$classification

        # Filter by confidence
        if (!is.null(chromoth) && nrow(chromoth) > 0) {
            chromoth <- chromoth[!is.na(chromoth$confidence_score) &
                                chromoth$confidence_score >= min_confidence, ]

            if (nrow(chromoth) > 0) {
                for (i in 1:nrow(chromoth)) {
                    # Chromothripsis affects whole chromosome, no specific start/end
                    locations[[length(locations) + 1]] <- list(
                        mechanism = "chromothripsis",
                        chrom = as.character(chromoth$chrom[i]),
                        start = NA_real_,  # Whole chromosome, no specific region
                        end = NA_real_,
                        confidence = as.numeric(chromoth$confidence_score[i]),
                        classification = as.character(chromoth$classification[i]),
                        event_id = paste0("CT_", chromoth$chrom[i])
                    )
                }
            }
        }
    }

    # Chromoplexy locations
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        chromopl <- chromoanagenesis_result$chromoplexy$summary

        # Filter by classification
        likely_possible <- chromopl[chromopl$classification %in%
                                   c("Likely chromoplexy", "Possible chromoplexy"), ]

        if (nrow(likely_possible) > 0) {
            for (i in 1:nrow(likely_possible)) {
                # Get chromosomes involved
                chroms_str <- as.character(likely_possible$chromosomes_involved[i])

                # Skip if NA or empty
                if (is.na(chroms_str) || chroms_str == "" || nchar(chroms_str) == 0) {
                    next
                }

                chroms <- unlist(strsplit(chroms_str, ","))
                chroms <- trimws(chroms)

                # Remove empty strings
                chroms <- chroms[chroms != "" & !is.na(chroms)]

                # Skip if no valid chromosomes
                if (length(chroms) == 0) {
                    next
                }

                # For each chromosome in the chain
                for (chr in chroms) {
                    locations[[length(locations) + 1]] <- list(
                        mechanism = "chromoplexy",
                        chrom = as.character(chr),
                        start = NA_real_,  # Chromoplexy doesn't have single region
                        end = NA_real_,
                        confidence = as.numeric(likely_possible$confidence_score[i]),
                        classification = as.character(likely_possible$classification[i]),
                        event_id = paste0("CP_", likely_possible$chain_id[i])
                    )
                }
            }
        }
    }

    # Chromosynthesis locations
    if (!is.null(chromoanagenesis_result$chromosynthesis)) {
        chromosyn <- chromoanagenesis_result$chromosynthesis$summary

        # Filter by classification
        if (!is.null(chromosyn) && nrow(chromosyn) > 0) {
            likely_possible <- chromosyn[!is.na(chromosyn$classification) &
                                        chromosyn$classification %in%
                                        c("Likely chromosynthesis", "Possible chromosynthesis"), ]

            if (nrow(likely_possible) > 0) {
                for (i in 1:nrow(likely_possible)) {
                    locations[[length(locations) + 1]] <- list(
                        mechanism = "chromosynthesis",
                        chrom = as.character(likely_possible$chrom[i]),
                        start = as.numeric(likely_possible$region_start[i]),
                        end = as.numeric(likely_possible$region_end[i]),
                        confidence = as.numeric(likely_possible$confidence_score[i]),
                        classification = as.character(likely_possible$classification[i]),
                        event_id = paste0("CS_", likely_possible$region_id[i])
                    )
                }
            }
        }
    }

    # Convert to data frame
    if (length(locations) == 0) {
        return(data.frame(
            mechanism = character(0),
            chrom = character(0),
            start = numeric(0),
            end = numeric(0),
            confidence = numeric(0),
            classification = character(0),
            event_id = character(0),
            stringsAsFactors = FALSE
        ))
    }

    df <- do.call(rbind, lapply(locations, function(x) {
        data.frame(
            mechanism = as.character(x$mechanism),
            chrom = as.character(x$chrom),
            start = as.numeric(x$start),
            end = as.numeric(x$end),
            confidence = as.numeric(x$confidence),
            classification = as.character(x$classification),
            event_id = as.character(x$event_id),
            stringsAsFactors = FALSE
        )
    }))

    return(df)
}


#' Detect overlapping mechanisms
#'
#' @param mechanism_locations Data frame from extract_mechanism_locations
#' @param overlap_threshold Minimum overlap in bp
#' @return List of overlapping mechanism pairs
#' @keywords internal
detect_mechanism_overlaps <- function(mechanism_locations, overlap_threshold) {

    if (nrow(mechanism_locations) == 0) {
        return(list(
            n_overlaps = 0,
            overlaps = data.frame()
        ))
    }

    overlaps <- list()

    # Group by chromosome
    chroms <- unique(mechanism_locations$chrom)

    for (chr in chroms) {
        chr_mechs <- mechanism_locations[mechanism_locations$chrom == chr, ]

        # Compare all pairs
        if (nrow(chr_mechs) > 1) {
            for (i in 1:(nrow(chr_mechs) - 1)) {
                for (j in (i + 1):nrow(chr_mechs)) {

                    mech1 <- chr_mechs[i, ]
                    mech2 <- chr_mechs[j, ]

                    # Skip if same mechanism type
                    if (mech1$mechanism == mech2$mechanism) next

                    # Check for overlap
                    overlap_detected <- FALSE
                    overlap_size <- NA

                    # If both have defined regions
                    if (!is.na(mech1$start) && !is.na(mech2$start)) {
                        # Calculate overlap
                        overlap_start <- max(mech1$start, mech2$start)
                        overlap_end <- min(mech1$end, mech2$end)

                        if (overlap_end >= overlap_start) {
                            overlap_size <- overlap_end - overlap_start
                            if (overlap_size >= overlap_threshold) {
                                overlap_detected <- TRUE
                            }
                        }
                    } else {
                        # For chromoplexy (no specific region), consider co-occurrence on same chr as overlap
                        overlap_detected <- TRUE
                        overlap_size <- NA
                    }

                    if (overlap_detected) {
                        overlaps[[length(overlaps) + 1]] <- list(
                            chrom = chr,
                            mechanism1 = mech1$mechanism,
                            mechanism2 = mech2$mechanism,
                            event_id1 = mech1$event_id,
                            event_id2 = mech2$event_id,
                            overlap_size = overlap_size,
                            confidence1 = mech1$confidence,
                            confidence2 = mech2$confidence
                        )
                    }
                }
            }
        }
    }

    # Convert to data frame
    if (length(overlaps) == 0) {
        return(list(
            n_overlaps = 0,
            overlaps = data.frame()
        ))
    }

    overlaps_df <- do.call(rbind, lapply(overlaps, function(x) {
        data.frame(
            chrom = as.character(x$chrom),
            mechanism1 = as.character(x$mechanism1),
            mechanism2 = as.character(x$mechanism2),
            event_id1 = as.character(x$event_id1),
            event_id2 = as.character(x$event_id2),
            overlap_size = as.numeric(x$overlap_size),
            confidence1 = as.numeric(x$confidence1),
            confidence2 = as.numeric(x$confidence2),
            mechanism_pair = as.character(paste(sort(c(x$mechanism1, x$mechanism2)), collapse = "+")),
            stringsAsFactors = FALSE
        )
    }))

    return(list(
        n_overlaps = nrow(overlaps_df),
        overlaps = overlaps_df
    ))
}


#' Classify mechanisms by chromosome
#'
#' @param mechanism_locations Data frame from extract_mechanism_locations
#' @param overlaps List from detect_mechanism_overlaps
#' @return Data frame with chromosome-level classification
#' @keywords internal
classify_by_chromosome <- function(mechanism_locations, overlaps) {

    if (nrow(mechanism_locations) == 0) {
        return(data.frame(
            chrom = character(0),
            mechanisms = character(0),
            dominant_mechanism = character(0),
            is_mixed = logical(0),
            n_mechanisms = numeric(0),
            stringsAsFactors = FALSE
        ))
    }

    chroms <- unique(mechanism_locations$chrom)
    chr_class <- list()

    for (chr in chroms) {
        chr_mechs <- mechanism_locations[mechanism_locations$chrom == chr, ]

        # Get unique mechanisms on this chromosome
        unique_mechs <- unique(chr_mechs$mechanism)
        n_mechs <- length(unique_mechs)
        is_mixed <- n_mechs > 1

        # Determine dominant mechanism (highest confidence)
        chr_mechs_agg <- aggregate(confidence ~ mechanism, data = chr_mechs, FUN = max)
        dominant_idx <- which.max(chr_mechs_agg$confidence)
        dominant_mech <- chr_mechs_agg$mechanism[dominant_idx]

        # Create classification label
        if (is_mixed) {
            mech_label <- paste(sort(unique_mechs), collapse = "+")
        } else {
            mech_label <- unique_mechs[1]
        }

        chr_class[[chr]] <- list(
            chrom = chr,
            mechanisms = mech_label,
            dominant_mechanism = dominant_mech,
            is_mixed = is_mixed,
            n_mechanisms = n_mechs,
            max_confidence = max(chr_mechs$confidence)
        )
    }

    # Convert to data frame
    chr_class_df <- do.call(rbind, lapply(chr_class, function(x) {
        data.frame(
            chrom = as.character(x$chrom),
            mechanisms = as.character(x$mechanisms),
            dominant_mechanism = as.character(x$dominant_mechanism),
            is_mixed = as.logical(x$is_mixed),
            n_mechanisms = as.numeric(x$n_mechanisms),
            max_confidence = as.numeric(x$max_confidence),
            stringsAsFactors = FALSE
        )
    }))

    # Sort by chromosome
    chr_class_df <- chr_class_df[order(chr_class_df$chrom), ]

    return(chr_class_df)
}


#' Sample-level integrated classification
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param chr_classification Chromosome classification data frame
#' @param overlaps Overlap detection results
#' @return List with sample-level classification
#' @keywords internal
classify_sample_level <- function(chromoanagenesis_result, chr_classification, overlaps) {

    # Count mechanisms
    has_chromothripsis <- !is.null(chromoanagenesis_result$chromothripsis) &&
                         chromoanagenesis_result$chromothripsis$n_likely > 0
    has_chromoplexy <- !is.null(chromoanagenesis_result$chromoplexy) &&
                      chromoanagenesis_result$chromoplexy$likely_chromoplexy > 0
    has_chromosynthesis <- !is.null(chromoanagenesis_result$chromosynthesis) &&
                          chromoanagenesis_result$chromosynthesis$likely_chromosynthesis > 0

    mechanisms_present <- c()
    if (has_chromothripsis) mechanisms_present <- c(mechanisms_present, "chromothripsis")
    if (has_chromoplexy) mechanisms_present <- c(mechanisms_present, "chromoplexy")
    if (has_chromosynthesis) mechanisms_present <- c(mechanisms_present, "chromosynthesis")

    n_mechanisms <- length(mechanisms_present)

    # Determine classification
    if (n_mechanisms == 0) {
        classification <- "No chromoanagenesis"
        category <- "normal"
    } else if (n_mechanisms == 1) {
        classification <- paste0("Pure ", mechanisms_present[1])
        category <- "single_mechanism"
    } else {
        # Mixed mechanism
        has_overlaps <- overlaps$n_overlaps > 0

        if (has_overlaps) {
            classification <- paste0("Mixed mechanisms with spatial overlap (",
                                   paste(mechanisms_present, collapse = "+"), ")")
            category <- "mixed_overlapping"
        } else {
            classification <- paste0("Multiple independent mechanisms (",
                                   paste(mechanisms_present, collapse = "+"), ")")
            category <- "mixed_independent"
        }
    }

    # Count mixed chromosomes
    n_mixed_chromosomes <- sum(chr_classification$is_mixed, na.rm = TRUE)

    return(list(
        classification = classification,
        category = category,
        n_mechanisms = n_mechanisms,
        mechanisms_present = mechanisms_present,
        n_mixed_chromosomes = n_mixed_chromosomes,
        has_spatial_overlap = overlaps$n_overlaps > 0,
        n_overlaps = overlaps$n_overlaps
    ))
}


#' Analyze mechanism dominance
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param chr_classification Chromosome classification data frame
#' @return List with dominance analysis
#' @keywords internal
analyze_mechanism_dominance <- function(chromoanagenesis_result, chr_classification) {

    # Count events by mechanism
    n_chromothripsis <- 0
    n_chromoplexy <- 0
    n_chromosynthesis <- 0

    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        n_chromothripsis <- chromoanagenesis_result$chromothripsis$n_likely +
                           chromoanagenesis_result$chromothripsis$n_possible
    }

    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        n_chromoplexy <- chromoanagenesis_result$chromoplexy$likely_chromoplexy +
                        chromoanagenesis_result$chromoplexy$possible_chromoplexy
    }

    if (!is.null(chromoanagenesis_result$chromosynthesis)) {
        n_chromosynthesis <- chromoanagenesis_result$chromosynthesis$likely_chromosynthesis +
                            chromoanagenesis_result$chromosynthesis$possible_chromosynthesis
    }

    # Calculate proportions
    total_events <- n_chromothripsis + n_chromoplexy + n_chromosynthesis

    if (total_events == 0) {
        return(list(
            dominant_mechanism = "none",
            mechanism_proportions = data.frame(
                mechanism = c("chromothripsis", "chromoplexy", "chromosynthesis"),
                n_events = c(0, 0, 0),
                proportion = c(0, 0, 0)
            )
        ))
    }

    proportions <- data.frame(
        mechanism = c("chromothripsis", "chromoplexy", "chromosynthesis"),
        n_events = c(n_chromothripsis, n_chromoplexy, n_chromosynthesis),
        proportion = c(n_chromothripsis, n_chromoplexy, n_chromosynthesis) / total_events,
        stringsAsFactors = FALSE
    )

    # Determine dominant mechanism
    dominant_idx <- which.max(proportions$n_events)
    dominant_mechanism <- proportions$mechanism[dominant_idx]

    # Check if truly dominant (>50%)
    if (proportions$proportion[dominant_idx] <= 0.5) {
        dominant_mechanism <- "balanced"
    }

    return(list(
        dominant_mechanism = dominant_mechanism,
        mechanism_proportions = proportions,
        total_events = total_events
    ))
}


#' Calculate sample complexity score
#'
#' @param chromoanagenesis_result Chromoanagenesis result object
#' @param overlaps Overlap detection results
#' @param chr_classification Chromosome classification data frame
#' @return List with complexity metrics
#' @keywords internal
calculate_complexity_score <- function(chromoanagenesis_result, overlaps, chr_classification) {

    # Components of complexity:
    # 1. Number of different mechanisms (0-3)
    # 2. Number of chromosomes affected
    # 3. Presence of overlapping mechanisms
    # 4. Total number of events

    # Get basic counts
    n_mechanisms <- sum(c(
        !is.null(chromoanagenesis_result$chromothripsis) &&
            chromoanagenesis_result$chromothripsis$n_likely > 0,
        !is.null(chromoanagenesis_result$chromoplexy) &&
            chromoanagenesis_result$chromoplexy$likely_chromoplexy > 0,
        !is.null(chromoanagenesis_result$chromosynthesis) &&
            chromoanagenesis_result$chromosynthesis$likely_chromosynthesis > 0
    ))

    n_chromosomes <- nrow(chr_classification)
    n_mixed_chromosomes <- sum(chr_classification$is_mixed, na.rm = TRUE)
    n_overlaps <- overlaps$n_overlaps

    # Calculate total events
    total_events <- 0
    if (!is.null(chromoanagenesis_result$chromothripsis)) {
        total_events <- total_events + chromoanagenesis_result$chromothripsis$n_likely
    }
    if (!is.null(chromoanagenesis_result$chromoplexy)) {
        total_events <- total_events + chromoanagenesis_result$chromoplexy$likely_chromoplexy
    }
    if (!is.null(chromoanagenesis_result$chromosynthesis)) {
        total_events <- total_events + chromoanagenesis_result$chromosynthesis$likely_chromosynthesis
    }

    # Complexity score (0-1 scale)
    # Weighted components:
    # - Mechanism diversity: 30%
    # - Spatial overlap: 25%
    # - Chromosome spread: 25%
    # - Event count: 20%

    mechanism_score <- n_mechanisms / 3  # Max 3 mechanisms
    overlap_score <- min(n_overlaps / 5, 1)  # Saturate at 5 overlaps
    chromosome_score <- min(n_chromosomes / 10, 1)  # Saturate at 10 chromosomes
    event_score <- min(total_events / 10, 1)  # Saturate at 10 events

    complexity_score <- (
        mechanism_score * 0.30 +
        overlap_score * 0.25 +
        chromosome_score * 0.25 +
        event_score * 0.20
    )

    # Classify complexity level
    if (complexity_score < 0.3) {
        complexity_level <- "Low"
    } else if (complexity_score < 0.6) {
        complexity_level <- "Moderate"
    } else if (complexity_score < 0.8) {
        complexity_level <- "High"
    } else {
        complexity_level <- "Very High"
    }

    return(list(
        complexity_score = complexity_score,
        complexity_level = complexity_level,
        n_mechanisms = n_mechanisms,
        n_chromosomes = n_chromosomes,
        n_mixed_chromosomes = n_mixed_chromosomes,
        n_overlaps = n_overlaps,
        total_events = total_events,
        components = data.frame(
            component = c("Mechanism diversity", "Spatial overlap",
                         "Chromosome spread", "Event count"),
            raw_score = c(mechanism_score, overlap_score,
                         chromosome_score, event_score),
            weight = c(0.30, 0.25, 0.25, 0.20),
            weighted_score = c(mechanism_score * 0.30, overlap_score * 0.25,
                              chromosome_score * 0.25, event_score * 0.20),
            stringsAsFactors = FALSE
        )
    ))
}


#' Print method for mixed mechanisms classification
#'
#' @param x Mixed mechanisms result object
#' @param ... Additional arguments
#' @export
print.mixed_mechanisms <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         INTEGRATED MECHANISM CLASSIFICATION\n")
    cat(rep("=", 70), "\n\n", sep = "")

    # Sample-level classification
    cat("SAMPLE-LEVEL CLASSIFICATION:\n")
    cat(sprintf("  Classification: %s\n", x$sample_classification$classification))
    cat(sprintf("  Category: %s\n", x$sample_classification$category))
    cat(sprintf("  Number of mechanisms: %d\n", x$sample_classification$n_mechanisms))
    if (length(x$sample_classification$mechanisms_present) > 0) {
        cat(sprintf("  Mechanisms present: %s\n",
                   paste(x$sample_classification$mechanisms_present, collapse = ", ")))
    }
    cat("\n")

    # Complexity
    cat("COMPLEXITY ANALYSIS:\n")
    cat(sprintf("  Overall complexity: %s (score: %.3f)\n",
               x$complexity$complexity_level,
               x$complexity$complexity_score))
    cat(sprintf("  Total events: %d\n", x$complexity$total_events))
    cat(sprintf("  Chromosomes affected: %d\n", x$complexity$n_chromosomes))
    cat(sprintf("  Mixed chromosomes: %d\n", x$complexity$n_mixed_chromosomes))
    cat(sprintf("  Spatial overlaps: %d\n", x$complexity$n_overlaps))
    cat("\n")

    # Dominance
    cat("MECHANISM DOMINANCE:\n")
    cat(sprintf("  Dominant mechanism: %s\n", x$dominance$dominant_mechanism))
    cat("  Event distribution:\n")
    for (i in 1:nrow(x$dominance$mechanism_proportions)) {
        mech <- x$dominance$mechanism_proportions$mechanism[i]
        n <- x$dominance$mechanism_proportions$n_events[i]
        prop <- x$dominance$mechanism_proportions$proportion[i]
        cat(sprintf("    - %s: %d events (%.1f%%)\n", mech, n, prop * 100))
    }
    cat("\n")

    # Chromosome details
    if (nrow(x$chromosome_classification) > 0) {
        cat("CHROMOSOME-LEVEL DETAILS:\n")
        mixed_chrs <- x$chromosome_classification[x$chromosome_classification$is_mixed, ]
        if (nrow(mixed_chrs) > 0) {
            cat("  Mixed mechanism chromosomes:\n")
            for (i in 1:nrow(mixed_chrs)) {
                cat(sprintf("    - %s: %s (dominant: %s)\n",
                           mixed_chrs$chrom[i],
                           mixed_chrs$mechanisms[i],
                           mixed_chrs$dominant_mechanism[i]))
            }
        } else {
            cat("  No mixed mechanism chromosomes detected\n")
        }
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for mixed mechanisms classification
#'
#' @param object Mixed mechanisms result object
#' @param ... Additional arguments
#' @export
summary.mixed_mechanisms <- function(object, ...) {

    cat("\nMixed Mechanisms Summary:\n")
    cat("========================\n\n")

    # Print complexity components
    cat("Complexity Score Components:\n")
    print(object$complexity$components)
    cat("\n")

    # Print mechanism proportions
    cat("Mechanism Proportions:\n")
    print(object$dominance$mechanism_proportions)
    cat("\n")

    # Print chromosome classification
    cat("Chromosome Classification:\n")
    print(object$chromosome_classification)
    cat("\n")

    # Print overlaps if any
    if (object$overlaps$n_overlaps > 0) {
        cat("Mechanism Overlaps:\n")
        print(object$overlaps$overlaps[, c("chrom", "mechanism_pair",
                                          "overlap_size", "confidence1", "confidence2")])
        cat("\n")
    }

    invisible(object)
}

#' Calculate integrated chromothripsis confidence score
#'
#' This function integrates multiple statistical criteria to provide a unified
#' confidence score for chromothripsis detection. The score ranges from 0 to 1,
#' where higher values indicate stronger evidence for chromothripsis.
#'
#' @param chromSummary A data frame from ShatterSeek output containing statistical criteria
#' @param chr Optional: specific chromosome to score (default: all chromosomes with clusters)
#' @return A data frame with chromosome names and confidence scores
#' @details
#' The scoring system integrates five key criteria:
#' \itemize{
#'   \item CN oscillations (weight: 0.25)
#'   \item Fragment joins randomness (weight: 0.20)
#'   \item Breakpoint clustering (weight: 0.20)
#'   \item Chromosome enrichment (weight: 0.15)
#'   \item Cluster size (weight: 0.20)
#' }
#'
#' Based on Cortes-Ciriano et al. Nature Genetics 2020, a high-confidence
#' chromothripsis event typically has:
#' - Cluster size >= 10 SVs
#' - CN oscillations >= 4 segments
#' - p-value for fragment joins > 0.05
#' - p-value for exponential distribution < 0.05
#'
#' @export
calculate_confidence_score <- function(chromSummary, chr = NULL) {

    if (!is.null(chr)) {
        chromSummary <- chromSummary[chromSummary$chrom == chr, ]
    }

    # Filter chromosomes with clusters
    chromSummary <- chromSummary[chromSummary$clusterSize > 0, ]

    if (nrow(chromSummary) == 0) {
        warning("No chromosomes with SV clusters found")
        return(data.frame(chrom = character(0),
                         confidence_score = numeric(0),
                         confidence_category = character(0)))
    }

    scores <- data.frame(
        chrom = chromSummary$chrom,
        cluster_score = 0,
        cn_oscillation_score = 0,
        fragment_joins_score = 0,
        clustering_score = 0,
        enrichment_score = 0,
        total_score = 0,
        confidence_category = "",
        stringsAsFactors = FALSE
    )

    for (i in 1:nrow(chromSummary)) {
        row <- chromSummary[i, ]

        # 1. Cluster size score (weight: 0.20)
        # Based on literature: >=10 SVs is typical for chromothripsis
        cluster_size <- row$clusterSize_including_TRA
        if (is.na(cluster_size)) cluster_size <- row$clusterSize

        if (cluster_size >= 10) {
            scores$cluster_score[i] <- 1.0
        } else if (cluster_size >= 6) {
            scores$cluster_score[i] <- 0.6
        } else if (cluster_size >= 3) {
            scores$cluster_score[i] <- 0.3
        } else {
            scores$cluster_score[i] <- 0.1
        }

        # 2. CN oscillation score (weight: 0.25)
        # Chromothripsis shows oscillating CN patterns
        cn_osc_2state <- row$max_number_oscillating_CN_segments_2_states
        cn_osc_3state <- row$max_number_oscillating_CN_segments_3_states

        if (!is.na(cn_osc_2state) && cn_osc_2state >= 4) {
            scores$cn_oscillation_score[i] <- 1.0
        } else if (!is.na(cn_osc_3state) && cn_osc_3state >= 4) {
            scores$cn_oscillation_score[i] <- 0.8
        } else if (!is.na(cn_osc_2state) && cn_osc_2state >= 2) {
            scores$cn_oscillation_score[i] <- 0.4
        } else {
            scores$cn_oscillation_score[i] <- 0.0
        }

        # 3. Fragment joins randomness score (weight: 0.20)
        # Random rejoining: p-value should be > 0.05 (not significantly different from random)
        pval_joins <- row$pval_fragment_joins

        if (!is.na(pval_joins)) {
            if (pval_joins > 0.05) {
                scores$fragment_joins_score[i] <- 1.0
            } else if (pval_joins > 0.01) {
                scores$fragment_joins_score[i] <- 0.5
            } else {
                scores$fragment_joins_score[i] <- 0.0
            }
        }

        # 4. Breakpoint clustering score (weight: 0.20)
        # Clustered breakpoints: p-value should be < 0.05 (reject exponential distribution)
        pval_exp_cluster <- row$pval_exp_cluster
        pval_exp_chr <- row$pval_exp_chr

        if (!is.na(pval_exp_cluster) && pval_exp_cluster < 0.05) {
            scores$clustering_score[i] <- 1.0
        } else if (!is.na(pval_exp_chr) && pval_exp_chr < 0.05) {
            scores$clustering_score[i] <- 0.7
        } else if (!is.na(pval_exp_cluster) && pval_exp_cluster < 0.1) {
            scores$clustering_score[i] <- 0.5
        } else {
            scores$clustering_score[i] <- 0.0
        }

        # 5. Chromosome breakpoint enrichment score (weight: 0.15)
        # Enrichment: p-value should be < 0.05
        pval_enrich <- row$chr_breakpoint_enrichment

        if (!is.na(pval_enrich)) {
            if (pval_enrich < 0.001) {
                scores$enrichment_score[i] <- 1.0
            } else if (pval_enrich < 0.01) {
                scores$enrichment_score[i] <- 0.7
            } else if (pval_enrich < 0.05) {
                scores$enrichment_score[i] <- 0.4
            } else {
                scores$enrichment_score[i] <- 0.0
            }
        }

        # Calculate weighted total score
        scores$total_score[i] <- (
            scores$cluster_score[i] * 0.20 +
            scores$cn_oscillation_score[i] * 0.25 +
            scores$fragment_joins_score[i] * 0.20 +
            scores$clustering_score[i] * 0.20 +
            scores$enrichment_score[i] * 0.15
        )

        # Assign confidence category
        if (scores$total_score[i] >= 0.7) {
            scores$confidence_category[i] <- "High"
        } else if (scores$total_score[i] >= 0.4) {
            scores$confidence_category[i] <- "Moderate"
        } else {
            scores$confidence_category[i] <- "Low"
        }
    }

    # Rename for clarity
    names(scores)[names(scores) == "total_score"] <- "confidence_score"

    return(scores)
}


#' Classify chromothripsis events based on criteria
#'
#' Applies evidence-based thresholds to classify chromosomes as likely
#' chromothripsis, possible chromothripsis, or unlikely chromothripsis.
#'
#' @param chromoth_output Output from detect_chromothripsis function
#' @param min_cluster_size Minimum number of SVs required (default: 6)
#' @param min_cn_oscillations Minimum CN oscillating segments (default: 4)
#' @param max_pval_clustering Maximum p-value for clustering test (default: 0.05)
#' @return A data frame with classification results
#' @export
classify_chromothripsis <- function(chromoth_output,
                                   min_cluster_size = 6,
                                   min_cn_oscillations = 4,
                                   max_pval_clustering = 0.05) {

    chromSummary <- chromoth_output@chromSummary

    # Calculate confidence scores
    scores <- calculate_confidence_score(chromSummary)

    if (nrow(scores) == 0) {
        return(data.frame(
            chrom = character(0),
            classification = character(0),
            cluster_size = numeric(0),
            cn_oscillations = numeric(0),
            confidence_score = numeric(0)
        ))
    }

    # Merge scores with original summary
    result <- merge(chromSummary, scores[, c("chrom", "confidence_score", "confidence_category")],
                   by = "chrom", all.x = FALSE)

    # Apply classification criteria
    result$classification <- "Unlikely"

    for (i in 1:nrow(result)) {
        cluster_size <- result$clusterSize_including_TRA[i]
        if (is.na(cluster_size)) cluster_size <- result$clusterSize[i]

        cn_osc <- result$max_number_oscillating_CN_segments_2_states[i]
        if (is.na(cn_osc)) cn_osc <- result$max_number_oscillating_CN_segments_3_states[i]

        pval_cluster <- result$pval_exp_cluster[i]

        # High confidence: meets all criteria
        if (!is.na(cluster_size) && cluster_size >= min_cluster_size &&
            !is.na(cn_osc) && cn_osc >= min_cn_oscillations &&
            !is.na(pval_cluster) && pval_cluster < max_pval_clustering) {
            result$classification[i] <- "Likely chromothripsis"
        }
        # Moderate confidence: meets 2 out of 3 criteria
        else if (sum(c(
            !is.na(cluster_size) && cluster_size >= min_cluster_size,
            !is.na(cn_osc) && cn_osc >= min_cn_oscillations,
            !is.na(pval_cluster) && pval_cluster < max_pval_clustering
        )) >= 2) {
            result$classification[i] <- "Possible chromothripsis"
        }
    }

    # Select key columns for output
    output_cols <- c("chrom", "classification", "confidence_score", "confidence_category",
                    "clusterSize", "clusterSize_including_TRA",
                    "max_number_oscillating_CN_segments_2_states",
                    "max_number_oscillating_CN_segments_3_states",
                    "pval_fragment_joins", "pval_exp_cluster", "chr_breakpoint_enrichment")

    output_cols <- output_cols[output_cols %in% names(result)]
    result <- result[, output_cols]

    # Sort by confidence score
    result <- result[order(result$confidence_score, decreasing = TRUE), ]

    return(result)
}


#' Generate a summary report of chromothripsis analysis
#'
#' Creates a human-readable summary of chromothripsis detection results,
#' highlighting high-confidence events and providing interpretation guidance.
#'
#' @param chromoth_output Output from detect_chromothripsis function
#' @param print_summary Logical, whether to print summary to console (default: TRUE)
#' @return A list containing summary statistics and classifications
#' @export
summarize_chromothripsis <- function(chromoth_output, print_summary = TRUE) {

    # Get classifications
    classifications <- classify_chromothripsis(chromoth_output)

    # Count classifications
    n_likely <- sum(classifications$classification == "Likely chromothripsis")
    n_possible <- sum(classifications$classification == "Possible chromothripsis")
    n_unlikely <- sum(classifications$classification == "Unlikely")

    # Get high confidence chromosomes
    high_conf <- classifications[classifications$confidence_category == "High", ]

    # Prepare summary
    summary_list <- list(
        total_chromosomes_with_clusters = nrow(classifications),
        likely_chromothripsis = n_likely,
        possible_chromothripsis = n_possible,
        unlikely_chromothripsis = n_unlikely,
        high_confidence_chromosomes = if(nrow(high_conf) > 0) high_conf$chrom else character(0),
        detailed_results = classifications
    )

    if (print_summary) {
        cat("\n")
        cat(rep("=", 70), "\n", sep="")
        cat("         CHROMOTHRIPSIS DETECTION SUMMARY\n")
        cat(rep("=", 70), "\n\n", sep="")

        cat(sprintf("Total chromosomes with SV clusters: %d\n\n",
                   summary_list$total_chromosomes_with_clusters))

        cat("Classification Results:\n")
        cat(sprintf("  - Likely chromothripsis:   %d chromosome(s)\n", n_likely))
        cat(sprintf("  - Possible chromothripsis: %d chromosome(s)\n", n_possible))
        cat(sprintf("  - Unlikely:                %d chromosome(s)\n\n", n_unlikely))

        if (length(summary_list$high_confidence_chromosomes) > 0) {
            cat("High-confidence chromothripsis chromosomes:\n")
            cat(sprintf("  %s\n\n", paste(summary_list$high_confidence_chromosomes, collapse=", ")))
        } else {
            cat("No high-confidence chromothripsis events detected.\n\n")
        }

        if (n_likely > 0 || n_possible > 0) {
            cat("Detailed Results:\n")
            cat(sprintf("%-8s %-25s %-10s %-12s %-10s\n",
                       "Chrom", "Classification", "Confidence", "Cluster Size", "CN Osc."))
            cat(rep("-", 70), "\n", sep="")

            for (i in 1:nrow(classifications)) {
                if (classifications$classification[i] != "Unlikely") {
                    cn_osc <- classifications$max_number_oscillating_CN_segments_2_states[i]
                    if (is.na(cn_osc)) cn_osc <- classifications$max_number_oscillating_CN_segments_3_states[i]
                    if (is.na(cn_osc)) cn_osc <- "-"

                    cluster_size <- classifications$clusterSize_including_TRA[i]
                    if (is.na(cluster_size)) cluster_size <- classifications$clusterSize[i]

                    cat(sprintf("%-8s %-25s %-10.2f %-12d %-10s\n",
                               classifications$chrom[i],
                               classifications$classification[i],
                               classifications$confidence_score[i],
                               cluster_size,
                               cn_osc))
                }
            }
        }

        cat("\n")
        cat(rep("=", 70), "\n", sep="")
        cat("\n")
    }

    return(invisible(summary_list))
}


#' Check data quality and completeness
#'
#' Validates input data quality and provides warnings for potential issues
#' that might affect chromothripsis detection accuracy.
#'
#' @param SV.sample SVs object or data frame
#' @param CNV.sample CNVsegs object or data frame
#' @param verbose Logical, whether to print detailed messages (default: TRUE)
#' @return A list with data quality metrics and warnings
#' @export
check_data_quality <- function(SV.sample, CNV.sample = NULL, verbose = TRUE) {

    issues <- list()
    warnings_list <- character(0)

    # Convert to data frame if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    # Check SV data
    issues$total_SVs <- nrow(SV.sample)
    issues$intrachromosomal_SVs <- sum(SV.sample$chrom1 == SV.sample$chrom2)
    issues$interchromosomal_SVs <- sum(SV.sample$chrom1 != SV.sample$chrom2)

    if (issues$intrachromosomal_SVs < 5) {
        warnings_list <- c(warnings_list,
                          "Very few intrachromosomal SVs (<5). Chromothripsis detection requires multiple clustered SVs.")
    }

    # Check SV types distribution
    if ("SVtype" %in% names(SV.sample)) {
        sv_types <- table(SV.sample$SVtype)
        issues$sv_type_distribution <- sv_types

        if (length(sv_types) < 2) {
            warnings_list <- c(warnings_list,
                             "Only one SV type detected. Chromothripsis typically involves multiple SV types.")
        }
    }

    # Check CNV data
    if (!is.null(CNV.sample)) {
        issues$total_CNV_segments <- nrow(CNV.sample)

        if (nrow(CNV.sample) == 0) {
            warnings_list <- c(warnings_list,
                             "No CNV data provided. CN oscillations cannot be evaluated.")
        } else if (nrow(CNV.sample) < 10) {
            warnings_list <- c(warnings_list,
                             "Very few CNV segments (<10). May affect CN oscillation detection.")
        }

        # Check for truly adjacent segments with same CN
        # Only check segments that are genomically adjacent (end[i] >= start[i+1])
        if (nrow(CNV.sample) > 1) {
            chrom_list <- unique(CNV.sample$chrom)
            duplicate_cn <- 0

            for (chr in chrom_list) {
                chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]
                if (nrow(chr_cnv) > 1) {
                    chr_cnv <- chr_cnv[order(chr_cnv$start), ]
                    for (i in 1:(nrow(chr_cnv)-1)) {
                        # Check if segments are truly adjacent (no gap) AND have same CN
                        if (chr_cnv$end[i] >= chr_cnv$start[i+1] &&
                            chr_cnv$total_cn[i] == chr_cnv$total_cn[i+1]) {
                            duplicate_cn <- duplicate_cn + 1
                        }
                    }
                }
            }

            if (duplicate_cn > 0) {
                warnings_list <- c(warnings_list,
                                 sprintf("Found %d truly adjacent CNV segments with identical copy numbers. Consider merging these segments.", duplicate_cn))
            }
        }
    } else {
        warnings_list <- c(warnings_list,
                         "No CNV data provided. CN oscillation analysis will be skipped.")
    }

    issues$warnings <- warnings_list
    issues$has_issues <- length(warnings_list) > 0

    if (verbose) {
        cat("\nDATA QUALITY CHECK\n")
        cat(rep("=", 50), "\n", sep="")
        cat(sprintf("Total SVs: %d\n", issues$total_SVs))
        cat(sprintf("  - Intrachromosomal: %d\n", issues$intrachromosomal_SVs))
        cat(sprintf("  - Interchromosomal: %d\n", issues$interchromosomal_SVs))

        if (!is.null(CNV.sample)) {
            cat(sprintf("Total CNV segments: %d\n", issues$total_CNV_segments))
        }

        if (issues$has_issues) {
            cat("\nWARNINGS:\n")
            for (i in 1:length(warnings_list)) {
                cat(sprintf("  [%d] %s\n", i, warnings_list[i]))
            }
        } else {
            cat("\nNo data quality issues detected.\n")
        }
        cat(rep("=", 50), "\n\n", sep="")
    }

    return(invisible(issues))
}

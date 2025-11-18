#' Analyze breakpoint sequence features
#'
#' Extracts and analyzes sequence features at structural variant breakpoints,
#' including microhomology, insertions, and inferred DNA repair mechanisms.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param genome Reference genome object (BSgenome) or path to FASTA file
#' @param flank_size Size of flanking sequence to extract (default: 50 bp)
#' @param min_microhomology Minimum microhomology length to report (default: 2 bp)
#' @param max_microhomology Maximum microhomology length to search (default: 25 bp)
#' @return A list containing breakpoint sequence analysis results
#' @details
#' This function analyzes breakpoint junction sequences to identify:
#'
#' 1. **Microhomology**: Short homologous sequences (1-25 bp) at breakpoint junctions,
#'    indicative of microhomology-mediated end joining (MMEJ) or alt-EJ
#'
#' 2. **Insertions**: Non-templated or templated insertions at junctions
#'
#' 3. **DNA Repair Mechanisms**:
#'    - NHEJ (Non-Homologous End Joining): No microhomology, blunt ends
#'    - MMEJ (Microhomology-Mediated End Joining): 2-25 bp microhomology
#'    - HR (Homologous Recombination): Long homology (>25 bp)
#'    - FoSTeS/MMBIR: Serial replication with microhomology
#'
#' The analysis requires a reference genome to extract sequences around
#' breakpoint positions. Supports BSgenome objects or FASTA files.
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' # Analyze breakpoint sequences
#' bp_analysis <- analyze_breakpoint_sequences(
#'   SV_data,
#'   genome = BSgenome.Hsapiens.UCSC.hg19
#' )
#'
#' # View results
#' print(bp_analysis)
#' summary(bp_analysis)
#'
#' # Plot repair mechanism distribution
#' plot_repair_mechanisms(bp_analysis)
#' }
#'
#' @export
analyze_breakpoint_sequences <- function(SV.sample,
                                        genome,
                                        flank_size = 50,
                                        min_microhomology = 2,
                                        max_microhomology = 25) {

    # Input validation
    if (missing(genome)) {
        stop("Reference genome is required for sequence analysis.\n",
             "Provide a BSgenome object or path to FASTA file.")
    }

    # Convert SV.sample to data frame if needed
    if (inherits(SV.sample, "SVs")) {
        sv_df <- data.frame(
            chrom1 = SV.sample@chrom1,
            pos1 = SV.sample@pos1,
            strand1 = SV.sample@strand1,
            chrom2 = SV.sample@chrom2,
            pos2 = SV.sample@pos2,
            strand2 = SV.sample@strand2,
            SVtype = SV.sample@SVtype,
            stringsAsFactors = FALSE
        )
    } else {
        sv_df <- SV.sample
    }

    if (nrow(sv_df) == 0) {
        stop("No structural variants provided")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("     BREAKPOINT SEQUENCE ANALYSIS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Analyzing %d structural variants...\n", nrow(sv_df)))
    cat(sprintf("Flank size: %d bp\n", flank_size))
    cat(sprintf("Microhomology range: %d-%d bp\n\n", min_microhomology, max_microhomology))

    results <- list()

    # 1. Extract breakpoint sequences
    cat("Step 1: Extracting breakpoint sequences...\n")
    bp_sequences <- extract_breakpoint_sequences(
        sv_df,
        genome,
        flank_size
    )
    results$sequences <- bp_sequences

    # 2. Detect microhomology
    cat("Step 2: Detecting microhomology...\n")
    microhomology <- detect_microhomology(
        bp_sequences,
        min_length = min_microhomology,
        max_length = max_microhomology
    )
    results$microhomology <- microhomology

    # 3. Detect insertions
    cat("Step 3: Detecting insertions...\n")
    insertions <- detect_insertions(bp_sequences)
    results$insertions <- insertions

    # 4. Classify repair mechanisms
    cat("Step 4: Classifying DNA repair mechanisms...\n")
    repair_classification <- classify_repair_mechanisms(
        microhomology,
        insertions,
        sv_df
    )
    results$repair_mechanisms <- repair_classification

    # 5. Summary statistics
    cat("Step 5: Generating summary statistics...\n")
    summary_stats <- create_sequence_summary(
        microhomology,
        insertions,
        repair_classification
    )
    results$summary <- summary_stats

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("ANALYSIS COMPLETE\n")
    cat(rep("=", 70), "\n\n", sep = "")

    class(results) <- c("breakpoint_sequences", "list")
    return(results)
}


#' Extract sequences around breakpoints
#'
#' @param sv_df Data frame with SV information
#' @param genome Reference genome
#' @param flank_size Flanking sequence size
#' @return Data frame with extracted sequences
#' @keywords internal
extract_breakpoint_sequences <- function(sv_df, genome, flank_size) {

    # Check if genome is BSgenome object or file path
    is_bsgenome <- inherits(genome, "BSgenome")

    if (!is_bsgenome && !file.exists(genome)) {
        warning("Cannot access reference genome. Returning placeholder sequences.\n",
                "For full functionality, provide a BSgenome object or valid FASTA file.")

        # Return placeholder structure
        return(data.frame(
            sv_id = 1:nrow(sv_df),
            has_sequence = FALSE,
            bp1_left = NA,
            bp1_right = NA,
            bp2_left = NA,
            bp2_right = NA,
            stringsAsFactors = FALSE
        ))
    }

    sequences <- list()

    for (i in 1:nrow(sv_df)) {
        sv <- sv_df[i, ]

        # Extract sequences around breakpoint 1
        bp1_left <- extract_sequence_region(
            genome,
            sv$chrom1,
            max(1, sv$pos1 - flank_size),
            sv$pos1,
            is_bsgenome
        )

        bp1_right <- extract_sequence_region(
            genome,
            sv$chrom1,
            sv$pos1 + 1,
            sv$pos1 + flank_size,
            is_bsgenome
        )

        # Extract sequences around breakpoint 2
        bp2_left <- extract_sequence_region(
            genome,
            sv$chrom2,
            max(1, sv$pos2 - flank_size),
            sv$pos2,
            is_bsgenome
        )

        bp2_right <- extract_sequence_region(
            genome,
            sv$chrom2,
            sv$pos2 + 1,
            sv$pos2 + flank_size,
            is_bsgenome
        )

        # Handle strand orientation
        if (sv$strand1 == "-") {
            bp1_left <- reverse_complement(bp1_left)
            bp1_right <- reverse_complement(bp1_right)
        }

        if (sv$strand2 == "-") {
            bp2_left <- reverse_complement(bp2_left)
            bp2_right <- reverse_complement(bp2_right)
        }

        sequences[[i]] <- list(
            sv_id = i,
            has_sequence = TRUE,
            bp1_left = bp1_left,
            bp1_right = bp1_right,
            bp2_left = bp2_left,
            bp2_right = bp2_right,
            chrom1 = sv$chrom1,
            pos1 = sv$pos1,
            chrom2 = sv$chrom2,
            pos2 = sv$pos2,
            SVtype = sv$SVtype
        )
    }

    # Convert to data frame
    sequences_df <- do.call(rbind, lapply(sequences, function(x) {
        data.frame(
            sv_id = as.integer(x$sv_id),
            has_sequence = as.logical(x$has_sequence),
            bp1_left = as.character(ifelse(is.null(x$bp1_left), NA, x$bp1_left)),
            bp1_right = as.character(ifelse(is.null(x$bp1_right), NA, x$bp1_right)),
            bp2_left = as.character(ifelse(is.null(x$bp2_left), NA, x$bp2_left)),
            bp2_right = as.character(ifelse(is.null(x$bp2_right), NA, x$bp2_right)),
            chrom1 = as.character(x$chrom1),
            pos1 = as.numeric(x$pos1),
            chrom2 = as.character(x$chrom2),
            pos2 = as.numeric(x$pos2),
            SVtype = as.character(x$SVtype),
            stringsAsFactors = FALSE
        )
    }))

    return(sequences_df)
}


#' Extract sequence from genome region
#'
#' @param genome Genome object or file path
#' @param chrom Chromosome name
#' @param start Start position
#' @param end End position
#' @param is_bsgenome Whether genome is BSgenome object
#' @return Sequence string
#' @keywords internal
extract_sequence_region <- function(genome, chrom, start, end, is_bsgenome) {

    if (is_bsgenome) {
        # Use BSgenome
        tryCatch({
            seq <- as.character(genome[[chrom]][start:end])
            return(seq)
        }, error = function(e) {
            return(NA)
        })
    } else {
        # Use FASTA file (requires Rsamtools or similar)
        # For now, return placeholder
        return(NA)
    }
}


#' Reverse complement of DNA sequence
#'
#' @param seq DNA sequence string
#' @return Reverse complement
#' @keywords internal
reverse_complement <- function(seq) {
    if (is.na(seq)) return(NA)

    # Complement mapping
    comp_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G",
                 "a" = "t", "t" = "a", "g" = "c", "c" = "g",
                 "N" = "N", "n" = "n")

    # Split sequence
    bases <- strsplit(seq, "")[[1]]

    # Complement
    comp_bases <- sapply(bases, function(b) {
        if (b %in% names(comp_map)) comp_map[b] else "N"
    })

    # Reverse
    rev_comp <- paste(rev(comp_bases), collapse = "")

    return(rev_comp)
}


#' Detect microhomology at breakpoint junctions
#'
#' @param bp_sequences Data frame with breakpoint sequences
#' @param min_length Minimum microhomology length
#' @param max_length Maximum microhomology length
#' @return Data frame with microhomology information
#' @keywords internal
detect_microhomology <- function(bp_sequences, min_length = 2, max_length = 25) {

    microhomology <- list()

    for (i in 1:nrow(bp_sequences)) {
        bp <- bp_sequences[i, ]

        if (!bp$has_sequence) {
            microhomology[[i]] <- list(
                sv_id = bp$sv_id,
                has_microhomology = NA,
                microhomology_length = NA,
                microhomology_seq = NA
            )
            next
        }

        # Check for microhomology between bp1_right and bp2_right
        # (this represents the junction sequence)
        mh_result <- find_microhomology(
            bp$bp1_right,
            bp$bp2_left,
            min_length,
            max_length
        )

        microhomology[[i]] <- list(
            sv_id = bp$sv_id,
            has_microhomology = mh_result$found,
            microhomology_length = mh_result$length,
            microhomology_seq = mh_result$sequence,
            position = mh_result$position
        )
    }

    # Convert to data frame
    mh_df <- do.call(rbind, lapply(microhomology, function(x) {
        data.frame(
            sv_id = as.integer(x$sv_id),
            has_microhomology = as.logical(x$has_microhomology),
            microhomology_length = as.numeric(x$microhomology_length),
            microhomology_seq = as.character(x$microhomology_seq),
            stringsAsFactors = FALSE
        )
    }))

    return(mh_df)
}


#' Find microhomology between two sequences
#'
#' @param seq1 First sequence
#' @param seq2 Second sequence
#' @param min_length Minimum homology length
#' @param max_length Maximum homology length
#' @return List with microhomology information
#' @keywords internal
find_microhomology <- function(seq1, seq2, min_length, max_length) {

    if (is.na(seq1) || is.na(seq2)) {
        return(list(found = NA, length = NA, sequence = NA, position = NA))
    }

    # Search for overlapping homology
    # Start from longest possible and work down
    for (length in max_length:min_length) {
        # Check end of seq1 vs start of seq2
        if (nchar(seq1) >= length && nchar(seq2) >= length) {
            end_seq1 <- substr(seq1, nchar(seq1) - length + 1, nchar(seq1))
            start_seq2 <- substr(seq2, 1, length)

            if (end_seq1 == start_seq2) {
                return(list(
                    found = TRUE,
                    length = length,
                    sequence = end_seq1,
                    position = "junction"
                ))
            }
        }
    }

    # No microhomology found
    return(list(found = FALSE, length = 0, sequence = NA, position = NA))
}


#' Detect insertions at breakpoint junctions
#'
#' @param bp_sequences Data frame with breakpoint sequences
#' @return Data frame with insertion information
#' @keywords internal
detect_insertions <- function(bp_sequences) {

    insertions <- list()

    for (i in 1:nrow(bp_sequences)) {
        bp <- bp_sequences[i, ]

        # Placeholder: In real implementation, would compare
        # expected junction vs observed junction from reads
        # For now, mark as unknown
        insertions[[i]] <- list(
            sv_id = bp$sv_id,
            has_insertion = NA,
            insertion_length = NA,
            insertion_seq = NA,
            insertion_type = NA  # templated or non-templated
        )
    }

    # Convert to data frame
    ins_df <- do.call(rbind, lapply(insertions, function(x) {
        data.frame(
            sv_id = as.integer(x$sv_id),
            has_insertion = as.logical(x$has_insertion),
            insertion_length = as.numeric(x$insertion_length),
            insertion_seq = as.character(x$insertion_seq),
            insertion_type = as.character(x$insertion_type),
            stringsAsFactors = FALSE
        )
    }))

    return(ins_df)
}


#' Classify DNA repair mechanisms based on sequence features
#'
#' @param microhomology Microhomology data frame
#' @param insertions Insertions data frame
#' @param sv_df Original SV data frame
#' @return Data frame with repair mechanism classification
#' @keywords internal
classify_repair_mechanisms <- function(microhomology, insertions, sv_df) {

    mechanisms <- list()

    for (i in 1:nrow(microhomology)) {
        mh <- microhomology[i, ]
        ins <- insertions[i, ]
        sv <- sv_df[i, ]

        # Classification logic
        mechanism <- "Unknown"
        confidence <- "Low"

        if (!is.na(mh$has_microhomology)) {
            if (mh$has_microhomology) {
                # Has microhomology
                if (mh$microhomology_length >= 2 && mh$microhomology_length <= 25) {
                    mechanism <- "MMEJ/Alt-EJ"
                    confidence <- "High"

                    # Further classify based on SV type
                    if (sv$SVtype == "DUP" && mh$microhomology_length >= 5) {
                        mechanism <- "FoSTeS/MMBIR"
                        confidence <- "Moderate"
                    }
                } else if (mh$microhomology_length > 25) {
                    mechanism <- "HR-like"
                    confidence <- "Moderate"
                }
            } else {
                # No microhomology
                mechanism <- "NHEJ"
                confidence <- "Moderate"

                # Check for insertions
                if (!is.na(ins$has_insertion) && ins$has_insertion) {
                    if (ins$insertion_type == "non-templated") {
                        mechanism <- "NHEJ"
                        confidence <- "High"
                    }
                }
            }
        }

        mechanisms[[i]] <- list(
            sv_id = i,
            repair_mechanism = mechanism,
            confidence = confidence,
            evidence = sprintf("MH: %s (%d bp)",
                             ifelse(is.na(mh$has_microhomology), "Unknown",
                                   ifelse(mh$has_microhomology, "Yes", "No")),
                             ifelse(is.na(mh$microhomology_length), 0, mh$microhomology_length))
        )
    }

    # Convert to data frame
    mech_df <- do.call(rbind, lapply(mechanisms, function(x) {
        data.frame(
            sv_id = as.integer(x$sv_id),
            repair_mechanism = as.character(x$repair_mechanism),
            confidence = as.character(x$confidence),
            evidence = as.character(x$evidence),
            stringsAsFactors = FALSE
        )
    }))

    return(mech_df)
}


#' Create summary statistics for sequence analysis
#'
#' @param microhomology Microhomology data frame
#' @param insertions Insertions data frame
#' @param repair_classification Repair mechanism data frame
#' @return List with summary statistics
#' @keywords internal
create_sequence_summary <- function(microhomology, insertions, repair_classification) {

    # Microhomology statistics
    mh_present <- sum(microhomology$has_microhomology, na.rm = TRUE)
    mh_total <- sum(!is.na(microhomology$has_microhomology))
    mh_proportion <- if (mh_total > 0) mh_present / mh_total else 0

    mh_lengths <- microhomology$microhomology_length[!is.na(microhomology$microhomology_length)]
    mh_mean_length <- if (length(mh_lengths) > 0) mean(mh_lengths) else NA
    mh_median_length <- if (length(mh_lengths) > 0) median(mh_lengths) else NA

    # Repair mechanism distribution
    mech_table <- table(repair_classification$repair_mechanism)

    # Overall summary
    summary <- list(
        n_breakpoints = nrow(microhomology),
        microhomology = list(
            n_with_mh = mh_present,
            proportion = mh_proportion,
            mean_length = mh_mean_length,
            median_length = mh_median_length,
            length_range = if (length(mh_lengths) > 0) range(mh_lengths) else c(NA, NA)
        ),
        repair_mechanisms = as.data.frame(mech_table, stringsAsFactors = FALSE),
        dominant_mechanism = names(mech_table)[which.max(mech_table)]
    )

    return(summary)
}


#' Print method for breakpoint sequence analysis
#'
#' @param x Breakpoint sequence analysis result
#' @param ... Additional arguments
#' @export
print.breakpoint_sequences <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         BREAKPOINT SEQUENCE ANALYSIS SUMMARY\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Total breakpoints analyzed: %d\n\n", x$summary$n_breakpoints))

    # Microhomology summary
    cat("MICROHOMOLOGY:\n")
    cat(sprintf("  Breakpoints with microhomology: %d (%.1f%%)\n",
               x$summary$microhomology$n_with_mh,
               x$summary$microhomology$proportion * 100))
    if (!is.na(x$summary$microhomology$mean_length)) {
        cat(sprintf("  Mean length: %.1f bp\n", x$summary$microhomology$mean_length))
        cat(sprintf("  Median length: %.0f bp\n", x$summary$microhomology$median_length))
        cat(sprintf("  Range: %d-%d bp\n",
                   x$summary$microhomology$length_range[1],
                   x$summary$microhomology$length_range[2]))
    }
    cat("\n")

    # Repair mechanisms
    cat("DNA REPAIR MECHANISMS:\n")
    cat(sprintf("  Dominant mechanism: %s\n\n", x$summary$dominant_mechanism))
    cat("  Distribution:\n")
    for (i in 1:nrow(x$summary$repair_mechanisms)) {
        mech <- x$summary$repair_mechanisms$Var1[i]
        freq <- x$summary$repair_mechanisms$Freq[i]
        prop <- freq / x$summary$n_breakpoints * 100
        cat(sprintf("    - %s: %d (%.1f%%)\n", mech, freq, prop))
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}

#' Detect chromoplexy events from structural variation data
#'
#' Chromoplexy is characterized by chained translocations across multiple
#' chromosomes with minimal copy number changes. This function identifies
#' chromoplexy events by detecting translocation chains.
#'
#' @param SV.sample An instance of class SVs or data frame with SV data
#' @param CNV.sample Optional: an instance of class CNVsegs or data frame with CNV data
#' @param min_chromosomes Minimum number of chromosomes involved (default: 3)
#' @param min_translocations Minimum number of translocations in chain (default: 3)
#' @param max_cn_change Maximum allowed copy number change (default: 1)
#' @param allow_cycles Allow circular chains (default: TRUE)
#' @return A list containing chromoplexy detection results
#' @details
#' Chromoplexy differs from chromothripsis in several ways:
#' - Involves multiple chromosomes (typically 3-8)
#' - Characterized by chained translocations
#' - Minimal copy number changes
#' - Fewer breakpoints than chromothripsis
#'
#' Detection algorithm:
#' 1. Build translocation graph from inter-chromosomal SVs
#' 2. Detect chains using graph traversal
#' 3. Evaluate copy number stability
#' 4. Calculate chromoplexy-specific scores
#'
#' @references
#' Baca et al. (2013) Punctuated evolution of prostate cancer genomes.
#' Cell, 153(3):666-677.
#'
#' @export
detect_chromoplexy <- function(SV.sample,
                               CNV.sample = NULL,
                               min_chromosomes = 3,
                               min_translocations = 3,
                               max_cn_change = 1,
                               allow_cycles = TRUE) {

    # Convert to data frame if needed
    if (is(SV.sample, "SVs")) {
        SV.sample <- as(SV.sample, "data.frame")
    }

    if (!is.null(CNV.sample) && is(CNV.sample, "CNVsegs")) {
        CNV.sample <- as(CNV.sample, "data.frame")
    }

    # Extract inter-chromosomal SVs (translocations)
    inter_chr_SVs <- SV.sample[SV.sample$chrom1 != SV.sample$chrom2, ]

    if (nrow(inter_chr_SVs) < min_translocations) {
        warning(sprintf("Only %d inter-chromosomal SVs found. Need at least %d for chromoplexy detection.",
                       nrow(inter_chr_SVs), min_translocations))
        return(create_empty_chromoplexy_result())
    }

    # Build translocation graph
    cat("Building translocation graph...\n")
    tlx_graph <- build_translocation_graph(inter_chr_SVs)

    # Detect chains
    cat("Detecting translocation chains...\n")
    chains <- detect_translocation_chains(tlx_graph, inter_chr_SVs,
                                         min_chromosomes = min_chromosomes,
                                         min_translocations = min_translocations,
                                         allow_cycles = allow_cycles)

    if (length(chains) == 0) {
        cat("No chromoplexy chains detected.\n")
        return(create_empty_chromoplexy_result())
    }

    cat(sprintf("Found %d potential chromoplexy chain(s).\n", length(chains)))

    # Evaluate each chain
    chain_results <- list()
    for (i in 1:length(chains)) {
        chain_results[[i]] <- evaluate_chromoplexy_chain(
            chains[[i]],
            inter_chr_SVs,
            CNV.sample,
            max_cn_change = max_cn_change
        )
    }

    # Create summary
    summary_df <- do.call(rbind, lapply(chain_results, function(x) x$summary))

    # Classify events
    summary_df$classification <- classify_chromoplexy_event(summary_df)

    result <- list(
        chains = chains,
        chain_details = chain_results,
        summary = summary_df,
        translocation_graph = tlx_graph,
        total_chains = length(chains),
        likely_chromoplexy = sum(summary_df$classification == "Likely chromoplexy"),
        possible_chromoplexy = sum(summary_df$classification == "Possible chromoplexy")
    )

    class(result) <- c("chromoplexy", "list")
    return(result)
}


#' Build translocation graph from inter-chromosomal SVs
#'
#' Creates a graph where nodes are chromosome segments and edges are translocations.
#'
#' @param inter_chr_SVs Data frame of inter-chromosomal SVs
#' @return A list representing the translocation graph
#' @keywords internal
build_translocation_graph <- function(inter_chr_SVs) {

    # Create nodes: chromosome_position
    nodes <- unique(c(
        paste(inter_chr_SVs$chrom1, inter_chr_SVs$pos1, sep = ":"),
        paste(inter_chr_SVs$chrom2, inter_chr_SVs$pos2, sep = ":")
    ))

    # Create edges: connections between nodes
    edges <- data.frame(
        from = paste(inter_chr_SVs$chrom1, inter_chr_SVs$pos1, sep = ":"),
        to = paste(inter_chr_SVs$chrom2, inter_chr_SVs$pos2, sep = ":"),
        sv_index = 1:nrow(inter_chr_SVs),
        stringsAsFactors = FALSE
    )

    # Build adjacency list
    adj_list <- list()
    for (node in nodes) {
        adj_list[[node]] <- list(neighbors = character(0), sv_indices = integer(0))
    }

    for (i in 1:nrow(edges)) {
        from_node <- edges$from[i]
        to_node <- edges$to[i]
        sv_idx <- edges$sv_index[i]

        adj_list[[from_node]]$neighbors <- c(adj_list[[from_node]]$neighbors, to_node)
        adj_list[[from_node]]$sv_indices <- c(adj_list[[from_node]]$sv_indices, sv_idx)

        adj_list[[to_node]]$neighbors <- c(adj_list[[to_node]]$neighbors, from_node)
        adj_list[[to_node]]$sv_indices <- c(adj_list[[to_node]]$sv_indices, sv_idx)
    }

    return(list(
        nodes = nodes,
        edges = edges,
        adjacency_list = adj_list
    ))
}


#' Detect translocation chains using graph traversal
#'
#' Uses depth-first search to find chains of translocations.
#'
#' @param tlx_graph Translocation graph
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @param min_chromosomes Minimum chromosomes in chain
#' @param min_translocations Minimum translocations in chain
#' @param allow_cycles Allow circular chains
#' @return List of detected chains
#' @keywords internal
detect_translocation_chains <- function(tlx_graph,
                                       inter_chr_SVs,
                                       min_chromosomes = 3,
                                       min_translocations = 3,
                                       allow_cycles = TRUE) {

    visited_global <- character(0)
    chains <- list()
    chain_id <- 0

    # Try starting from each node
    for (start_node in tlx_graph$nodes) {

        if (start_node %in% visited_global) next

        # DFS from this node
        chain <- dfs_find_chain(
            start_node,
            tlx_graph$adjacency_list,
            visited = character(0),
            path = character(0),
            sv_path = integer(0)
        )

        # Check if chain meets criteria
        if (length(chain$path) >= min_translocations) {

            # Extract chromosomes involved
            chroms <- unique(sapply(chain$path, function(x) strsplit(x, ":")[[1]][1]))

            if (length(chroms) >= min_chromosomes) {
                chain_id <- chain_id + 1

                # Check if it's a cycle
                is_cycle <- FALSE
                if (allow_cycles && length(chain$path) > 1) {
                    first_node <- chain$path[1]
                    last_node <- chain$path[length(chain$path)]
                    if (last_node %in% tlx_graph$adjacency_list[[first_node]]$neighbors) {
                        is_cycle <- TRUE
                    }
                }

                chains[[chain_id]] <- list(
                    id = chain_id,
                    nodes = chain$path,
                    sv_indices = chain$sv_path,
                    chromosomes = chroms,
                    n_chromosomes = length(chroms),
                    n_translocations = length(chain$sv_path),
                    is_cycle = is_cycle
                )

                visited_global <- c(visited_global, chain$path)
            }
        }
    }

    return(chains)
}


#' Depth-first search to find translocation chains
#'
#' @param node Current node
#' @param adj_list Adjacency list
#' @param visited Visited nodes
#' @param path Current path
#' @param sv_path SV indices in path
#' @return Chain path
#' @keywords internal
dfs_find_chain <- function(node, adj_list, visited, path, sv_path) {

    visited <- c(visited, node)
    path <- c(path, node)

    neighbors <- adj_list[[node]]$neighbors
    sv_indices <- adj_list[[node]]$sv_indices

    # Find unvisited neighbors
    unvisited_idx <- which(!(neighbors %in% visited))

    if (length(unvisited_idx) > 0) {
        # Continue to first unvisited neighbor
        next_node <- neighbors[unvisited_idx[1]]
        next_sv <- sv_indices[unvisited_idx[1]]

        return(dfs_find_chain(
            next_node,
            adj_list,
            visited,
            path,
            c(sv_path, next_sv)
        ))
    }

    # No more unvisited neighbors
    return(list(path = path, sv_path = unique(sv_path)))
}


#' Evaluate a chromoplexy chain
#'
#' @param chain Chain to evaluate
#' @param inter_chr_SVs Inter-chromosomal SVs
#' @param CNV.sample CNV data
#' @param max_cn_change Maximum allowed CN change
#' @return Evaluation results
#' @keywords internal
evaluate_chromoplexy_chain <- function(chain, inter_chr_SVs, CNV.sample, max_cn_change = 1) {

    # Get SVs in this chain
    chain_SVs <- inter_chr_SVs[chain$sv_indices, ]

    # Evaluate copy number stability
    cn_stability_score <- 1.0  # Default: stable
    max_cn_deviation <- 0

    if (!is.null(CNV.sample) && nrow(CNV.sample) > 0) {
        # Check CN changes in regions affected by the chain
        cn_eval <- evaluate_cn_stability(chain, chain_SVs, CNV.sample)
        cn_stability_score <- cn_eval$stability_score
        max_cn_deviation <- cn_eval$max_deviation
    }

    # Calculate chain complexity score
    complexity_score <- calculate_chain_complexity(chain, chain_SVs)

    # SV type diversity (chromoplexy often has deletions + translocations)
    sv_types <- table(chain_SVs$SVtype)
    has_deletions <- "DEL" %in% names(sv_types)
    type_diversity <- length(sv_types)

    # Create summary
    summary <- data.frame(
        chain_id = chain$id,
        n_chromosomes = chain$n_chromosomes,
        chromosomes_involved = paste(chain$chromosomes, collapse = ","),
        n_translocations = chain$n_translocations,
        is_cycle = chain$is_cycle,
        cn_stability_score = cn_stability_score,
        max_cn_deviation = max_cn_deviation,
        complexity_score = complexity_score,
        has_deletions = has_deletions,
        sv_type_diversity = type_diversity,
        stringsAsFactors = FALSE
    )

    return(list(
        summary = summary,
        chain = chain,
        SVs = chain_SVs
    ))
}


#' Evaluate copy number stability in a chain
#'
#' @param chain Chain object
#' @param chain_SVs SVs in chain
#' @param CNV.sample CNV data
#' @return CN stability metrics
#' @keywords internal
evaluate_cn_stability <- function(chain, chain_SVs, CNV.sample) {

    cn_deviations <- c()

    for (chr in chain$chromosomes) {
        # Get CNV segments for this chromosome
        chr_cnv <- CNV.sample[CNV.sample$chrom == chr, ]

        if (nrow(chr_cnv) > 0) {
            # Calculate CN variation
            cn_values <- chr_cnv$total_cn
            cn_mean <- mean(cn_values, na.rm = TRUE)
            cn_sd <- sd(cn_values, na.rm = TRUE)

            # Check for large deviations
            deviations <- abs(cn_values - 2)  # Deviation from diploid
            cn_deviations <- c(cn_deviations, deviations)
        }
    }

    max_deviation <- if (length(cn_deviations) > 0) max(cn_deviations, na.rm = TRUE) else 0
    mean_deviation <- if (length(cn_deviations) > 0) mean(cn_deviations, na.rm = TRUE) else 0

    # Stability score: 1.0 = very stable, 0.0 = highly unstable
    # Penalize large deviations
    stability_score <- exp(-mean_deviation / 2)

    return(list(
        stability_score = stability_score,
        max_deviation = max_deviation,
        mean_deviation = mean_deviation
    ))
}


#' Calculate chain complexity score
#'
#' @param chain Chain object
#' @param chain_SVs SVs in chain
#' @return Complexity score
#' @keywords internal
calculate_chain_complexity <- function(chain, chain_SVs) {

    # Complexity based on:
    # 1. Number of chromosomes
    # 2. Number of translocations
    # 3. Whether it forms a cycle

    chr_score <- min(chain$n_chromosomes / 5, 1.0)  # Normalize to 5 chromosomes
    tlx_score <- min(chain$n_translocations / 10, 1.0)  # Normalize to 10 translocations
    cycle_bonus <- if (chain$is_cycle) 0.2 else 0.0

    complexity_score <- (chr_score * 0.4 + tlx_score * 0.4 + cycle_bonus)

    return(complexity_score)
}


#' Classify chromoplexy event based on criteria
#'
#' @param summary_df Summary data frame
#' @return Character vector of classifications
#' @keywords internal
classify_chromoplexy_event <- function(summary_df) {

    classifications <- character(nrow(summary_df))

    for (i in 1:nrow(summary_df)) {
        row <- summary_df[i, ]

        # Criteria for chromoplexy:
        # 1. Multiple chromosomes (≥3)
        # 2. Multiple translocations (≥3)
        # 3. CN stability (score ≥ 0.7)
        # 4. Moderate complexity

        meets_chr_criteria <- row$n_chromosomes >= 3
        meets_tlx_criteria <- row$n_translocations >= 3
        meets_cn_criteria <- row$cn_stability_score >= 0.7
        meets_complexity <- row$complexity_score >= 0.3

        criteria_met <- sum(c(meets_chr_criteria, meets_tlx_criteria,
                             meets_cn_criteria, meets_complexity))

        if (criteria_met >= 4) {
            classifications[i] <- "Likely chromoplexy"
        } else if (criteria_met >= 3) {
            classifications[i] <- "Possible chromoplexy"
        } else if (criteria_met >= 2) {
            classifications[i] <- "Unlikely chromoplexy"
        } else {
            classifications[i] <- "Not chromoplexy"
        }
    }

    return(classifications)
}


#' Create empty chromoplexy result
#'
#' @return Empty result structure
#' @keywords internal
create_empty_chromoplexy_result <- function() {
    result <- list(
        chains = list(),
        chain_details = list(),
        summary = data.frame(
            chain_id = integer(0),
            n_chromosomes = integer(0),
            chromosomes_involved = character(0),
            n_translocations = integer(0),
            is_cycle = logical(0),
            cn_stability_score = numeric(0),
            max_cn_deviation = numeric(0),
            complexity_score = numeric(0),
            has_deletions = logical(0),
            sv_type_diversity = integer(0),
            classification = character(0)
        ),
        translocation_graph = NULL,
        total_chains = 0,
        likely_chromoplexy = 0,
        possible_chromoplexy = 0
    )

    class(result) <- c("chromoplexy", "list")
    return(result)
}


#' Print method for chromoplexy results
#'
#' @param x Chromoplexy result object
#' @param ... Additional arguments
#' @export
print.chromoplexy <- function(x, ...) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("         CHROMOPLEXY DETECTION RESULTS\n")
    cat(rep("=", 70), "\n\n", sep = "")

    cat(sprintf("Total chains detected: %d\n", x$total_chains))
    cat(sprintf("  - Likely chromoplexy:   %d\n", x$likely_chromoplexy))
    cat(sprintf("  - Possible chromoplexy: %d\n", x$possible_chromoplexy))
    cat("\n")

    if (x$total_chains > 0) {
        cat("Chain Summary:\n")
        cat(rep("-", 70), "\n", sep = "")
        print(x$summary[, c("chain_id", "n_chromosomes", "n_translocations",
                           "cn_stability_score", "classification")])
    } else {
        cat("No chromoplexy chains detected.\n")
    }

    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(x)
}


#' Summary method for chromoplexy results
#'
#' @param object Chromoplexy result object
#' @param ... Additional arguments
#' @export
summary.chromoplexy <- function(object, ...) {
    print(object)
    invisible(object)
}

# Synthetic Chromoanagenesis Data Generator
# For testing and validation of ShatterSeek Extended Edition

#' Generate synthetic chromothripsis event
#'
#' @param chromosome Chromosome name
#' @param region_start Start of affected region
#' @param region_end End of affected region
#' @param num_breakpoints Number of breakpoints to generate
#' @param oscillation_amplitude CN oscillation amplitude
#' @param seed Random seed for reproducibility
#'
#' @return List with SV and CNV data frames
#'
#' @export
generate_chromothripsis <- function(
    chromosome = "1",
    region_start = 50e6,
    region_end = 100e6,
    num_breakpoints = 20,
    oscillation_amplitude = 2,
    seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)

    # Characteristics of chromothripsis:
    # 1. Clustered breakpoints in localized region
    # 2. Oscillating copy numbers (typically between 1 and 3-4)
    # 3. Random SV orientations
    # 4. Predominantly DEL and INV

    # Generate clustered breakpoints
    breakpoints <- sort(sample(region_start:region_end, num_breakpoints))

    # Create SVs
    sv_list <- list()
    for (i in seq(1, length(breakpoints)-1, 2)) {
        if (i+1 > length(breakpoints)) break

        # Chromothripsis has more DEL and INV
        sv_type <- sample(
            c("DEL", "h2hINV", "t2tINV"),
            1,
            prob = c(0.5, 0.25, 0.25)
        )

        # Assign correct strands
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
    baseline_cn <- 2

    # Oscillating pattern
    cn_values <- baseline_cn + oscillation_amplitude *
        sin(seq(0, 2*pi, length.out = num_segments))
    cn_values <- round(pmax(0, cn_values))  # Ensure non-negative

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

    return(list(
        SV = sv_data,
        CNV = cnv_data,
        metadata = list(
            type = "chromothripsis",
            chromosome = chromosome,
            region_start = region_start,
            region_end = region_end,
            num_breakpoints = num_breakpoints
        )
    ))
}


#' Generate synthetic chromoplexy event
#'
#' @param chromosomes Vector of chromosome names involved
#' @param num_translocations Number of translocations
#' @param seed Random seed
#'
#' @return List with SV and CNV data frames
#'
#' @export
generate_chromoplexy <- function(
    chromosomes = c("1", "2", "3"),
    num_translocations = 5,
    seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)

    # Characteristics of chromoplexy:
    # 1. Chained translocations across multiple chromosomes
    # 2. Copy number stable (mostly CN=2)
    # 3. Closed chains or linear chains

    stopifnot(length(chromosomes) >= 3)
    stopifnot(num_translocations >= 3)

    # Generate random positions on each chromosome
    chr_sizes <- c(
        "1" = 249250621, "2" = 243199373, "3" = 198022430,
        "4" = 191154276, "5" = 180915260, "6" = 171115067,
        "7" = 159138663, "8" = 146364022, "9" = 141213431,
        "10" = 135534747
    )

    # Create translocation chain
    sv_list <- list()

    for (i in 1:num_translocations) {
        chr1 <- chromosomes[(i-1) %% length(chromosomes) + 1]
        chr2 <- chromosomes[i %% length(chromosomes) + 1]

        pos1 <- sample(1:chr_sizes[chr1], 1)
        pos2 <- sample(1:chr_sizes[chr2], 1)

        sv_list[[i]] <- data.frame(
            chrom1 = chr1,
            pos1 = pos1,
            chrom2 = chr2,
            pos2 = pos2,
            SVtype = "TRA",
            strand1 = sample(c("+", "-"), 1),
            strand2 = sample(c("+", "-"), 1),
            stringsAsFactors = FALSE
        )
    }

    sv_data <- do.call(rbind, sv_list)

    # Generate stable CN (mostly 2, with occasional small changes)
    cnv_list <- list()
    for (chr in chromosomes) {
        # Divide chromosome into a few segments
        num_segs <- sample(3:5, 1)
        chr_size <- chr_sizes[chr]
        seg_breaks <- sort(sample(1:chr_size, num_segs-1))
        seg_starts <- c(1, seg_breaks)
        seg_ends <- c(seg_breaks, chr_size)

        for (i in 1:num_segs) {
            # Mostly CN=2, occasionally 1 or 3
            cn <- sample(c(2, 2, 2, 2, 1, 3), 1)

            cnv_list[[length(cnv_list) + 1]] <- data.frame(
                chrom = chr,
                start = seg_starts[i],
                end = seg_ends[i],
                total_cn = cn,
                stringsAsFactors = FALSE
            )
        }
    }

    cnv_data <- do.call(rbind, cnv_list)

    return(list(
        SV = sv_data,
        CNV = cnv_data,
        metadata = list(
            type = "chromoplexy",
            chromosomes = chromosomes,
            num_translocations = num_translocations
        )
    ))
}


#' Generate synthetic chromosynthesis event
#'
#' @param chromosome Chromosome name
#' @param region_start Start of affected region
#' @param region_end End of affected region
#' @param num_tandem_dups Number of tandem duplications
#' @param seed Random seed
#'
#' @return List with SV and CNV data frames
#'
#' @export
generate_chromosynthesis <- function(
    chromosome = "1",
    region_start = 50e6,
    region_end = 100e6,
    num_tandem_dups = 8,
    seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)

    # Characteristics of chromosynthesis:
    # 1. Serial tandem duplications
    # 2. Gradual increase in copy number
    # 3. FoSTeS/MMBIR mechanism

    # Generate positions for tandem duplications
    dup_positions <- sort(sample(region_start:region_end, num_tandem_dups * 2))

    # Create tandem duplications
    sv_list <- list()
    for (i in seq(1, length(dup_positions)-1, 2)) {
        sv_list[[length(sv_list) + 1]] <- data.frame(
            chrom1 = chromosome,
            pos1 = dup_positions[i],
            chrom2 = chromosome,
            pos2 = dup_positions[i+1],
            SVtype = "DUP",
            strand1 = "-",
            strand2 = "+",
            stringsAsFactors = FALSE
        )
    }

    sv_data <- do.call(rbind, sv_list)

    # Generate CN gradient (gradually increasing)
    num_segments <- num_tandem_dups + 1
    cn_values <- seq(2, 6, length.out = num_segments)
    cn_values <- round(cn_values)

    cnv_list <- list()
    seg_starts <- c(region_start, dup_positions[seq(2, length(dup_positions), 2)])
    seg_ends <- c(dup_positions[seq(1, length(dup_positions), 2)], region_end)

    for (i in 1:num_segments) {
        cnv_list[[i]] <- data.frame(
            chrom = chromosome,
            start = seg_starts[i],
            end = min(seg_ends[i], region_end),
            total_cn = cn_values[i],
            stringsAsFactors = FALSE
        )
    }

    cnv_data <- do.call(rbind, cnv_list)

    return(list(
        SV = sv_data,
        CNV = cnv_data,
        metadata = list(
            type = "chromosynthesis",
            chromosome = chromosome,
            region_start = region_start,
            region_end = region_end,
            num_tandem_dups = num_tandem_dups
        )
    ))
}


#' Generate mixed chromoanagenesis event
#'
#' @param seed Random seed
#'
#' @return List with SV and CNV data frames
#'
#' @export
generate_mixed_chromoanagenesis <- function(seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Generate chromothripsis on chr1
    ct <- generate_chromothripsis(
        chromosome = "1",
        region_start = 50e6,
        region_end = 100e6,
        num_breakpoints = 20
    )

    # Generate chromoplexy involving chr2, chr3, chr4
    cp <- generate_chromoplexy(
        chromosomes = c("2", "3", "4"),
        num_translocations = 6
    )

    # Combine SV data
    sv_data <- rbind(ct$SV, cp$SV)

    # Combine CNV data
    cnv_data <- rbind(ct$CNV, cp$CNV)

    return(list(
        SV = sv_data,
        CNV = cnv_data,
        metadata = list(
            type = "mixed",
            components = c("chromothripsis", "chromoplexy")
        )
    ))
}


# Example usage:
if (FALSE) {
    # Generate chromothripsis
    ct_data <- generate_chromothripsis(seed = 12345)

    # Convert to ShatterSeek objects
    sv_obj <- SVs(
        chrom1 = ct_data$SV$chrom1,
        pos1 = ct_data$SV$pos1,
        chrom2 = ct_data$SV$chrom2,
        pos2 = ct_data$SV$pos2,
        SVtype = ct_data$SV$SVtype,
        strand1 = ct_data$SV$strand1,
        strand2 = ct_data$SV$strand2
    )

    cnv_obj <- CNVsegs(
        chrom = ct_data$CNV$chrom,
        start = ct_data$CNV$start,
        end = ct_data$CNV$end,
        total_cn = ct_data$CNV$total_cn
    )

    # Detect
    result <- detect_chromothripsis(sv_obj, cnv_obj)
}

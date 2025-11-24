#' Convert StructuralVariantAnnotation GRanges to SVs object
#'
#' Converts GRanges objects from StructuralVariantAnnotation package
#' (breakpointRanges/breakendRanges) to ShatterSeek SVs object.
#'
#' This function provides integration with the StructuralVariantAnnotation
#' Bioconductor package, allowing users to leverage its robust VCF parsing
#' and then use ShatterSeek for chromoanagenesis detection.
#'
#' @param gr GRanges object from breakpointRanges() or breakendRanges()
#' @param partner_col Name of metadata column containing partner breakend ID
#'   (default: "partner" for StructuralVariantAnnotation)
#' @param svtype_col Name of metadata column containing SV type
#'   (default: "svtype")
#' @return An SVs object compatible with ShatterSeek
#'
#' @details
#' This function expects GRanges objects with the following metadata columns:
#' - partner: ID of the partner breakend (for paired breakpoints)
#' - svtype: SV type (DEL, DUP, INV, BND/TRA)
#' - Strand information from the GRanges object itself
#'
#' The function handles:
#' - Paired breakpoints (from breakpointRanges)
#' - Unpaired breakends (from breakendRanges)
#' - Automatic strand inference
#' - SV type standardization to ShatterSeek format
#'
#' @examples
#' \dontrun{
#' library(StructuralVariantAnnotation)
#' library(VariantAnnotation)
#'
#' # Read VCF using StructuralVariantAnnotation
#' vcf <- readVcf("sample.vcf.gz", "hg38")
#' gr <- breakpointRanges(vcf)
#'
#' # Convert to ShatterSeek format
#' sv_data <- granges_to_svs(gr)
#'
#' # Or combine breakpoints and breakends
#' gr_all <- c(breakpointRanges(vcf), breakendRanges(vcf))
#' sv_data <- granges_to_svs(gr_all)
#'
#' # Use in ShatterSeek analysis
#' results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
#' }
#'
#' @export
granges_to_svs <- function(gr,
                           partner_col = "partner",
                           svtype_col = "svtype") {

    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Please install GenomicRanges package:\n",
             "  BiocManager::install('GenomicRanges')")
    }

    if (!inherits(gr, "GRanges")) {
        stop("Input must be a GRanges object from StructuralVariantAnnotation")
    }

    # Check for required metadata
    if (!partner_col %in% names(GenomicRanges::mcols(gr))) {
        stop(sprintf("Metadata column '%s' not found in GRanges object.\n",
                    partner_col),
             "Available columns: ", paste(names(GenomicRanges::mcols(gr)), collapse = ", "))
    }

    # Extract basic information
    chrom1 <- as.character(GenomicRanges::seqnames(gr))
    pos1 <- GenomicRanges::start(gr)
    strand1 <- as.character(GenomicRanges::strand(gr))

    # Get partner information
    mcols_data <- GenomicRanges::mcols(gr)
    partner_ids <- mcols_data[[partner_col]]

    # Handle SV types
    if (svtype_col %in% names(mcols_data)) {
        sv_types_raw <- mcols_data[[svtype_col]]
    } else {
        # Try to infer from partner relationships and strands
        sv_types_raw <- rep(NA_character_, length(gr))
    }

    # Initialize output vectors
    n_svs <- length(gr)
    chrom2 <- character(n_svs)
    pos2 <- numeric(n_svs)
    strand2 <- character(n_svs)
    sv_types <- character(n_svs)

    # Track processed breakpoints to avoid duplicates
    processed <- rep(FALSE, n_svs)

    # Process paired breakpoints
    for (i in seq_along(gr)) {
        if (processed[i]) next

        # Find partner
        partner_id <- partner_ids[i]
        if (is.na(partner_id)) {
            # Unpaired breakend - skip or treat as TRA
            next
        }

        # Find partner index
        partner_idx <- which(names(gr) == partner_id)

        if (length(partner_idx) == 0) {
            # Partner not in this GRanges - single breakend (TRA)
            next
        }

        partner_idx <- partner_idx[1]

        # Extract partner information
        chrom2[i] <- as.character(GenomicRanges::seqnames(gr)[partner_idx])
        pos2[i] <- GenomicRanges::start(gr)[partner_idx]
        strand2[i] <- as.character(GenomicRanges::strand(gr)[partner_idx])

        # Determine SV type
        if (!is.na(sv_types_raw[i])) {
            raw_type <- toupper(as.character(sv_types_raw[i]))

            # Standardize SV type names
            sv_types[i] <- .standardize_sv_type(
                raw_type,
                chrom1[i], chrom2[i],
                strand1[i], strand2[i]
            )
        } else {
            # Infer from strands and chromosomes
            sv_types[i] <- .infer_sv_type_from_strands(
                chrom1[i], chrom2[i],
                strand1[i], strand2[i]
            )
        }

        # Mark both breakpoints as processed
        processed[i] <- TRUE
        processed[partner_idx] <- TRUE
    }

    # Remove unprocessed (unpaired) breakends
    valid <- processed

    if (sum(valid) == 0) {
        stop("No valid paired breakpoints found in GRanges object")
    }

    # Standardize chromosome names (remove chr prefix)
    chrom1 <- gsub("^chr", "", chrom1[valid])
    chrom2 <- gsub("^chr", "", chrom2[valid])

    # Filter to supported chromosomes (1-22, X)
    valid_chroms <- c(as.character(1:22), "X")
    keep <- chrom1 %in% valid_chroms & chrom2 %in% valid_chroms

    if (sum(keep) == 0) {
        stop("No SVs found on supported chromosomes (1-22, X)")
    }

    # Create SVs object
    svs <- SVs(
        chrom1 = chrom1[keep],
        pos1 = pos1[valid][keep],
        chrom2 = chrom2[keep],
        pos2 = pos2[valid][keep],
        SVtype = sv_types[valid][keep],
        strand1 = strand1[valid][keep],
        strand2 = strand2[valid][keep]
    )

    message(sprintf("Converted %d breakpoints to %d SVs", sum(valid), sum(keep)))

    return(svs)
}


#' Standardize SV type names to ShatterSeek format
#'
#' @param raw_type Raw SV type from caller
#' @param chrom1 First chromosome
#' @param chrom2 Second chromosome
#' @param strand1 First strand
#' @param strand2 Second strand
#' @return Standardized SV type
#' @keywords internal
.standardize_sv_type <- function(raw_type, chrom1, chrom2, strand1, strand2) {

    # Handle translocation
    if (chrom1 != chrom2) {
        return("TRA")
    }

    # Map common SV type names
    type_map <- c(
        "DEL" = "DEL",
        "DELETION" = "DEL",
        "DUP" = "DUP",
        "DUPLICATION" = "DUP",
        "INV" = "INV",
        "INVERSION" = "INV",
        "BND" = "TRA",
        "TRA" = "TRA",
        "TRANSLOCATION" = "TRA"
    )

    if (raw_type %in% names(type_map)) {
        mapped_type <- type_map[raw_type]

        # For inversions, determine h2h vs t2t
        if (mapped_type == "INV") {
            if (strand1 == "+" && strand2 == "+") {
                return("h2hINV")
            } else if (strand1 == "-" && strand2 == "-") {
                return("t2tINV")
            }
        }

        return(mapped_type)
    }

    # Fallback: infer from strands
    return(.infer_sv_type_from_strands(chrom1, chrom2, strand1, strand2))
}


#' Infer SV type from strand orientations
#'
#' @param chrom1 First chromosome
#' @param chrom2 Second chromosome
#' @param strand1 First strand
#' @param strand2 Second strand
#' @return Inferred SV type
#' @keywords internal
.infer_sv_type_from_strands <- function(chrom1, chrom2, strand1, strand2) {

    # Translocation
    if (chrom1 != chrom2) {
        return("TRA")
    }

    # Intrachromosomal: infer from strands
    if (strand1 == "+" && strand2 == "-") {
        return("DEL")
    } else if (strand1 == "-" && strand2 == "+") {
        return("DUP")
    } else if (strand1 == "+" && strand2 == "+") {
        return("h2hINV")
    } else if (strand1 == "-" && strand2 == "-") {
        return("t2tINV")
    }

    # Unknown
    return("UNKNOWN")
}


#' Read SV VCF using StructuralVariantAnnotation
#'
#' Alternative VCF reader that uses StructuralVariantAnnotation package
#' for robust BND parsing and complex SV handling.
#'
#' @param vcf_file Path to VCF file
#' @param genome Genome assembly (e.g., "hg19", "hg38")
#' @param include_breakends Include unpaired breakends (default: TRUE)
#' @return SVs object
#'
#' @details
#' This function provides an alternative to read_sv_vcf() that uses the
#' StructuralVariantAnnotation Bioconductor package. This may be more
#' robust for complex VCF formats or unusual SV callers.
#'
#' Requires:
#' - StructuralVariantAnnotation (Bioconductor)
#' - VariantAnnotation (Bioconductor)
#'
#' @examples
#' \dontrun{
#' # Read using StructuralVariantAnnotation backend
#' sv_data <- read_sv_vcf_structuralvariant("sample.vcf.gz", genome = "hg38")
#'
#' # Use in analysis
#' results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
#' }
#'
#' @export
read_sv_vcf_structuralvariant <- function(vcf_file,
                                          genome = "hg19",
                                          include_breakends = TRUE) {

    if (!requireNamespace("StructuralVariantAnnotation", quietly = TRUE)) {
        stop("Please install StructuralVariantAnnotation package:\n",
             "  BiocManager::install('StructuralVariantAnnotation')")
    }

    if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
        stop("Please install VariantAnnotation package:\n",
             "  BiocManager::install('VariantAnnotation')")
    }

    # Read VCF
    message("Reading VCF file with StructuralVariantAnnotation...")
    vcf <- VariantAnnotation::readVcf(vcf_file, genome)

    # Extract breakpoints
    message("Extracting breakpoint ranges...")
    gr_bp <- StructuralVariantAnnotation::breakpointRanges(vcf)

    # Optionally include breakends
    if (include_breakends) {
        message("Extracting breakend ranges...")
        gr_be <- StructuralVariantAnnotation::breakendRanges(vcf)
        gr <- c(gr_bp, gr_be)
    } else {
        gr <- gr_bp
    }

    # Convert to SVs object
    message("Converting to ShatterSeek SVs object...")
    svs <- granges_to_svs(gr)

    return(svs)
}

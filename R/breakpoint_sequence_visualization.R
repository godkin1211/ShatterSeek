#' Plot DNA repair mechanism distribution
#'
#' Creates a bar plot showing the distribution of inferred DNA repair
#' mechanisms across all breakpoints.
#'
#' @param breakpoint_result Result from analyze_breakpoint_sequences()
#' @param sample_name Sample name for plot title
#' @return A ggplot object
#' @export
plot_repair_mechanisms <- function(breakpoint_result,
                                  sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    if (!inherits(breakpoint_result, "breakpoint_sequences")) {
        stop("Input must be a breakpoint_sequences result object")
    }

    mech_data <- breakpoint_result$summary$repair_mechanisms
    colnames(mech_data) <- c("mechanism", "count")

    # Calculate proportions
    mech_data$proportion <- mech_data$count / sum(mech_data$count)

    # Color scheme for mechanisms
    mech_colors <- c(
        "NHEJ" = "#E41A1C",
        "MMEJ/Alt-EJ" = "#377EB8",
        "FoSTeS/MMBIR" = "#4DAF4A",
        "HR-like" = "#984EA3",
        "Unknown" = "#999999"
    )

    # Create plot
    p <- ggplot2::ggplot(mech_data, ggplot2::aes(x = mechanism, y = count, fill = mechanism))

    # Add bars
    p <- p + ggplot2::geom_col(width = 0.7, color = "black", size = 0.5)

    # Add count labels
    p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d\n(%.1f%%)", count, proportion * 100)),
        vjust = -0.5,
        size = 4,
        fontface = "bold"
    )

    # Colors
    p <- p + ggplot2::scale_fill_manual(
        values = mech_colors,
        name = "Repair Mechanism"
    )

    # Labels
    title_text <- sprintf("DNA Repair Mechanism Distribution%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Dominant mechanism: %s | Total breakpoints: %d",
        breakpoint_result$summary$dominant_mechanism,
        breakpoint_result$summary$n_breakpoints
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Repair Mechanism",
        y = "Number of Breakpoints"
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none"
        )

    # Set y-axis limits
    max_y <- max(mech_data$count) * 1.2
    p <- p + ggplot2::ylim(0, max_y)

    return(p)
}


#' Plot microhomology length distribution
#'
#' Creates a histogram showing the distribution of microhomology lengths
#' at breakpoint junctions.
#'
#' @param breakpoint_result Result from analyze_breakpoint_sequences()
#' @param sample_name Sample name for plot title
#' @param binwidth Histogram bin width (default: 1 bp)
#' @return A ggplot object
#' @export
plot_microhomology_distribution <- function(breakpoint_result,
                                           sample_name = "",
                                           binwidth = 1) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    mh_data <- breakpoint_result$microhomology

    # Filter to only those with microhomology
    mh_data <- mh_data[!is.na(mh_data$microhomology_length) &
                      mh_data$microhomology_length > 0, ]

    if (nrow(mh_data) == 0) {
        message("No microhomology detected to plot.")
        return(NULL)
    }

    # Create plot
    p <- ggplot2::ggplot(mh_data, ggplot2::aes(x = microhomology_length))

    # Add histogram
    p <- p + ggplot2::geom_histogram(
        binwidth = binwidth,
        fill = "#377EB8",
        color = "black",
        alpha = 0.7
    )

    # Add mean and median lines
    mean_mh <- mean(mh_data$microhomology_length, na.rm = TRUE)
    median_mh <- median(mh_data$microhomology_length, na.rm = TRUE)

    p <- p + ggplot2::geom_vline(
        xintercept = mean_mh,
        color = "red",
        linetype = "dashed",
        size = 1
    )

    p <- p + ggplot2::geom_vline(
        xintercept = median_mh,
        color = "darkgreen",
        linetype = "dashed",
        size = 1
    )

    # Add labels for mean and median
    max_y <- max(table(mh_data$microhomology_length))
    p <- p + ggplot2::annotate(
        "text",
        x = mean_mh,
        y = max_y * 0.9,
        label = sprintf("Mean: %.1f bp", mean_mh),
        color = "red",
        hjust = -0.1,
        size = 4,
        fontface = "bold"
    )

    p <- p + ggplot2::annotate(
        "text",
        x = median_mh,
        y = max_y * 0.8,
        label = sprintf("Median: %.0f bp", median_mh),
        color = "darkgreen",
        hjust = -0.1,
        size = 4,
        fontface = "bold"
    )

    # Labels
    title_text <- sprintf("Microhomology Length Distribution%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    subtitle_text <- sprintf(
        "Breakpoints with microhomology: %d/%d (%.1f%%)",
        nrow(mh_data),
        breakpoint_result$summary$n_breakpoints,
        breakpoint_result$summary$microhomology$proportion * 100
    )

    p <- p + ggplot2::labs(
        title = title_text,
        subtitle = subtitle_text,
        x = "Microhomology Length (bp)",
        y = "Number of Breakpoints"
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 11),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}


#' Plot repair mechanism by SV type
#'
#' Creates a stacked bar plot showing the relationship between
#' SV types and DNA repair mechanisms.
#'
#' @param breakpoint_result Result from analyze_breakpoint_sequences()
#' @param sv_data Original SV data to get SV types
#' @param sample_name Sample name for plot title
#' @return A ggplot object
#' @export
plot_repair_by_svtype <- function(breakpoint_result,
                                  sv_data,
                                  sample_name = "") {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plotting.")
    }

    # Combine repair mechanism with SV type
    repair_mech <- breakpoint_result$repair_mechanisms

    # Get SV types
    if (inherits(sv_data, "SVs")) {
        sv_types <- sv_data@SVtype
    } else {
        sv_types <- sv_data$SVtype
    }

    combined_data <- data.frame(
        SVtype = sv_types,
        repair_mechanism = repair_mech$repair_mechanism,
        stringsAsFactors = FALSE
    )

    # Create contingency table
    cont_table <- table(combined_data$SVtype, combined_data$repair_mechanism)
    plot_data <- as.data.frame(cont_table, stringsAsFactors = FALSE)
    colnames(plot_data) <- c("SVtype", "mechanism", "count")

    # Filter out zero counts
    plot_data <- plot_data[plot_data$count > 0, ]

    if (nrow(plot_data) == 0) {
        message("No data to plot.")
        return(NULL)
    }

    # Color scheme
    mech_colors <- c(
        "NHEJ" = "#E41A1C",
        "MMEJ/Alt-EJ" = "#377EB8",
        "FoSTeS/MMBIR" = "#4DAF4A",
        "HR-like" = "#984EA3",
        "Unknown" = "#999999"
    )

    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = SVtype, y = count, fill = mechanism))

    # Add stacked bars
    p <- p + ggplot2::geom_col(
        position = "stack",
        color = "black",
        size = 0.3
    )

    # Colors
    p <- p + ggplot2::scale_fill_manual(
        values = mech_colors,
        name = "Repair Mechanism"
    )

    # Labels
    title_text <- sprintf("Repair Mechanisms by SV Type%s",
                         ifelse(sample_name != "", paste0(" - ", sample_name), ""))

    p <- p + ggplot2::labs(
        title = title_text,
        x = "SV Type",
        y = "Number of Breakpoints"
    )

    # Theme
    p <- p + ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "right"
        )

    return(p)
}


#' Plot comprehensive breakpoint sequence report
#'
#' Combines multiple visualizations into a comprehensive report.
#'
#' @param breakpoint_result Result from analyze_breakpoint_sequences()
#' @param sv_data Original SV data
#' @param sample_name Sample name for plot title
#' @return A combined plot object
#' @export
plot_breakpoint_report <- function(breakpoint_result,
                                  sv_data,
                                  sample_name = "") {

    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for creating reports.")
    }

    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required for creating reports.")
    }

    # Create individual plots
    p1 <- plot_repair_mechanisms(breakpoint_result, sample_name)
    p2 <- plot_microhomology_distribution(breakpoint_result, sample_name)
    p3 <- plot_repair_by_svtype(breakpoint_result, sv_data, sample_name)

    # Combine plots
    plots_to_combine <- list()
    if (!is.null(p1)) plots_to_combine[[length(plots_to_combine) + 1]] <- p1
    if (!is.null(p2)) plots_to_combine[[length(plots_to_combine) + 1]] <- p2
    if (!is.null(p3)) plots_to_combine[[length(plots_to_combine) + 1]] <- p3

    if (length(plots_to_combine) == 0) {
        message("No plots could be generated.")
        return(NULL)
    }

    # Arrange plots
    if (length(plots_to_combine) == 3) {
        combined <- gridExtra::grid.arrange(
            plots_to_combine[[1]],
            plots_to_combine[[2]],
            plots_to_combine[[3]],
            ncol = 2,
            nrow = 2,
            layout_matrix = rbind(c(1, 2), c(3, 3)),
            top = grid::textGrob(
                sprintf("Breakpoint Sequence Analysis Report%s",
                       ifelse(sample_name != "", paste0(" - ", sample_name), "")),
                gp = grid::gpar(fontsize = 16, fontface = "bold")
            )
        )
    } else {
        combined <- gridExtra::grid.arrange(
            grobs = plots_to_combine,
            ncol = 2,
            top = grid::textGrob(
                sprintf("Breakpoint Sequence Analysis Report%s",
                       ifelse(sample_name != "", paste0(" - ", sample_name), "")),
                gp = grid::gpar(fontsize = 16, fontface = "bold")
            )
        )
    }

    return(combined)
}

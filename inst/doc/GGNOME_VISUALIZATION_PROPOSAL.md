# gGnome-Style Visualization Implementation for ShatterSeek Extended Edition

## Executive Summary

This document describes the implementation of gGnome-style integrated genome graph views in ShatterSeek Extended Edition, enhancing visualization capabilities while maintaining compatibility with the existing ggplot2-based framework.

## Background

### gGnome/gTrack Features
- **Graph-based representation** of structural variations
- **Track-based layout** for multi-dimensional data
- **Interactive exploration** capabilities
- **Genome-wide and region-specific** views

### ShatterSeek Extended Edition Implementation
ShatterSeek Extended Edition has implemented enhanced visualization capabilities that were originally proposed:
1. Integrated genome-wide overview
2. Graph-theoretic SV network visualization
3. Multi-track integrated region views

## Proposal: Three-Phase Implementation Plan

---

## Phase 1: Genome-Wide Dashboard (Priority: HIGH)

### Estimated Time: 1-2 weeks

### Goal
Create an integrated dashboard providing a complete overview of all chromoanagenesis events across the genome.

### Implementation

#### Function 1: `plot_genome_dashboard()`

```r
#' Comprehensive genome-wide chromoanagenesis dashboard
#'
#' @param chromoanag_result Result from detect_chromoanagenesis()
#' @param sample_name Sample identifier
#' @return Grid of plots showing genome-wide patterns
#' @export
plot_genome_dashboard <- function(chromoanag_result, sample_name = "") {

    # Panel 1: Ideogram with mechanism overlay
    p1 <- plot_genome_ideogram(chromoanag_result)

    # Panel 2: SV density heatmap by chromosome
    p2 <- plot_sv_density_heatmap(chromoanag_result)

    # Panel 3: CN profile across all chromosomes
    p3 <- plot_genome_wide_cn(chromoanag_result)

    # Panel 4: Mechanism distribution sunburst/treemap
    p4 <- plot_mechanism_distribution(chromoanag_result)

    # Combine into dashboard
    gridExtra::grid.arrange(p1, p2, p3, p4,
                           layout_matrix = rbind(c(1,1,2),
                                                c(3,3,2),
                                                c(4,4,4)))
}
```

#### Components:

**1.1 Genome Ideogram with Events**
```r
plot_genome_ideogram <- function(chromoanag_result) {
    # Karyogram-style visualization
    # Color-coded chromosomes by:
    # - No events (gray)
    # - Single mechanism (color by type)
    # - Mixed mechanisms (striped/gradient)

    # Implementation using ggplot2:
    # - geom_rect() for chromosome bands
    # - geom_point() or geom_segment() for events
    # - scale_color_manual() for mechanism types
}
```

**1.2 SV Density Heatmap**
```r
plot_sv_density_heatmap <- function(chromoanag_result) {
    # Chromosome × SV Type heatmap
    # Shows density of each SV type per chromosome
    # Uses geom_tile() with color gradient
}
```

**1.3 Genome-wide CN Profile**
```r
plot_genome_wide_cn <- function(chromoanag_result) {
    # Linear representation of all chromosomes
    # CN on y-axis, genomic position on x-axis
    # Highlight chromoanagenesis regions
}
```

**Benefits:**
- Quick overview of sample complexity
- Identifies chromosomes of interest
- Publication-quality figure
- Minimal new dependencies

---

## Phase 2: Integrated Track View (Priority: MEDIUM)

### Estimated Time: 2-3 weeks

### Goal
Create IGV/UCSC Genome Browser-style track visualization focused on chromoanagenesis features.

### Implementation

#### Function 2: `plot_region_tracks()`

```r
#' Multi-track visualization of genomic region
#'
#' @param chromoanag_result Chromoanagenesis results
#' @param chr Chromosome to display
#' @param start Start position (optional, auto if NULL)
#' @param end End position (optional, auto if NULL)
#' @param tracks Vector of tracks to display
#' @return Multi-panel plot with aligned tracks
#' @export
plot_region_tracks <- function(chromoanag_result,
                              chr,
                              start = NULL,
                              end = NULL,
                              tracks = c("ideogram", "genes", "cnv",
                                        "svs", "mechanisms", "stats")) {

    # Auto-detect region if not specified
    if (is.null(start) || is.null(end)) {
        region <- get_chromoanagenesis_region(chromoanag_result, chr)
        start <- region$start
        end <- region$end
    }

    plots <- list()

    # Track 1: Chromosome ideogram (context)
    if ("ideogram" %in% tracks) {
        plots$ideogram <- plot_ideogram_track(chr, start, end)
    }

    # Track 2: Gene annotations
    if ("genes" %in% tracks) {
        plots$genes <- plot_gene_track(chr, start, end)
    }

    # Track 3: Copy number variations
    if ("cnv" %in% tracks) {
        plots$cnv <- plot_cnv_track(chromoanag_result, chr, start, end)
    }

    # Track 4: Structural variations
    if ("svs" %in% tracks) {
        plots$svs <- plot_sv_track(chromoanag_result, chr, start, end)
    }

    # Track 5: Chromoanagenesis mechanisms
    if ("mechanisms" %in% tracks) {
        plots$mechanisms <- plot_mechanism_track(chromoanag_result, chr, start, end)
    }

    # Track 6: Statistical criteria
    if ("stats" %in% tracks) {
        plots$stats <- plot_stats_track(chromoanag_result, chr)
    }

    # Align all tracks by genomic position
    align_tracks(plots, xlim = c(start, end))
}
```

#### Track Components:

**2.1 Gene Track**
```r
plot_gene_track <- function(chr, start, end, genome = "hg19") {
    # Option 1: Use TxDb annotation packages
    # Option 2: Load from external gene annotation file
    # Display as rectangles with gene names
}
```

**2.2 CNV Track**
```r
plot_cnv_track <- function(chromoanag_result, chr, start, end) {
    # Step plot of copy number
    # Highlight oscillating regions
    # Color by CN state
}
```

**2.3 SV Track**
```r
plot_sv_track <- function(chromoanag_result, chr, start, end) {
    # Arc diagram for intrachromosomal SVs
    # Arrows/lines to chromosome edges for interchromosomal
    # Color and shape by SV type
}
```

**2.4 Mechanism Track**
```r
plot_mechanism_track <- function(chromoanag_result, chr, start, end) {
    # Annotation bars showing mechanism classifications
    # Color by mechanism type
    # Transparency by confidence
}
```

**Benefits:**
- Detailed region exploration
- Contextualize chromoanagenesis with genes
- Compare multiple samples side-by-side
- Identify breakpoint hotspots

---

## Phase 3: Graph Network Visualization (Priority: LOW)

### Estimated Time: 3-4 weeks

### Goal
Implement graph-theoretic visualization of SV networks similar to gGnome's genome graph representation.

### Implementation

#### Function 3: `plot_sv_network()`

```r
#' Network graph visualization of structural variations
#'
#' @param chromoanag_result Chromoanagenesis results
#' @param chr Chromosome(s) to include (NULL = all)
#' @param layout Layout algorithm ("fr", "kk", "circle", "tree")
#' @param show_cn Include copy number as node attributes
#' @return ggraph plot object
#' @export
plot_sv_network <- function(chromoanag_result,
                            chr = NULL,
                            layout = "fr",
                            show_cn = TRUE) {

    # Build graph structure
    sv_graph <- build_sv_graph(chromoanag_result, chr)

    # Create ggraph visualization
    library(ggraph)
    library(igraph)

    ggraph(sv_graph, layout = layout) +
        # Nodes = genomic segments
        geom_node_point(aes(color = mechanism, size = cn)) +
        # Edges = SVs
        geom_edge_link(aes(color = sv_type),
                      arrow = arrow(length = unit(2, 'mm'))) +
        # Labels
        geom_node_text(aes(label = region), repel = TRUE) +
        # Theme
        theme_graph() +
        labs(title = "Structural Variation Network")
}
```

#### Components:

**3.1 Graph Construction**
```r
build_sv_graph <- function(chromoanag_result, chr = NULL) {
    # Nodes: genomic segments between breakpoints
    nodes <- extract_genomic_segments(chromoanag_result, chr)

    # Edges: structural variations connecting segments
    edges <- extract_sv_connections(chromoanag_result, chr)

    # Build igraph object
    graph <- igraph::graph_from_data_frame(
        d = edges,
        vertices = nodes,
        directed = TRUE
    )

    # Add attributes
    V(graph)$cn <- get_copy_number(nodes)
    V(graph)$mechanism <- get_mechanism_annotation(nodes)
    E(graph)$sv_type <- get_sv_types(edges)

    return(graph)
}
```

**3.2 Network Layouts**
```r
# Fruchterman-Reingold (force-directed)
layout = "fr"

# Kamada-Kawai (stress minimization)
layout = "kk"

# Circular (show cycles)
layout = "circle"

# Tree (hierarchical)
layout = "tree"

# Custom genomic layout (preserve genomic order)
layout = "genomic"  # Custom implementation
```

**3.3 Interactive Version**
```r
plot_sv_network_interactive <- function(...) {
    # Use plotly or networkD3 for interactivity
    # Click nodes to see details
    # Hover for tooltip information
    # Zoom and pan
}
```

**Benefits:**
- Identify complex SV patterns
- Detect chromoplexy chains visually
- Understand genome graph topology
- Research tool for mechanism discovery

---

## Technical Implementation Details

### Dependencies

#### Required (already in ShatterSeek):
- ggplot2
- grid
- gridExtra
- GenomicRanges

#### New Dependencies:
- **Phase 1**: None
- **Phase 2**:
  - TxDb.Hsapiens.UCSC.hg19.knownGene (optional, for genes)
  - AnnotationDbi (optional)
- **Phase 3**:
  - igraph
  - ggraph
  - plotly or networkD3 (for interactivity)

### Code Organization

```
R/
├── chromothripsis_detection.R
├── chromothripsis_visualization.R
├── chromoplexy_detection.R
├── chromoplexy_visualization.R
├── chromosynthesis_detection.R
├── chromosynthesis_visualization.R
├── integrated_classifier.R
├── integrated_classifier_visualization.R
├── breakpoint_sequence_analysis.R
├── breakpoint_sequence_visualization.R
├── genome_dashboard_visualization.R        # NEW - Phase 1
├── region_tracks_visualization.R           # NEW - Phase 2
├── sv_network_visualization.R              # NEW - Phase 3
└── visualization_utils.R                   # NEW - Helper functions
```

### Helper Functions

```r
# visualization_utils.R

#' Get genomic coordinates of chromoanagenesis region
get_chromoanagenesis_region <- function(result, chr)

#' Align multiple ggplot tracks by x-axis
align_tracks <- function(plot_list, xlim)

#' Format genomic positions for axis labels
format_genomic_position <- function(pos)

#' Extract color palette for mechanisms
get_mechanism_colors <- function()

#' Convert genomic ranges to plot coordinates
gr_to_plot_coords <- function(gr, xlim)
```

---

## Comparison with gGnome

### What We Gain:

| Feature | gGnome | ShatterSeek Enhanced |
|---------|--------|---------------------|
| **Chromothripsis Detection** | ❌ | ✅ |
| **Chromoplexy Detection** | ❌ | ✅ |
| **Chromosynthesis Detection** | ❌ | ✅ |
| **Integrated Classification** | ❌ | ✅ |
| **Breakpoint Mechanism Analysis** | ❌ | ✅ |
| **Graph Visualization** | ✅ | ✅ (Phase 3) |
| **Track-based View** | ✅ | ✅ (Phase 2) |
| **Genome-wide Dashboard** | Partial | ✅ (Phase 1) |
| **Interactive Plots** | ✅ | Optional (Phase 3) |
| **Dependency Weight** | Heavy | Light-Medium |

### What We Keep:

- **Consistency**: All using ggplot2 grammar
- **Integration**: Works seamlessly with existing functions
- **Simplicity**: No learning curve for current users
- **Publication-ready**: High-quality static plots
- **Flexibility**: Easy to customize and extend

---

## Example Usage

### Phase 1 Usage:

```r
library(ShatterSeek)

# Run analysis
results <- detect_chromoanagenesis(SV_data, CN_data)

# Create dashboard
dashboard <- plot_genome_dashboard(results, sample_name = "Patient_123")
print(dashboard)

# Save
ggsave("chromoanagenesis_dashboard.pdf", dashboard, width = 16, height = 12)
```

### Phase 2 Usage:

```r
# Focus on chromosome with chromothripsis
tracks <- plot_region_tracks(
    chromoanag_result = results,
    chr = "chr21",
    tracks = c("genes", "cnv", "svs", "mechanisms")
)
print(tracks)
```

### Phase 3 Usage:

```r
# Visualize chromoplexy as network
chromoplexy_chains <- results$chromoplexy

network <- plot_sv_network(
    chromoanag_result = results,
    chr = c("chr2", "chr5", "chr8"),  # Chromosomes in chain
    layout = "fr",
    show_cn = TRUE
)
print(network)

# Interactive version
network_interactive <- plot_sv_network_interactive(results)
htmlwidgets::saveWidget(network_interactive, "sv_network.html")
```

---

## Prioritization and Timeline

### Immediate (Phase 1): 1-2 weeks
✅ **High Impact, Low Effort**
- Genome-wide dashboard
- Uses existing infrastructure
- Provides immediate value

### Near-term (Phase 2): 2-3 weeks
⚠️ **Medium Impact, Medium Effort**
- Track-based region view
- Requires gene annotation integration
- Enhances exploratory analysis

### Future (Phase 3): 3-4 weeks
🔵 **Lower Priority**
- Graph network visualization
- More specialized use case
- For advanced users and research

---

## Recommendations

### Recommended Approach: **Phase 1 Immediately + Phase 2 Conditionally**

1. **Implement Phase 1 now**
   - Quick wins with high user value
   - No new dependencies
   - Improves publication figures

2. **Gather user feedback**
   - Do users need region-level views?
   - Is gene annotation important?
   - What resolution is needed?

3. **Decide on Phase 2/3 based on demand**
   - If users want detailed exploration → Phase 2
   - If users want network analysis → Phase 3
   - Both can be optional extensions

### Alternative: Export to gGnome

If graph-based analysis is critical, consider:

```r
#' Export ShatterSeek results to gGnome format
#' @export
export_to_gGnome <- function(chromoanag_result) {
    # Convert to gGraph object
    # Let users leverage gGnome/gTrack directly
    # Maintains compatibility without heavy integration
}
```

This provides a bridge to gGnome ecosystem without reinventing the wheel.

---

## Conclusion

The proposed three-phase approach allows ShatterSeek to:

1. **Enhance** visualization capabilities incrementally
2. **Maintain** consistency with existing ggplot2 framework
3. **Avoid** heavy dependencies on experimental packages
4. **Provide** gGnome-style insights with ShatterSeek-specific focus

**Recommendation**: Start with Phase 1 (genome dashboard) as it provides maximum value with minimal effort, then evaluate Phase 2/3 based on user needs and feedback.

---

## Appendix: Code Skeleton for Phase 1

See separate file: `R/genome_dashboard_visualization.R` for implementation template.

---

**Author**: ShatterSeek Development Team
**Date**: 2025
**Status**: Proposal - Awaiting approval

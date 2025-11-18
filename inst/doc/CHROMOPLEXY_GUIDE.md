# Chromoplexy Detection in ShatterSeek

This guide describes the new chromoplexy detection functionality added to ShatterSeek, enabling comprehensive chromoanagenesis analysis.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Chromoplexy vs Chromothripsis](#chromoplexy-vs-chromothripsis)
3. [Detection Algorithm](#detection-algorithm)
4. [Usage Examples](#usage-examples)
5. [Interpreting Results](#interpreting-results)
6. [Visualization](#visualization)
7. [References](#references)

---

## Introduction

**Chromoplexy** is a complex chromosomal rearrangement pattern characterized by chained translocations across multiple chromosomes with minimal copy number changes. It represents a distinct mechanism of chromoanagenesis, different from chromothripsis.

ShatterSeek now includes dedicated functionality to detect and visualize chromoplexy events, complementing its existing chromothripsis detection capabilities.

### Key Features

- **Translocation chain detection**: Identifies chains of inter-chromosomal rearrangements
- **Multi-chromosome analysis**: Detects events spanning 3-8 chromosomes
- **Copy number integration**: Evaluates CN stability as a key criterion
- **Cycle detection**: Identifies circular chromoplexy patterns
- **Integrated analysis**: Simultaneous detection of both chromothripsis and chromoplexy

---

## Chromoplexy vs Chromothripsis

| Feature | Chromothripsis | Chromoplexy |
|---------|---------------|-------------|
| **Mechanism** | Chromosome shattering | Chained translocations |
| **Chromosomes** | 1-2 (localized) | 3-8 (multiple) |
| **Breakpoints** | Many (10+) | Moderate (3-10) |
| **Copy Number** | Oscillating (2-3 states) | Stable (minimal change) |
| **SV Types** | DEL, DUP, INV (random) | TRA + DEL (chained) |
| **Pattern** | Clustered, interleaved | Linear chains or cycles |
| **Timing** | Single catastrophic event | Successive events |

### Biological Context

**Chromothripsis**: Results from a single catastrophic event where one or two chromosomes shatter and are randomly reassembled, creating oscillating copy number patterns.

**Chromoplexy**: Involves a series of chained translocations across multiple chromosomes, often associated with deletions at breakpoints, but maintaining overall copy number stability.

---

## Detection Algorithm

### Overview

The chromoplexy detection algorithm consists of four main steps:

1. **Build Translocation Graph**
2. **Detect Chains**
3. **Evaluate CN Stability**
4. **Calculate Scores and Classify**

### Step 1: Build Translocation Graph

Creates a graph representation where:
- **Nodes**: Chromosome positions (breakpoints)
- **Edges**: Translocations connecting breakpoints

```
Example:
chr1:1000000 --[TRA]--> chr5:2000000
chr5:3000000 --[TRA]--> chr12:500000
chr12:1000000 --[TRA]--> chr1:1500000 (forms a cycle)
```

### Step 2: Detect Chains

Uses **Depth-First Search (DFS)** to find connected paths:
- Traverse the graph from each unvisited node
- Build chains by following translocation connections
- Identify both linear chains and circular chains (cycles)

**Criteria for a valid chain**:
- ≥ 3 translocations (default, customizable)
- ≥ 3 chromosomes involved (default, customizable)

### Step 3: Evaluate CN Stability

For each chain, assess copy number changes in affected regions:

**CN Stability Score** = exp(-mean_deviation / 2)

Where:
- mean_deviation = average deviation from diploid CN (2)
- Score ranges from 0 (unstable) to 1 (stable)
- Chromoplexy typically has scores > 0.7

### Step 4: Calculate Scores and Classify

**Complexity Score** (0-1):
- Chromosome participation: 40%
- Translocation count: 40%
- Cycle bonus: 20%

**Classification Criteria**:

| Classification | Criteria |
|----------------|----------|
| **Likely chromoplexy** | ≥4 of: Multiple chromosomes (≥3), Multiple translocations (≥3), CN stable (≥0.7), Moderate complexity (≥0.3) |
| **Possible chromoplexy** | 3 of above criteria |
| **Unlikely chromoplexy** | 2 of above criteria |
| **Not chromoplexy** | <2 criteria |

---

## Usage Examples

### Basic Chromoplexy Detection

```r
library(ShatterSeek)

# Prepare data
SV_data <- SVs(chrom1=..., pos1=..., chrom2=..., pos2=..., ...)
CN_data <- CNVsegs(chrom=..., start=..., end=..., total_cn=...)

# Detect chromoplexy
results <- detect_chromoplexy(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    min_chromosomes = 3,
    min_translocations = 3,
    max_cn_change = 1
)

# View results
print(results)
```

### Comprehensive Chromoanagenesis Analysis

```r
# Detect both chromothripsis and chromoplexy
all_results <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    genome = "hg19",
    detect_chromothripsis = TRUE,
    detect_chromoplexy = TRUE
)

# View integrated results
print(all_results)
summary(all_results)
```

### Custom Thresholds

```r
# Stringent detection
stringent <- detect_chromoplexy(
    SV_data, CN_data,
    min_chromosomes = 4,        # More chromosomes
    min_translocations = 5,     # More translocations
    max_cn_change = 0.5         # Stricter CN stability
)

# Lenient detection
lenient <- detect_chromoplexy(
    SV_data, CN_data,
    min_chromosomes = 2,
    min_translocations = 2,
    max_cn_change = 2
)
```

---

## Interpreting Results

### Result Structure

```r
chromoplexy_result <- detect_chromoplexy(...)

# Access components:
chromoplexy_result$total_chains              # Number of chains detected
chromoplexy_result$likely_chromoplexy        # Number of likely events
chromoplexy_result$possible_chromoplexy      # Number of possible events
chromoplexy_result$summary                   # Data frame with chain details
chromoplexy_result$chain_details             # Detailed information per chain
```

### Summary Table Columns

| Column | Description |
|--------|-------------|
| `chain_id` | Unique identifier for each chain |
| `n_chromosomes` | Number of chromosomes in chain |
| `chromosomes_involved` | Comma-separated chromosome list |
| `n_translocations` | Number of translocations in chain |
| `is_cycle` | Whether chain forms a closed loop |
| `cn_stability_score` | Copy number stability (0-1) |
| `max_cn_deviation` | Maximum deviation from diploid |
| `complexity_score` | Overall complexity metric (0-1) |
| `has_deletions` | Deletions present in chain |
| `sv_type_diversity` | Number of different SV types |
| `classification` | Likely/Possible/Unlikely chromoplexy |

### Key Indicators of Chromoplexy

✅ **Strong Evidence**:
- ≥3 chromosomes involved
- ≥3 chained translocations
- CN stability score ≥ 0.7
- Forms a cycle (circular pattern)
- Associated deletions at breakpoints

⚠️ **Moderate Evidence**:
- 2-3 chromosomes
- 2-3 translocations
- CN stability 0.5-0.7
- Linear chain

❌ **Weak Evidence**:
- Only 2 chromosomes
- Only 1-2 translocations
- CN stability < 0.5
- Isolated events

---

## Visualization

### Linear Chain Plot

```r
# Plot all chains
plots <- plot_chromoplexy(
    chromoplexy_result = results,
    sample_name = "Sample001",
    genome = "hg19"
)

# Plot specific chain
plot_chromoplexy(results, chain_id = 1, sample_name = "Sample001")
```

**Features**:
- Horizontal chromosome tracks
- Curved arcs representing translocations
- Breakpoint markers
- Chain statistics in subtitle

### Circular Plot (for cycles)

```r
# Plot circular chromoplexy pattern
circular_plot <- plot_chromoplexy_circular(
    chromoplexy_result = results,
    chain_id = 1,
    sample_name = "Sample001"
)
```

**Features**:
- Chromosomes arranged in circle
- Straight lines showing connections
- Clear visualization of cycle structure

---

## Integrated Analysis Workflow

### Complete Pipeline

```r
# 1. Load and prepare data
data(DO17373)
SV_data <- SVs(...)
CN_data <- CNVsegs(...)

# 2. Run comprehensive analysis
results <- detect_chromoanagenesis(SV_data, CN_data)

# 3. View results
print(results)
summary(results)

# 4. Compare mechanisms
if (results$chromothripsis$n_likely > 0 &&
    results$chromoplexy$likely_chromoplexy > 0) {
    cat("Both mechanisms detected - complex chromoanagenesis!\n")
}

# 5. Visualize
# Chromothripsis
for (chr in results$chromothripsis$classification$chrom) {
    plot_chromothripsis(results$chromothripsis$detection_output,
                       chr=chr)
}

# Chromoplexy
chromoplexy_plots <- plot_chromoplexy(results$chromoplexy)
```

### Export Results

```r
# Export chromoplexy summary
write.csv(results$chromoplexy$summary,
         "chromoplexy_results.csv",
         row.names = FALSE)

# Export integrated summary
write.csv(results$integrated_summary,
         "chromoanagenesis_summary.csv",
         row.names = FALSE)
```

---

## Performance and Limitations

### Performance

- **Speed**: Fast graph-based algorithm
- **Scalability**: Handles samples with hundreds of SVs
- **Memory**: Moderate usage (proportional to SV count)

### Limitations

1. **Requires inter-chromosomal SVs**: At least 3 translocations needed
2. **Short-read limitations**: Cannot resolve all complex rearrangements
3. **CN dependency**: Best results with high-quality CNV calls
4. **Chain ambiguity**: Multiple valid chain interpretations possible

### Best Practices

✅ **DO**:
- Use whole-genome sequencing data
- Include both SV and CNV data
- Verify results with visual inspection
- Use appropriate thresholds for your data type

❌ **DON'T**:
- Apply to exome-only data (insufficient coverage)
- Rely solely on automated classification
- Ignore data quality warnings
- Use extremely lenient thresholds

---

## Clinical and Research Applications

### Cancer Genomics

Chromoplexy is frequently observed in:
- **Prostate cancer** (highly prevalent)
- **Breast cancer** (ER+ subtype)
- **Lymphomas**
- **Multiple myeloma**

### Significance

- **Prognostic marker**: Associated with specific cancer subtypes
- **Therapeutic targets**: May indicate actionable rearrangements
- **Evolution studies**: Reveals clonal architecture
- **Mechanism insights**: Different from chromothripsis

---

## References

1. **Baca SC, et al. (2013)** Punctuated evolution of prostate cancer genomes. *Cell*, 153(3):666-677.
   - Original description of chromoplexy

2. **Cortes-Ciriano I, et al. (2020)** Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nature Genetics*, 52(3):331-341.
   - ShatterSeek publication, chromothripsis detection

3. **Notta F, et al. (2016)** A renewed model of pancreatic cancer evolution based on genomic rearrangement patterns. *Nature*, 538(7625):378-382.
   - Chromoanagenesis in pancreatic cancer

4. **Umbreit NT, et al. (2020)** Mechanisms generating cancer genome complexity from a single cell division error. *Science*, 368(6488):eaba0712.
   - Mechanisms of chromothripsis and chromoplexy

---

## Support and Contact

For questions, issues, or feature requests related to chromoplexy detection:

- **GitHub Issues**: [ShatterSeek repository]
- **Documentation**: This guide and `inst/examples/chromoplexy_usage.R`
- **Original Authors**: See ShatterSeek README

---

## Version History

- **v1.0**: Initial chromoplexy detection implementation
  - Graph-based chain detection
  - CN stability evaluation
  - Visualization functions
  - Integrated chromoanagenesis analysis

---

*Last updated: 2025*

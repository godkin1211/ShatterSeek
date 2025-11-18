# ShatterSeek

**ShatterSeek** is a comprehensive R package for detecting and analyzing **chromoanagenesis** events from next-generation sequencing (NGS) data. Chromoanagenesis represents catastrophic genomic rearrangements including:

- **Chromothripsis**: Chromosome shattering with random reassembly
- **Chromoplexy**: Chained translocations across multiple chromosomes
- **Chromosynthesis**: Serial replication-based rearrangements (FoSTeS/MMBIR)

ShatterSeek takes copy number (CN) and structural variation (SV) calls as input, making it compatible with virtually any CN and SV caller.

## Key Features

### Core Detection Modules
- **Chromothripsis Detection**: Identifies localized chromosome shattering events with oscillating copy numbers
- **Chromoplexy Detection**: Detects inter-chromosomal translocation chains with stable copy numbers
- **Chromosynthesis Detection**: Identifies replication-based complex rearrangements with CN gradients

### Advanced Analysis
- **Integrated Classification**: Analyzes mixed mechanism events where multiple chromoanagenesis types co-occur
- **Breakpoint Sequence Analysis**: Infers DNA repair mechanisms (NHEJ, MMEJ, HR, FoSTeS/MMBIR) from junction sequences
- **Confidence Scoring**: Evidence-based classification with quantitative confidence metrics
- **Complexity Assessment**: Multi-dimensional scoring of genomic instability

### Visualization
- **Comprehensive Plotting**: Publication-quality visualizations for all chromoanagenesis mechanisms
- **Interactive Reports**: Integrated dashboards combining multiple analysis views
- **Mechanism Landscape**: Genome-wide visualization of chromoanagenesis distribution

### Performance
- Fast analysis: ~20-30s per sample
- Scalable for large cohort studies
- Modular architecture for easy extension

## Scientific Background

ShatterSeek was originally developed and validated using ~2,600 whole-genome sequencing datasets from The Pan-Cancer Analysis of Whole Genomes (PCAWG) project.

**Original Publication:**
*Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing*
[Cortes-Ciriano et al. Nature Genetics, 2020](https://www.nature.com/articles/s41588-019-0576-7)

**Interactive Browser:** http://compbio.med.harvard.edu/chromothripsis/

**Extended Capabilities:**
The current version extends beyond chromothripsis to provide comprehensive chromoanagenesis analysis, including chromoplexy and chromosynthesis detection modules with integrated classification and mechanistic insights.

## Prerequisites

ShatterSeek requires R (>= 3.0.1) and depends on:
- **Core**: methods, BiocGenerics, S4Vectors, IRanges, GenomicRanges
- **Analysis**: graph, MASS, foreach
- **Visualization**: ggplot2, grid, gridExtra

**Optional** (for breakpoint sequence analysis):
- BSgenome packages (e.g., BSgenome.Hsapiens.UCSC.hg19)

## Installation

### From GitHub
```bash
$ git clone https://github.com/parklab/ShatterSeek.git
$ R CMD INSTALL ShatterSeek
```

### Using devtools
```R
library(devtools)
install_github("parklab/ShatterSeek")
```

## Quick Start

### Basic Chromothripsis Detection

```R
library(ShatterSeek)
data(DO17373)

# Prepare data
SV_data <- SVs(
    chrom1 = as.character(SV_DO17373$chrom1),
    pos1 = as.numeric(SV_DO17373$start1),
    chrom2 = as.character(SV_DO17373$chrom2),
    pos2 = as.numeric(SV_DO17373$end2),
    SVtype = as.character(SV_DO17373$svclass),
    strand1 = as.character(SV_DO17373$strand1),
    strand2 = as.character(SV_DO17373$strand2)
)

CN_data <- CNVsegs(
    chrom = as.character(SCNA_DO17373$chromosome),
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

# Run chromothripsis detection
chromothripsis <- shatterseek(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg19"
)

# View results
chromothripsis@chromSummary
```

### Comprehensive Chromoanagenesis Analysis

```R
# Detect all three chromoanagenesis mechanisms
results <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    genome = "hg19",
    detect_chromothripsis = TRUE,
    detect_chromoplexy = TRUE,
    detect_chromosynthesis = TRUE
)

# View comprehensive results
print(results)
summary(results)
```

### Mixed Mechanism Classification

```R
# Classify complex mixed mechanism events
mixed_class <- classify_mixed_mechanisms(
    chromoanagenesis_result = results,
    overlap_threshold = 1e6,
    min_confidence = 0.3
)

# View classification
print(mixed_class)

# Generate integrated report
plot_mechanism_report(mixed_class, sample_name = "Sample_01")
```

### Breakpoint Sequence Analysis

```R
# Requires reference genome
library(BSgenome.Hsapiens.UCSC.hg19)

# Analyze breakpoint sequences and infer repair mechanisms
bp_analysis <- analyze_breakpoint_sequences(
    SV.sample = SV_data,
    genome = BSgenome.Hsapiens.UCSC.hg19,
    flank_size = 50,
    min_microhomology = 2,
    max_microhomology = 25
)

# View results
print(bp_analysis)

# Visualize repair mechanisms
plot_repair_mechanisms(bp_analysis, sample_name = "Sample_01")
plot_microhomology_distribution(bp_analysis, sample_name = "Sample_01")

# Generate comprehensive report
plot_breakpoint_report(bp_analysis, SV_data, sample_name = "Sample_01")
```

## Detailed Usage

### 1. Chromothripsis Detection

```R
# Standard detection
chromothripsis <- shatterseek(SV_data, CN_data, genome = "hg19")

# Calculate confidence scores
scores <- calculate_confidence_score(chromothripsis@chromSummary)

# Classify chromosomes
classification <- classify_chromothripsis(
    chromothripsis,
    min_cluster_size = 6,
    min_oscillations = 4,
    max_pval_clustering = 0.05
)

# Visualize specific chromosome
plots <- plot_chromothripsis(
    ShatterSeek_output = chromothripsis,
    chr = "chr21",
    sample_name = "Sample_01",
    genome = "hg19"
)
```

### 2. Chromoplexy Detection

```R
# Detect chromoplexy
chromoplexy_result <- detect_chromoplexy(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    min_chromosomes = 3,
    min_translocations = 3,
    max_cn_change = 1
)

# View results
print(chromoplexy_result)

# Visualize translocation chains
plot_chromoplexy(chromoplexy_result, chain_id = 1)

# Circular plot for cycles
plot_chromoplexy_circular(chromoplexy_result, chain_id = 1)
```

### 3. Chromosynthesis Detection

```R
# Detect chromosynthesis
chromosyn_result <- detect_chromosynthesis(
    CNV.sample = CN_data,
    SV.sample = SV_data,
    min_segments = 5,
    gradient_threshold = 0.5,
    min_tandem_dups = 3
)

# View results
print(chromosyn_result)

# Visualize CN gradients
plot_chromosynthesis(chromosyn_result, region_id = 1)
plot_chromosynthesis_heatmap(chromosyn_result, region_id = 1)
plot_cn_gradient_scatter(chromosyn_result, region_id = 1)
```

### 4. Integrated Analysis Workflow

```R
# Complete analysis pipeline
# Step 1: Detect all mechanisms
results <- detect_chromoanagenesis(SV_data, CN_data)

# Step 2: Classify mixed mechanisms
mixed_class <- classify_mixed_mechanisms(results)

# Step 3: Analyze breakpoint sequences
bp_analysis <- analyze_breakpoint_sequences(
    SV_data,
    genome = BSgenome.Hsapiens.UCSC.hg19
)

# Step 4: Integrated interpretation
cat("\n=== Integrated Results ===\n")
cat("Classification:", mixed_class$sample_classification$classification, "\n")
cat("Complexity:", mixed_class$complexity$complexity_level, "\n")
cat("Dominant repair:", bp_analysis$summary$dominant_mechanism, "\n")

# Step 5: Generate comprehensive reports
plot_mechanism_report(mixed_class, "Sample_01")
plot_breakpoint_report(bp_analysis, SV_data, "Sample_01")
```

## Important Notes

### Data Preparation

1. **Minimum cluster size**: The `min.Size` parameter in `shatterseek()` defaults to 1. This ensures detection of chromothripsis events that might involve linked clusters across multiple chromosomes.

2. **CNV segment merging**: ShatterSeek expects adjacent CNV segments to have different copy numbers. Merge adjacent segments with identical copy numbers using:

```R
## d is a data.frame with columns: chr, start, end, total_cn
library(GenomicRanges)
dd <- d
dd$total_cn[dd$total_cn == 0] <- 150000
dd$total_cn[is.na(dd$total_cn)] <- 0
dd <- as(dd, "GRanges")
cov <- coverage(dd, weight = dd$total_cn)
dd1 <- as(cov, "GRanges")
dd1 <- as.data.frame(dd1)
dd1 <- dd1[dd1$score != 0, ]
dd1 <- dd1[, c(1,2,3,6)]
names(dd1) <- names(d)[1:4]
dd1$total_cn[dd1$total_cn == 150000] <- 0
d <- dd1
```

3. **SV type validation**: ShatterSeek validates SV type-strand consistency:
   - DEL: +/- strand combination
   - DUP: -/+ strand combination
   - h2hINV: +/+ strand combination
   - t2tINV: -/- strand combination

4. **Chromosome naming**: Supports both "chr1" and "1" formats. Automatically detects and handles both conventions.

## Clinical Applications

### Treatment Selection
- **MMEJ-dominant samples**: Consider PARP inhibitors
- **FoSTeS/MMBIR patterns**: ATR/CHK1 inhibitors for replication stress
- **HR-deficient tumors**: Platinum-based chemotherapy

### Risk Stratification
- **Complexity score**: Quantifies overall genomic chaos
- **Mixed mechanisms**: Indicates severe instability and poor prognosis
- **Repair mechanism profile**: Predicts therapy response

### Research Applications
- Cohort-level chromoanagenesis characterization
- Mechanism evolution studies
- Biomarker discovery
- Therapeutic target identification

## Documentation

- **Tutorial**: See `./inst/tutorial/tutorial.pdf` for detailed guide
- **Improvements**: See `./inst/doc/IMPROVEMENTS.md` for new features
- **Advanced Features**: See `./inst/doc/ADVANCED_FEATURES_GUIDE.md` for detailed explanations
- **Examples**: See `./inst/examples/` for usage scripts

## Output Files

ShatterSeek generates multiple output types:

1. **Data frames**: Quantitative results for each chromosome/region
2. **Classification objects**: S3/S4 objects with structured results
3. **Plots**: ggplot2 objects for publication-quality figures
4. **Reports**: Combined multi-panel visualizations

## Performance Tips

- For large cohorts, process samples in parallel using `foreach` or `parallel`
- Cache intermediate results for iterative analysis
- Use `min_confidence` thresholds to filter low-quality events
- Breakpoint sequence analysis is optional (requires reference genome)

## Troubleshooting

### Common Issues

**Problem**: "Arguments imply differing number of rows"
**Solution**: Ensure CNV segments are properly formatted and merged

**Problem**: Visualization fails with chromosome naming errors
**Solution**: Use consistent chromosome naming ("chr1" or "1" throughout)

**Problem**: No chromoplexy/chromosynthesis detected
**Solution**: Adjust detection thresholds or check input data quality

### Debug Functions

```R
# Debug chromoanagenesis results
debug_chromoanagenesis_result(results)

# Debug mixed mechanisms classification
debug_mixed_mechanisms(mixed_class)

# Check data quality
quality <- check_data_quality(SV_data, CN_data)
print(quality)
```

## Citation

If you use ShatterSeek in your research, please cite:

**Cortes-Ciriano, I., Lee, J.J., Xi, R. et al.** Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat Genet* **52**, 331–341 (2020). https://doi.org/10.1038/s41588-019-0576-7

For chromoplexy and chromosynthesis features, please also cite relevant method papers:
- **Chromoplexy**: Baca, S.C. et al. *Cell* (2013)
- **Chromosynthesis**: Liu, P. et al. *Cell* (2011)

## License

ShatterSeek is **free for academic use only**.

For non-academic or commercial use, please contact:
**Dr. Sonalee Barthakur**
Harvard University Office of Technology Development
Email: hms_otd@harvard.edu

## Contact

**For scientific questions:**
- Isidro Cortes Ciriano: isidrolauscher@gmail.com or icortes@ebi.ac.uk
- Peter J Park: peter_park@hms.harvard.edu

**For technical issues:**
- GitHub Issues: https://github.com/parklab/ShatterSeek/issues

## Acknowledgments

Development of the extended chromoanagenesis analysis modules was supported by contributions from the genomics community and feedback from active users.

---

**Version**: 2.0+ (with extended chromoanagenesis analysis)
**Last Updated**: 2025

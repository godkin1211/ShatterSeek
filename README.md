# ShatterSeek - Extended Edition

**ShatterSeek Extended** is a comprehensive R package for detecting and analyzing **chromoanagenesis** events from next-generation sequencing (NGS) data. Built upon the foundational chromothripsis detection framework developed by Cortes-Ciriano et al. (Nature Genetics, 2020), this extended version provides a complete chromoanagenesis analysis suite including:

- **Chromothripsis**: Chromosome shattering with random reassembly (original implementation + enhancements)
- **Chromoplexy**: Chained translocations across multiple chromosomes (NEW)
- **Chromosynthesis**: Serial replication-based rearrangements (FoSTeS/MMBIR) (NEW)

ShatterSeek Extended takes copy number (CN) and structural variation (SV) calls as input, making it compatible with virtually any CN and SV caller. The extended version adds **advanced visualization**, **direct VCF support with BND parsing**, **breakpoint sequence analysis**, and **integrated multi-mechanism classification**.

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
- **Genome-Wide Dashboard**: Multi-panel view with chromosome ideogram, SV density, CN profile, and mechanism distribution
- **gGnome-style Regional Plots**: Detailed regional views showing SV arcs and copy number profiles
- **Circos Plots**: Circular genome-wide visualization for multi-chromosome chromoanagenesis events
- **Mechanism Landscape**: Genome-wide visualization of chromoanagenesis distribution

### Input Flexibility
- **VCF Support**: Direct input from popular SV callers (Manta, Delly, GRIDSS, LUMPY, Sniffles)
- **CNV VCF Support**: Automatic parsing from CNV callers (CNVkit, GATK, Control-FREEC, Canvas)
- **Auto-detection**: Intelligent caller format recognition
- **Multi-sample VCF**: Extract specific samples from cohort VCFs

### Performance
- Fast analysis: ~20-30s per sample
- Scalable for large cohort studies
- Modular architecture for easy extension

## Scientific Background

### Original Framework
ShatterSeek was originally developed and validated using ~2,600 whole-genome sequencing datasets from The Pan-Cancer Analysis of Whole Genomes (PCAWG) project for chromothripsis detection.

**Original Publication:**
*Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing*
[Cortes-Ciriano et al. Nature Genetics, 2020](https://www.nature.com/articles/s41588-019-0576-7)

**Interactive Browser:** http://compbio.med.harvard.edu/chromothripsis/

### Extended Edition (Version 2.0+)

This **Extended Edition** represents a major enhancement of the original ShatterSeek framework with the following significant additions:

**New Detection Modules:**
- **Chromoplexy Detection**: Complete implementation of translocation chain detection across multiple chromosomes with confidence scoring
- **Chromosynthesis Detection**: Replication-based rearrangement analysis with CN gradient detection
- **Mixed Mechanism Classification**: Integrated analysis of co-occurring chromoanagenesis events

**Advanced Visualization System:**
- **Genome-Wide Dashboards**: Multi-panel views with anatomically accurate chromosome ideograms (p-arm, centromere, q-arm)
- **gGnome-style Regional Plots**: Detailed regional visualization with SV arcs and copy number profiles
- **Circos Plots**: Circular genome-wide visualization optimized for multi-chromosome events with four-track layout (chromosome ideogram, copy number, chromoanagenesis regions, SV links)

**Enhanced Input/Output:**
- **Direct VCF Support**: Automatic parsing from major SV callers (Manta, Delly, DRAGEN, GRIDSS, LUMPY, Sniffles)
- **BND/Translocation Parsing**: Robust parsing of breakend records from VCF ALT fields using VCF 4.2+ bracket notation
- **CNV VCF Support**: Direct input from CNVkit, GATK, Control-FREEC, Canvas

**Mechanistic Analysis:**
- **Breakpoint Sequence Analysis**: DNA repair mechanism inference (NHEJ, MMEJ, HR, FoSTeS/MMBIR)
- **Confidence Scoring**: Evidence-based quantitative metrics for all detection modules
- **Complexity Assessment**: Multi-dimensional genomic instability scoring

The extended edition maintains full backward compatibility with the original ShatterSeek API while adding comprehensive multi-mechanism chromoanagenesis analysis capabilities.

## Prerequisites

ShatterSeek requires R (>= 3.0.1) and depends on:
- **Core**: methods, BiocGenerics, S4Vectors, IRanges, GenomicRanges
- **Analysis**: graph, MASS, foreach
- **Visualization**: ggplot2, grid, gridExtra

**Optional Dependencies**:
- **VCF support**: VariantAnnotation or vcfR (for direct VCF input)
- **Circos plots**: circlize (for circular genome-wide visualization)
- **Breakpoint sequence analysis**: BSgenome packages (e.g., BSgenome.Hsapiens.UCSC.hg19)

## Installation

### From GitHub
```bash
$ git clone https://github.com/godkin1211/ShatterSeek.git
$ R CMD INSTALL ShatterSeek
```

### Using devtools
```R
library(devtools)
install_github("godkin1211/ShatterSeek")
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
chromothripsis <- detect_chromothripsis(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg19"
)

# View results
chromothripsis@chromSummary

# Classify chromothripsis events (High/Low confidence)
classification <- classify_chromothripsis(chromothripsis)
print(classification)

# Generate detailed summary report
summarize_chromothripsis(chromothripsis)
```

**Classification Levels** (based on Cortes-Ciriano et al. 2020):

- **High confidence**: ≥6 intrachromosomal SVs, ≥7 oscillating CN segments, fragment joins test passed (p>0.05), clustering test passed
- **Low confidence**: ≥6 intrachromosomal SVs, 4-6 oscillating CN segments, fragment joins test passed, clustering test passed
- **Not chromothripsis**: Does not meet above criteria

See `CHROMOTHRIPSIS_CLASSIFICATION.md` for detailed classification criteria.

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

## VCF Input Support

ShatterSeek now supports **direct VCF input** from popular SV and CNV callers, eliminating the need for manual data conversion.

### Supported SV Callers
- **Manta**: Illumina short-read SV caller
- **Delly**: Pan-SV caller with precise breakpoint detection
- **DRAGEN**: Illumina pipeline with comprehensive BND/TRA support
- **GRIDSS**: High-sensitivity breakpoint detection
- **LUMPY**: Probabilistic SV detection
- **Sniffles**: Long-read (PacBio/ONT) SV caller
- **Generic**: Any VCF with standard SVTYPE/END fields

**BND/Translocation Support**: ShatterSeek automatically parses breakend (BND) records from VCF ALT fields using standard VCF 4.2+ bracket notation (e.g., `[chr13:49291490[`, `]chr5:15886744]`), enabling accurate detection of inter-chromosomal events crucial for chromoplexy analysis.

### Supported CNV Callers
- **CNVkit**: Hybrid capture and WGS CNV detection
- **GATK**: Genome Analysis Toolkit CNV caller
- **Control-FREEC**: Copy number and allelic content detection
- **Canvas**: Illumina CNV caller
- **Generic**: Any VCF with CN field

### Basic VCF Workflow

```R
# Read VCF files directly (caller auto-detection)
sv_data <- read_sv_vcf("sample.sv.vcf.gz")
cnv_data <- read_cnv_vcf("sample.cnv.vcf.gz")

# Run chromoanagenesis detection
results <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    genome = "hg38"
)
```

### Advanced VCF Options

```R
# Specify caller explicitly and set filters
sv_data <- read_sv_vcf(
    vcf_file = "sample.manta.vcf.gz",
    caller = "manta",              # Explicitly specify caller
    min_sv_size = 1000,            # Filter SVs < 1kb
    include_tra = TRUE             # Include translocations
)

# CNV with custom field
cnv_data <- read_cnv_vcf(
    vcf_file = "sample.cnv.vcf.gz",
    caller = "cnvkit",
    cn_field = "CN",               # Specify CN field name
    merge_adjacent = TRUE          # Merge adjacent segments with same CN
)

# Multi-sample VCF
sv_data <- read_sv_vcf(
    vcf_file = "cohort.vcf.gz",
    sample_name = "SAMPLE_001"     # Extract specific sample
)
```

### Complete VCF-based Analysis Example

```R
library(ShatterSeek)

# 1. Read data from VCF files
sv_data <- read_sv_vcf(
    vcf_file = "patient_001.manta.vcf.gz",
    caller = "manta",
    min_sv_size = 1000
)

cnv_data <- read_cnv_vcf(
    vcf_file = "patient_001.cnvkit.vcf.gz",
    caller = "cnvkit"
)

# 2. Quality check
quality <- check_data_quality(sv_data, cnv_data)

# 3. Detect all chromoanagenesis types
results <- detect_chromoanagenesis(
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    genome = "hg38",
    verbose = TRUE
)

# 4. Classify mixed mechanisms
mixed <- classify_mixed_mechanisms(results)

# 5. Visualize
plot_mechanism_landscape(mixed, "Patient_001")
plot_complexity_breakdown(mixed, "Patient_001")
```

### VCF Installation Requirements

VCF support requires one of the following packages:

**Option 1: VariantAnnotation (recommended)**
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
```

**Option 2: vcfR (lightweight alternative)**
```R
install.packages("vcfR")
```

For more examples, see `inst/examples/vcf_workflow.R`.

## Detailed Usage

### 1. Chromothripsis Detection

```R
# Standard detection
chromothripsis <- detect_chromothripsis(SV_data, CN_data, genome = "hg19")

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

### 5. Advanced Visualization

ShatterSeek provides multiple visualization approaches for comprehensive chromoanagenesis analysis.

#### Genome-Wide Dashboard

Multi-panel overview showing chromosome ideogram with true chromosome shapes, SV density heatmap, copy number profile, and mechanism distribution.

```R
# Generate genome-wide dashboard
plot_genome_dashboard(
    chromoanagenesis_result = results,
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    sample_name = "Sample_01"
)
```

The dashboard includes:
- **Chromosome Ideogram**: Displays chromoanagenesis events with anatomically accurate chromosome representation (p-arm, centromere, q-arm)
- **SV Density Heatmap**: Shows concentration of structural variants across the genome
- **Copy Number Profile**: Genome-wide CN landscape with color-coded gains/losses
- **Mechanism Distribution**: Pie chart showing relative proportions of each mechanism

#### gGnome-style Regional Visualization

Detailed view of chromoanagenesis regions with SV arcs and copy number profiles, inspired by gGnome gWalks visualization.

```R
# Visualize specific chromoanagenesis region
plots <- plot_chromoanagenesis_region(
    chromoanagenesis_result = results,
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    chrom = "1",                      # Chromosome of interest
    mechanism = "chromothripsis",     # or "chromoplexy", "chromosynthesis"
    include_ideogram = TRUE
)

# Display combined plot
print(plots)

# Or arrange multiple regions
combined <- arrange_regional_plots(
    list(plots1, plots2, plots3),
    ncol = 1
)
```

Regional plots show:
- **SV Arcs**: Curved lines connecting breakpoints with arc curvature indicating SV span
- **Copy Number Profile**: Detailed CN changes with intelligent Y-axis scaling
- **Regional Ideogram**: Chromosome context with highlighted region

#### Circos Plots for Multi-Chromosome Events

Circular genome-wide visualization ideal for chromoplexy and complex multi-chromosome events.

```R
# Generate circos plot
plot_chromoanagenesis_circos(
    chromoanagenesis_result = results,
    SV.sample = sv_data,
    CNV.sample = cnv_data,
    sample_name = "Sample_01",
    mechanisms = c("chromothripsis", "chromoplexy", "chromosynthesis")
)
```

Circos plot structure (from outer to inner):
1. **Chromosome Ideogram**: Chromosome labels and boundaries
2. **Copy Number Track**: Color-coded CN profile (blue = loss, red = gain, gray = neutral)
3. **Chromoanagenesis Regions**: Highlights detected events (red = chromothripsis, blue = chromoplexy, green = chromosynthesis)
4. **SV Links**: Arcs connecting breakpoints, with different colors for SV types (DEL, DUP, INV, TRA)

**Note**: Circos plots require the `circlize` package:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("circlize")
```

## Important Notes

### Data Preparation

1. **Minimum cluster size**: The `min.Size` parameter in `detect_chromothripsis()` defaults to 1. This ensures detection of chromothripsis events that might involve linked clusters across multiple chromosomes.

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

### Primary Citation

If you use **ShatterSeek Extended** in your research, please cite:

**For the extended edition (chromoplexy, chromosynthesis, visualization, VCF support):**
> ShatterSeek Extended Edition (v2.0+): Comprehensive chromoanagenesis detection and analysis.
> GitHub: https://github.com/godkin1211/ShatterSeek
> DOI: [To be assigned upon publication]

**For the original chromothripsis detection framework:**
> Cortes-Ciriano, I., Lee, J.J., Xi, R. et al. Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat Genet* **52**, 331–341 (2020). https://doi.org/10.1038/s41588-019-0576-7

### Additional Citations

For specific features, please also cite the relevant method papers:

- **Chromoplexy concept**: Baca, S.C. et al. Punctuated evolution of prostate cancer genomes. *Cell* **153**(3), 666-677 (2013). https://doi.org/10.1016/j.cell.2013.03.021
- **Chromosynthesis/FoSTeS**: Liu, P. et al. Chromosome catastrophes involve replication mechanisms generating complex genomic rearrangements. *Cell* **146**(6), 889-903 (2011). https://doi.org/10.1016/j.cell.2011.07.042
- **Breakpoint repair mechanisms**: Tubio, J.M.C. & Estivill, X. Cancer: When catastrophe strikes a cell. *Nature* **470**, 476-477 (2011). https://doi.org/10.1038/470476a

### Recommended Citation Format

**For publications using the extended features:**
```
We performed chromoanagenesis analysis using ShatterSeek Extended Edition
(v2.0+), which extends the original chromothripsis detection framework
(Cortes-Ciriano et al., 2020) with chromoplexy and chromosynthesis detection,
advanced visualization, and direct VCF support.
```

**For publications using only chromothripsis detection:**
```
We performed chromothripsis analysis using ShatterSeek (Cortes-Ciriano et al., 2020).
```

## License

### Software License

**ShatterSeek Extended Edition** is released under the **GNU General Public License v3.0 or later (GPL-3+)**.

This means you are free to:
- ✓ Use the software for any purpose (academic or commercial)
- ✓ Study and modify the source code
- ✓ Distribute copies of the software
- ✓ Distribute modified versions

**Under the following terms:**
- Source code must be made available when distributing the software
- Modifications must be released under the same license
- Changes made to the code must be documented
- No warranty is provided

See the full license text at: https://www.gnu.org/licenses/gpl-3.0.html

### Original Framework License Notice

The original ShatterSeek chromothripsis detection framework (versions < 2.0) was developed at Harvard University and is subject to additional licensing terms for commercial use.

**For commercial applications of the original chromothripsis detection framework**, please contact:
- **Dr. Sonalee Barthakur**
- Harvard University Office of Technology Development
- Email: hms_otd@harvard.edu

**Note:** The extended features (chromoplexy detection, chromosynthesis detection, advanced visualization, VCF parsing, breakpoint analysis) added in version 2.0+ are independent implementations under GPL-3+ and do not require additional commercial licensing.

### Academic Use

ShatterSeek Extended is **strongly recommended for academic and research use**. We encourage researchers to:
- Contribute improvements back to the community
- Report bugs and suggest features via GitHub Issues
- Cite both the original framework and extended edition in publications
- Share analysis pipelines and custom workflows

### Disclaimer

THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY ARISING FROM THE USE OF THE SOFTWARE.

## Contact

### For ShatterSeek Extended Edition (v2.0+)

**Technical issues and feature requests:**
- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues

**Collaboration and contributions:**
- Email: n28111021@gs.ncku.edu.tw

### For Original ShatterSeek Framework

**Scientific questions about chromothripsis detection:**
- Isidro Cortes Ciriano: isidrolauscher@gmail.com or icortes@ebi.ac.uk
- Peter J Park: peter_park@hms.harvard.edu

**Commercial licensing:**
- Harvard University Office of Technology Development: hms_otd@harvard.edu

## Acknowledgments

### Extended Edition Development

**ShatterSeek Extended Edition (v2.0+)** represents a significant expansion of the original framework, adding:
- Chromoplexy and chromosynthesis detection modules
- Advanced visualization system (genome dashboards, gGnome-style plots, circos plots)
- Direct VCF input with robust BND/translocation parsing
- Breakpoint sequence analysis and repair mechanism inference
- Integrated multi-mechanism classification

Development of the extended edition was motivated by the need for comprehensive chromoanagenesis analysis tools that go beyond single-mechanism detection and provide production-ready VCF support for modern clinical pipelines.

### Original Framework

The original ShatterSeek chromothripsis detection framework was developed by:
- **Isidro Cortes-Ciriano** (European Bioinformatics Institute, EMBL-EBI)
- **Ruibin Xi** (Peking University)
- **Peter J. Park** (Harvard Medical School)

And validated using data from:
- The Pan-Cancer Analysis of Whole Genomes (PCAWG) Consortium
- ~2,600 whole-genome sequencing datasets across 38 tumor types

### Community

We thank the genomics community for feedback, bug reports, and feature requests that have shaped the extended edition. Special thanks to users who provided VCF format examples from diverse variant callers, enabling robust multi-caller support.

---

**Current Version**: 2.0.0 (Extended Edition)
**Release Date**: January 2025
**Original ShatterSeek**: Cortes-Ciriano et al., 2020
**Extended Features**: Chromoplexy, chromosynthesis, advanced visualization, VCF support

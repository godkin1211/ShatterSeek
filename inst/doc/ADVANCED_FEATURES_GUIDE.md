# ShatterSeek Extended Edition - Advanced Features Guide

## Table of Contents
1. [Integrated Classifier for Mixed Mechanisms](#integrated-classifier)
2. [Breakpoint Sequence Analysis](#breakpoint-sequence-analysis)
3. [Workflow Examples](#workflow-examples)
4. [Interpretation Guidelines](#interpretation-guidelines)

---

## Integrated Classifier for Mixed Mechanisms

### Overview

The integrated classifier identifies and characterizes **mixed chromoanagenesis events**, where multiple catastrophic mechanisms co-occur in the same sample. This is clinically important as mixed mechanisms often indicate:

- More severe genomic instability
- Multiple catastrophic events
- Higher therapeutic resistance
- Poorer prognosis

### Key Features

#### 1. Spatial Overlap Detection
Identifies chromosomes where multiple mechanisms operate on the same region:
- Chromothripsis + Chromoplexy
- Chromothripsis + Chromosynthesis
- All three mechanisms simultaneously

#### 2. Chromosome-Level Classification
For each chromosome, determines:
- Which mechanisms are present
- Dominant mechanism
- Whether mechanisms overlap spatially

#### 3. Sample-Level Classification
Overall classification into categories:
- **Pure mechanism**: Single mechanism type only
- **Mixed independent**: Multiple mechanisms on different chromosomes
- **Mixed overlapping**: Mechanisms co-occurring on same chromosomes

#### 4. Complexity Scoring
Multi-component score (0-1 scale) considering:
- **Mechanism diversity** (30%): Number of different mechanisms
- **Spatial overlap** (25%): Number of overlapping events
- **Chromosome spread** (25%): Number of affected chromosomes
- **Event count** (20%): Total catastrophic events

Complexity levels:
- **Low** (< 0.3): Simple, single-mechanism events
- **Moderate** (0.3-0.6): Multiple events or limited overlap
- **High** (0.6-0.8): Extensive genomic chaos
- **Very High** (> 0.8): Extreme complexity, multiple overlapping mechanisms

### Usage

```r
library(ShatterSeek)

# Step 1: Run comprehensive analysis
results <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data,
    genome = "hg19"
)

# Step 2: Classify mixed mechanisms
mixed_class <- classify_mixed_mechanisms(
    chromoanagenesis_result = results,
    overlap_threshold = 1e6,  # 1 Mb minimum overlap
    min_confidence = 0.3      # Include events with ≥30% confidence
)

# Step 3: Examine results
print(mixed_class)
summary(mixed_class)
```

### Visualization Functions

#### 1. `plot_mechanism_landscape()`
Genome-wide view showing mechanism distribution across chromosomes.

**Features:**
- Color-coded mechanisms
- Confidence scores as point sizes
- Easy identification of affected chromosomes

```r
p <- plot_mechanism_landscape(
    mixed_mechanisms_result = mixed_class,
    sample_name = "Sample_01",
    show_confidence = TRUE
)
print(p)
```

#### 2. `plot_chromosome_mechanisms()`
Heatmap showing presence/absence of each mechanism per chromosome.

**Features:**
- Red borders highlight chromosomes with mixed mechanisms
- Quick overview of mechanism distribution
- Identifies chromosomes requiring detailed analysis

```r
p <- plot_chromosome_mechanisms(
    mixed_mechanisms_result = mixed_class,
    sample_name = "Sample_01"
)
print(p)
```

#### 3. `plot_complexity_breakdown()`
Bar chart decomposing the complexity score into components.

**Features:**
- Shows contribution of each component
- Helps understand what drives complexity
- Raw vs weighted scores

```r
p <- plot_complexity_breakdown(
    mixed_mechanisms_result = mixed_class,
    sample_name = "Sample_01"
)
print(p)
```

#### 4. `plot_mechanism_dominance()`
Pie chart showing proportion of each mechanism.

**Features:**
- Visual representation of mechanism balance
- Identifies dominant mechanisms
- Useful for quick assessment

```r
p <- plot_mechanism_dominance(
    mixed_mechanisms_result = mixed_class,
    sample_name = "Sample_01"
)
print(p)
```

#### 5. `plot_mechanism_report()`
Comprehensive report combining all visualizations.

```r
report <- plot_mechanism_report(
    mixed_mechanisms_result = mixed_class,
    sample_name = "Sample_01"
)

# Save to file
ggsave("mechanism_report.pdf", report, width = 14, height = 10)
```

### Result Interpretation

#### Sample Classification Categories

**1. No chromoanagenesis**
- No catastrophic events detected
- Normal or simple rearrangements only

**2. Pure [mechanism] only**
- Single mechanism type
- Examples: "Pure chromothripsis only", "Pure chromoplexy only"
- Interpretation: Single catastrophic event

**3. Multiple independent mechanisms**
- Different mechanisms on different chromosomes
- No spatial overlap
- Interpretation: Sequential or independent catastrophic events

**4. Mixed mechanisms with spatial overlap**
- Multiple mechanisms on same chromosome/region
- Spatial overlap detected
- Interpretation: Complex, overlapping catastrophic events
- Clinical significance: Highest complexity, poorest prognosis

#### Complexity Interpretation

| Level | Score | Interpretation | Clinical Implications |
|-------|-------|----------------|----------------------|
| Low | < 0.3 | Simple event | Standard treatment may work |
| Moderate | 0.3-0.6 | Moderate complexity | Consider combination therapy |
| High | 0.6-0.8 | Extensive chaos | High risk, intensive monitoring |
| Very High | > 0.8 | Extreme complexity | Experimental approaches needed |

---

## Breakpoint Sequence Analysis

### Overview

Breakpoint sequence analysis examines the DNA sequences at structural variant junctions to infer the **DNA repair mechanisms** that created them. Different repair pathways leave characteristic signatures:

- **NHEJ** (Non-Homologous End Joining): Blunt ends, no microhomology
- **MMEJ** (Microhomology-Mediated End Joining): 2-25 bp microhomology
- **HR** (Homologous Recombination): Long homology (>25 bp)
- **FoSTeS/MMBIR**: Serial replication with microhomology, associated with chromosynthesis

### Requirements

To run breakpoint sequence analysis, you need a **reference genome**:

```r
# Option 1: Use BSgenome package (recommended)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# Option 2: Use FASTA file
genome <- "/path/to/reference/genome.fa"
```

### Usage

```r
# Analyze breakpoint sequences
bp_analysis <- analyze_breakpoint_sequences(
    SV.sample = SV_data,
    genome = genome,
    flank_size = 50,           # Extract 50 bp around breakpoints
    min_microhomology = 2,     # Minimum 2 bp to report
    max_microhomology = 25     # Maximum search length
)

# View results
print(bp_analysis)
summary(bp_analysis)
```

### What Gets Analyzed

For each structural variant breakpoint:

1. **Sequence Extraction**
   - Flanking sequences around both breakpoints
   - Respect strand orientation
   - Default: 50 bp each side

2. **Microhomology Detection**
   - Search for homologous sequences at junction
   - Range: 2-25 bp (configurable)
   - Records sequence and position

3. **Insertion Detection**
   - Non-templated insertions
   - Templated insertions
   - Classification of insertion type

4. **Repair Mechanism Classification**
   - Based on sequence features
   - Confidence level assigned
   - Evidence documented

### Output Structure

```r
bp_analysis$sequences          # Extracted sequences
bp_analysis$microhomology      # Microhomology results
bp_analysis$insertions         # Insertion analysis
bp_analysis$repair_mechanisms  # Classified mechanisms
bp_analysis$summary           # Summary statistics
```

### Visualization Functions

#### 1. `plot_repair_mechanisms()`
Distribution of inferred DNA repair mechanisms.

```r
p <- plot_repair_mechanisms(
    breakpoint_result = bp_analysis,
    sample_name = "Sample_01"
)
print(p)
```

**Interpretation:**
- High NHEJ: Standard DSB repair, common in normal cells
- High MMEJ: Alternative repair, common in cancer
- High FoSTeS/MMBIR: Serial replication, indicates chromosynthesis

#### 2. `plot_microhomology_distribution()`
Histogram of microhomology lengths.

```r
p <- plot_microhomology_distribution(
    breakpoint_result = bp_analysis,
    sample_name = "Sample_01"
)
print(p)
```

**Interpretation:**
- Peak at 2-5 bp: Classic MMEJ
- Peak at 5-10 bp: FoSTeS/MMBIR
- Uniform distribution: Mixed repair processes

#### 3. `plot_repair_by_svtype()`
Relationship between SV types and repair mechanisms.

```r
p <- plot_repair_by_svtype(
    breakpoint_result = bp_analysis,
    sv_data = SV_data,
    sample_name = "Sample_01"
)
print(p)
```

**Expected patterns:**
- Deletions: Often NHEJ or MMEJ
- Duplications: FoSTeS/MMBIR (especially tandem)
- Translocations: NHEJ or MMEJ
- Inversions: NHEJ

#### 4. `plot_breakpoint_report()`
Comprehensive report combining all breakpoint visualizations.

```r
report <- plot_breakpoint_report(
    breakpoint_result = bp_analysis,
    sv_data = SV_data,
    sample_name = "Sample_01"
)

ggsave("breakpoint_report.pdf", report, width = 12, height = 10)
```

### Repair Mechanism Interpretation

#### NHEJ (Non-Homologous End Joining)
- **Signature**: No or minimal microhomology (0-1 bp)
- **Process**: Direct ligation of broken ends
- **Accuracy**: Error-prone, can cause small indels
- **Context**: Primary DSB repair in mammalian cells
- **Clinical**: Normal repair, not specific to cancer

#### MMEJ/Alt-EJ (Microhomology-Mediated End Joining)
- **Signature**: 2-25 bp microhomology
- **Process**: End resection, annealing at microhomology
- **Accuracy**: Error-prone, causes deletions
- **Context**: Backup pathway, upregulated in cancer
- **Clinical**: Associated with genomic instability

#### FoSTeS/MMBIR (Fork Stalling and Template Switching / Microhomology-Mediated Break-Induced Replication)
- **Signature**: Microhomology + tandem duplications
- **Process**: Serial template switching during replication
- **Accuracy**: Creates complex rearrangements
- **Context**: Chromosynthesis mechanism
- **Clinical**: Indicates replication stress

#### HR-like (Homologous Recombination)
- **Signature**: Long homology (>25 bp)
- **Process**: Template-directed accurate repair
- **Accuracy**: High fidelity
- **Context**: Rare in structural variants
- **Clinical**: May indicate HR pathway activation

---

## Workflow Examples

### Example 1: Complete Analysis of Single Sample

```r
library(ShatterSeek)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load data
data(DO17373)

# Prepare data objects
SV_data <- SVs(chrom1 = SV_DO17373$chrom1, ...)
CN_data <- CNVsegs(chrom = SCNA_DO17373$chromosome, ...)

# Step 1: Comprehensive chromoanagenesis detection
chromoanag <- detect_chromoanagenesis(
    SV.sample = SV_data,
    CNV.sample = CN_data
)

# Step 2: Mixed mechanism classification
mixed_class <- classify_mixed_mechanisms(chromoanag)

# Step 3: Breakpoint sequence analysis
bp_analysis <- analyze_breakpoint_sequences(
    SV.sample = SV_data,
    genome = BSgenome.Hsapiens.UCSC.hg19
)

# Step 4: Create comprehensive reports
plot_mechanism_report(mixed_class, "Sample_01")
plot_breakpoint_report(bp_analysis, SV_data, "Sample_01")

# Step 5: Integrate findings
cat("\n=== Integrated Results ===\n")
cat("Chromoanagenesis:", mixed_class$sample_classification$classification, "\n")
cat("Complexity:", mixed_class$complexity$complexity_level, "\n")
cat("Dominant repair:", bp_analysis$summary$dominant_mechanism, "\n")
```

### Example 2: Batch Processing Multiple Samples

```r
# Define sample list
samples <- list(
    Sample_01 = list(SV = SV_data_1, CNV = CN_data_1),
    Sample_02 = list(SV = SV_data_2, CNV = CN_data_2),
    Sample_03 = list(SV = SV_data_3, CNV = CN_data_3)
)

# Batch analysis
results <- lapply(names(samples), function(sample_name) {
    cat("\nProcessing:", sample_name, "\n")

    # Chromoanagenesis
    chromoanag <- detect_chromoanagenesis(
        SV.sample = samples[[sample_name]]$SV,
        CNV.sample = samples[[sample_name]]$CNV,
        verbose = FALSE
    )

    # Mixed mechanisms
    mixed <- classify_mixed_mechanisms(chromoanag)

    # Breakpoints (if genome available)
    # bp <- analyze_breakpoint_sequences(
    #     SV.sample = samples[[sample_name]]$SV,
    #     genome = genome
    # )

    list(
        chromoanag = chromoanag,
        mixed = mixed,
        complexity = mixed$complexity$complexity_level
    )
})
names(results) <- names(samples)

# Summary table
summary_df <- data.frame(
    Sample = names(results),
    Classification = sapply(results, function(x) x$mixed$sample_classification$classification),
    Complexity = sapply(results, function(x) x$complexity),
    N_Mechanisms = sapply(results, function(x) x$mixed$sample_classification$n_mechanisms)
)
print(summary_df)
```

### Example 3: Focused Analysis on Mixed Chromosomes

```r
# Identify chromosomes with mixed mechanisms
mixed_chroms <- mixed_class$chromosome_classification[
    mixed_class$chromosome_classification$is_mixed,
]

if (nrow(mixed_chroms) > 0) {
    cat("Analyzing mixed mechanism chromosomes:\n")

    for (i in 1:nrow(mixed_chroms)) {
        chr <- mixed_chroms$chrom[i]
        mechs <- mixed_chroms$mechanisms[i]

        cat(sprintf("\n%s: %s\n", chr, mechs))

        # Filter SVs on this chromosome
        chr_svs <- SV_data[SV_data@chrom1 == chr | SV_data@chrom2 == chr, ]

        # Detailed breakpoint analysis for this chromosome
        chr_bp <- analyze_breakpoint_sequences(
            SV.sample = chr_svs,
            genome = genome
        )

        # Check repair mechanisms
        print(table(chr_bp$repair_mechanisms$repair_mechanism))
    }
}
```

---

## Interpretation Guidelines

### Clinical Decision Framework

#### 1. Assess Overall Pattern

**Question**: What type of chromoanagenesis is present?

- **Pure chromothripsis**: Single catastrophic event, localized
  - *Action*: Consider standard targeted therapy
  - *Prognosis*: Moderate

- **Pure chromoplexy**: Chained translocations
  - *Action*: Check for specific cancer types (prostate, etc.)
  - *Prognosis*: Depends on cancer type

- **Mixed mechanisms**: Multiple catastrophic events
  - *Action*: Aggressive treatment, close monitoring
  - *Prognosis*: Generally poor

#### 2. Evaluate Complexity

**Question**: How complex is the genomic chaos?

| Complexity | Genomic Impact | Therapeutic Strategy |
|------------|---------------|---------------------|
| Low | Limited | Standard protocols |
| Moderate | Intermediate | Enhanced monitoring |
| High | Severe | Combination therapy |
| Very High | Extreme | Experimental approaches |

#### 3. Examine Repair Mechanisms

**Question**: Which repair pathways are active?

- **High NHEJ**: Normal repair, not cancer-specific
- **High MMEJ**: Backup pathway activation, instability
- **High FoSTeS/MMBIR**: Replication stress, chromosynthesis
- **Mixed mechanisms**: Multiple pathway dysregulation

**Therapeutic implications**:
- MMEJ-dominant: Consider PARP inhibitors
- HR-deficient: Platinum-based chemotherapy
- FoSTeS/MMBIR: Replication stress inhibitors

#### 4. Integration

Combine findings for comprehensive assessment:

```
Sample Assessment Framework:

1. Mechanism Type: [chromothripsis/chromoplexy/mixed]
2. Spatial Pattern: [localized/distributed/overlapping]
3. Complexity: [Low/Moderate/High/Very High]
4. Dominant Repair: [NHEJ/MMEJ/FoSTeS/mixed]
5. Clinical Risk: [Low/Intermediate/High/Very High]

Risk Calculation:
- Pure + Low complexity + NHEJ = Low risk
- Mixed + High complexity + MMEJ = High risk
- Overlapping + Very High + FoSTeS = Very High risk
```

### Research Applications

#### 1. Mechanistic Studies
- Compare repair mechanisms across cancer types
- Identify pathway dysregulation patterns
- Study temporal evolution of chromoanagenesis

#### 2. Therapeutic Response
- Correlate mechanism type with treatment response
- Identify biomarkers for pathway inhibitors
- Predict therapy resistance

#### 3. Evolutionary Analysis
- Track mechanism changes over time
- Identify clonal vs subclonal events
- Study progression patterns

---

## Troubleshooting

### Common Issues

**1. No reference genome available**
```r
# Error: Cannot access reference genome
# Solution: Install BSgenome package
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
```

**2. Low confidence scores**
```r
# Issue: Most events have low confidence
# Solution: Adjust confidence threshold
mixed_class <- classify_mixed_mechanisms(
    chromoanag,
    min_confidence = 0.2  # Lower threshold
)
```

**3. No mixed mechanisms detected**
```r
# Issue: All mechanisms are independent
# Solution: Adjust overlap threshold
mixed_class <- classify_mixed_mechanisms(
    chromoanag,
    overlap_threshold = 5e5  # Smaller threshold (500 kb)
)
```

**4. Memory issues with large datasets**
```r
# Solution: Process chromosomes separately
for (chr in unique(SV_data@chrom1)) {
    chr_sv <- SV_data[SV_data@chrom1 == chr, ]
    # Analyze chr_sv
}
```

---

## References

### Chromoanagenesis Mechanisms

1. **Chromothripsis**: Stephens et al., Cell 2011
2. **Chromoplexy**: Baca et al., Cell 2013
3. **Chromosynthesis**: Liu et al., Cell 2011

### DNA Repair Mechanisms

1. **NHEJ**: Lieber, Annu Rev Biochem 2010
2. **MMEJ**: McVey & Lee, Trends Genet 2008
3. **FoSTeS/MMBIR**: Lee et al., Cell 2007; Hastings et al., Nat Rev Genet 2009

### Clinical Implications

1. **Cancer Genomics**: Cortés-Ciriano et al., Nat Genet 2020
2. **Therapeutic Response**: Zhang et al., Nat Med 2015

---

## Support

For questions or issues:
- GitHub: https://github.com/parklab/ShatterSeek
- Email: Support contact information
- Documentation: Run `?function_name` in R for detailed help

---

**Version**: 1.2.0
**Last Updated**: 2025
**Authors**: ShatterSeek Development Team

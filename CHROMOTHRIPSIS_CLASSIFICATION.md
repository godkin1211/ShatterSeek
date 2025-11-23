# Chromothripsis Classification Criteria

## Overview

ShatterSeek Extended Edition implements the chromothripsis classification criteria from the original ShatterSeek publication (Cortes-Ciriano et al., Nature Genetics 2020), based on the framework proposed by Korbel and Campbell (Cell 2013).

## Classification Levels

Chromothripsis events are classified into three categories:

1. **High confidence chromothripsis**
2. **Low confidence chromothripsis**
3. **Not chromothripsis**

## Detailed Criteria

### High Confidence Chromothripsis (Type 1: Single/Few Chromosomes)

A chromothripsis event is classified as **high confidence** if it meets ALL of the following criteria:

- **≥ 6 interleaved intrachromosomal SVs**
- **≥ 7 contiguous CN segments oscillating between 2 states**
- **Fragment joins test passed** (p-value > 0.05, indicating random fragment joining consistent with NHEJ repair)
- **Clustering test passed**: EITHER
  - Chromosomal breakpoint enrichment (p < 0.05), OR
  - Exponential distribution of breakpoints test (p < 0.05)

### High Confidence Chromothripsis (Type 2: Multi-Chromosomal)

A chromothripsis event is classified as **high confidence** if it meets ALL of the following criteria:

- **≥ 3 interleaved intrachromosomal SVs**
- **≥ 4 interchromosomal SVs (translocations)**
- **≥ 7 contiguous CN segments oscillating between 2 states**
- **Fragment joins test passed** (p-value > 0.05)

### Low Confidence Chromothripsis

A chromothripsis event is classified as **low confidence** if it meets ALL of the following criteria:

- **≥ 6 interleaved intrachromosomal SVs**
- **4-6 contiguous CN segments oscillating between 2 states** (note the upper bound!)
- **Fragment joins test passed** (p-value > 0.05)
- **Clustering test passed**: EITHER
  - Chromosomal breakpoint enrichment (p < 0.05), OR
  - Exponential distribution of breakpoints test (p < 0.05)

### Not Chromothripsis

Events that do not meet the criteria for high or low confidence classification.

## Key Statistical Tests

### 1. Fragment Joins Test (Chi-squared test)

**Purpose**: Tests whether the distribution of SV types (DEL, DUP, h2hINV, t2tINV) deviates from a uniform distribution.

**Interpretation**:
- **p-value > 0.05**: Consistent with random fragment joining (PASS - characteristic of chromothripsis)
- **p-value ≤ 0.05**: Non-random fragment joining (FAIL - not chromothripsis)

**Biological meaning**: Chromothripsis involves random stitching of shattered DNA fragments through non-homologous end joining (NHEJ), resulting in approximately equal proportions of different SV types.

### 2. Chromosomal Breakpoint Enrichment Test (Binomial test)

**Purpose**: Tests whether a chromosome is significantly enriched for breakpoints compared to the genome-wide background.

**Interpretation**:
- **p-value < 0.05**: Chromosome is enriched for breakpoints (PASS)
- **p-value ≥ 0.05**: No significant enrichment (FAIL)

**Biological meaning**: Massive chromothripsis events concentrate breakpoints in specific chromosomes against a quiescent genomic background.

### 3. Exponential Distribution of Breakpoints Test (Kolmogorov-Smirnov test)

**Purpose**: Tests whether breakpoints are clustered (non-exponentially distributed).

**Interpretation**:
- **p-value < 0.05**: Breakpoints are clustered (PASS)
- **p-value ≥ 0.05**: Breakpoints follow exponential distribution (FAIL)

**Biological meaning**: Chromothripsis breakpoints are spatially clustered in localized genomic regions.

## CN Oscillation Counting

**2-state oscillation**: Alternating pattern between exactly 2 different copy number states (e.g., 1 → 2 → 1 → 2).

**3-state oscillation**: Pattern involving 3 copy number states where adjacent segments differ by at most 1 (e.g., 1 → 2 → 1 → 2 → 3).

**Counting method**: The maximum number of **contiguous** segments that satisfy the oscillation pattern.

Example:
```
CN: 2 1 2 1 2 1 3 3 2 1 2
     ↑─────────↑           (6 segments oscillating between CN=1 and CN=2)
                   ↑───↑   (4 segments: 3→3→2→1→2)

Maximum 2-state oscillation: 6 segments
```

## Important Implementation Notes

### Fragment Joins Test is CRITICAL

**This test was previously missing from the classification logic!**

The fragment joins test is **mandatory** for both high and low confidence classification. A p-value > 0.05 indicates random fragment joining, which is a hallmark of chromothripsis.

Without this test, events with non-random SV type distributions (e.g., biased toward duplications due to replicative mechanisms) could be incorrectly classified as chromothripsis.

### Upper Bound on CN Oscillations for Low Confidence

Low confidence chromothripsis requires **4-6 oscillating segments** (with upper bound).

Events with ≥7 oscillating segments that meet other criteria should be classified as **high confidence**, not low confidence.

This distinction is important for:
1. Proper stratification of event severity
2. Clinical interpretation (high confidence events may have different prognostic implications)

### Type 1 vs Type 2 High Confidence

**Type 1** (Single chromosome):
- Localized catastrophic shattering
- Many intrachromosomal SVs (≥6)
- Requires clustering test (enrichment OR exponential distribution)

**Type 2** (Multi-chromosomal):
- Involves multiple chromosomes
- Fewer intrachromosomal SVs needed (≥3)
- Requires many translocations (≥4)
- Does NOT require clustering test (translocations themselves indicate multi-chromosome involvement)

## Usage Example

```R
library(ShatterSeek)

# Run chromothripsis detection
chromoth_result <- detect_chromothripsis(SV_data, CN_data, genome = "hg19")

# Classify events
classification <- classify_chromothripsis(chromoth_result)

# View results
print(classification)

# Generate summary report
summarize_chromothripsis(chromoth_result)
```

Expected output:
```
         CHROMOTHRIPSIS DETECTION SUMMARY
======================================================================

Total chromosomes with SV clusters: 3

Classification Results (Cortes-Ciriano et al. 2020 criteria):
  - High confidence: 1 chromosome(s)
  - Low confidence:  1 chromosome(s)
  - Not chromothripsis: 1 chromosome(s)

High-confidence chromothripsis chromosomes:
  2

Low-confidence chromothripsis chromosomes:
  3

Detailed Results:
Chrom    Classification       CN Osc.    Intra SVs  TRA      FragJoin p
----------------------------------------------------------------------
2        High confidence      11         8          7        0.865
3        Low confidence       5          6          2        0.572

Note: Fragment joins p-value > 0.05 indicates random fragment joining
      (characteristic of chromothripsis)
======================================================================
```

## References

1. **Korbel, J. O. & Campbell, P. J.** Criteria for inference of chromothripsis in cancer genomes. *Cell* **152**, 1226-1236 (2013).

2. **Cortes-Ciriano, I. et al.** Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat. Genet.* **52**, 331-341 (2020).

3. **Stephens, P. J. et al.** Massive genomic rearrangement acquired in a single catastrophic event during cancer development. *Cell* **144**, 27-40 (2011).

## Change Log

### Version 2.0.1 (Current)

- **FIXED**: Fragment joins test now correctly enforced in classification logic
- **CHANGED**: Classification levels renamed from "Likely/Possible/Unlikely" to "High confidence/Low confidence/Not chromothripsis" to match original publication
- **ADDED**: Upper bound (≤6) for CN oscillations in low confidence category
- **IMPROVED**: Clear distinction between Type 1 and Type 2 high confidence criteria
- **ENHANCED**: Documentation of statistical test interpretations

### Version 2.0.0

- Original Extended Edition release with integrated chromoanagenesis detection

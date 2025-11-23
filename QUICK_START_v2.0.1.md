# Quick Start Guide - ShatterSeek v2.0.1

## What Changed?

### Classification System Update

**Previous (v2.0.0)**:
- Likely chromothripsis
- Possible chromothripsis
- Unlikely

**Current (v2.0.1)**:
- **High confidence** chromothripsis
- **Low confidence** chromothripsis
- **Not chromothripsis**

### Key Improvements

1. **Fragment joins test now enforced** (p > 0.05 required)
2. **CN oscillation upper bound** for low confidence (4-6 segments)
3. **Fully compliant** with Cortes-Ciriano et al. 2020 criteria

## Installation

After pulling the latest changes:

```bash
cd /Users/godkin/Projects/ShatterSeek
R CMD INSTALL .
```

Or from R:

```R
devtools::install()
```

## Basic Usage

```R
library(ShatterSeek)

# Load example data
data(DO17373)

# Prepare SV data
SV_data <- SVs(
    chrom1 = as.character(SV_DO17373$chrom1),
    pos1 = as.numeric(SV_DO17373$start1),
    chrom2 = as.character(SV_DO17373$chrom2),
    pos2 = as.numeric(SV_DO17373$end2),
    SVtype = as.character(SV_DO17373$svclass),
    strand1 = as.character(SV_DO17373$strand1),
    strand2 = as.character(SV_DO17373$strand2)
)

# Prepare CNV data
CN_data <- CNVsegs(
    chrom = as.character(SCNA_DO17373$chromosome),
    start = SCNA_DO17373$start,
    end = SCNA_DO17373$end,
    total_cn = SCNA_DO17373$total_cn
)

# Run detection
chromoth_result <- detect_chromothripsis(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg19"
)

# Classify events (NEW - no parameters needed)
classification <- classify_chromothripsis(chromoth_result)

# View results
print(classification)

# Generate detailed report
summarize_chromothripsis(chromoth_result)
```

## Example Output

```
         CHROMOTHRIPSIS DETECTION SUMMARY
======================================================================

Total chromosomes with SV clusters: 8

Classification Results (Cortes-Ciriano et al. 2020 criteria):
  - High confidence: 4 chromosome(s)
  - Low confidence:  0 chromosome(s)
  - Not chromothripsis: 4 chromosome(s)

High-confidence chromothripsis chromosomes:
  2, 21, 3, X

Detailed Results:
Chrom    Classification       CN Osc.    Intra SVs  TRA      FragJoin p
----------------------------------------------------------------------
2        High confidence      13         15         7        0.865
21       High confidence      70         785        602      0.855
3        High confidence      14         10         6        0.572
X        High confidence      624        595        1        0.698

Note: Fragment joins p-value > 0.05 indicates random fragment joining
      (characteristic of chromothripsis)
======================================================================
```

## Understanding the Classification

### High Confidence Chromothripsis

**Type 1** (Single/few chromosomes):
- ✓ ≥6 intrachromosomal SVs
- ✓ ≥7 CN oscillating segments
- ✓ Fragment joins p > 0.05 (random joining)
- ✓ Clustering test passed

**Type 2** (Multi-chromosomal):
- ✓ ≥3 intrachromosomal SVs
- ✓ ≥4 translocations
- ✓ ≥7 CN oscillating segments
- ✓ Fragment joins p > 0.05

### Low Confidence Chromothripsis

- ✓ ≥6 intrachromosomal SVs
- ✓ 4-6 CN oscillating segments (note the upper limit!)
- ✓ Fragment joins p > 0.05
- ✓ Clustering test passed

### Not Chromothripsis

Events that don't meet the above criteria.

## Key Statistical Tests

### Fragment Joins Test (Chi-squared)

**What it tests**: Are SV types (DEL, DUP, INV) randomly distributed?

**Interpretation**:
- **p > 0.05**: Random (PASS - characteristic of chromothripsis) ✓
- **p ≤ 0.05**: Non-random (FAIL - likely other mechanism) ✗

**Why it matters**: Chromothripsis involves random stitching through NHEJ, producing ~equal proportions of SV types.

### Chromosomal Enrichment Test (Binomial)

**What it tests**: Is a chromosome enriched for breakpoints?

**Interpretation**:
- **p < 0.05**: Enriched (PASS) ✓
- **p ≥ 0.05**: Not enriched (FAIL) ✗

### Exponential Distribution Test (KS test)

**What it tests**: Are breakpoints clustered?

**Interpretation**:
- **p < 0.05**: Clustered (PASS) ✓
- **p ≥ 0.05**: Not clustered (FAIL) ✗

## Testing Your Installation

Run the validation script:

```bash
cd tests/validation
Rscript test_new_classification.R
```

You should see:
- ✓ All classifications use valid categories
- ✓ All chromothripsis events pass fragment joins test
- ✓ CN oscillation requirements validated
- ✓ Test completed successfully

## Troubleshooting

### Error: "Invalid classification categories found"

**Cause**: Package not reinstalled after code changes.

**Solution**:
```bash
R CMD INSTALL .
```

### Old classification names still appearing

**Cause**: Old package version loaded in R session.

**Solution**:
```R
# Restart R session
.rs.restartR()  # RStudio
# OR quit and relaunch R

# Then reload
library(ShatterSeek)
```

### Fragment joins test always fails

**Check**: Are you using the correct SV type encoding?
- DEL: +/-
- DUP: -/+
- h2hINV: +/+
- t2tINV: -/-

## Migration from v2.0.0

### Code Changes

**Old**:
```R
classification <- classify_chromothripsis(
    chromoth_output,
    min_cluster_size = 6,
    min_cn_oscillations = 4,
    max_pval_clustering = 0.05
)
```

**New**:
```R
# Parameters removed - uses tutorial criteria
classification <- classify_chromothripsis(chromoth_output)
```

### Output Changes

**Old field names**:
- "Likely chromothripsis"
- "Possible chromothripsis"
- "Unlikely"

**New field names**:
- "High confidence"
- "Low confidence"
- "Not chromothripsis"

## Further Reading

- **CHROMOTHRIPSIS_CLASSIFICATION.md**: Detailed classification criteria
- **CHANGES_v2.0.1.md**: Complete changelog
- **README.md**: General usage guide

## Support

- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues
- Email: n28111021@gs.ncku.edu.tw

## References

1. Cortes-Ciriano et al. (2020) Nat. Genet. 52:331-341
2. Korbel & Campbell (2013) Cell 152:1226-1236

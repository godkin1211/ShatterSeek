# Diagnostic Guide - Troubleshooting Classification Changes

## Overview

If your samples previously showed chromothripsis but now show "Not chromothripsis" with v2.0.1, this guide will help you diagnose why and understand what changed.

## What Changed in v2.0.1?

The major changes that could affect your samples:

1. **Fragment joins test now enforced** (CRITICAL)
   - Previous versions calculated this test but **did not enforce it**
   - Now MANDATORY: p-value > 0.05 required for chromothripsis classification
   - This is the most likely reason for missing detections

2. **CN oscillation upper bound for low confidence**
   - Low confidence now requires 4-6 oscillations (not just ≥4)
   - Events with ≥7 oscillations must be high confidence

3. **Stricter criteria alignment**
   - Full compliance with Cortes-Ciriano et al. 2020 criteria
   - All thresholds match original publication

## Quick Diagnostic Workflow

### Step 1: Install Diagnostic Tools

```bash
cd /Users/godkin/Projects/ShatterSeek
R CMD INSTALL .
```

### Step 2: Run Diagnostic on Your Sample

```R
library(ShatterSeek)

# Load your SV and CNV data
SV_data <- SVs(...)  # Your SV data
CN_data <- CNVsegs(...)  # Your CNV data

# Run detection
chromoth_result <- detect_chromothripsis(
    SV.sample = SV_data,
    seg.sample = CN_data,
    genome = "hg19"  # or "hg38"
)

# Classify
classification <- classify_chromothripsis(chromoth_result)

# Run diagnostic analysis
diagnostic <- diagnose_chromothripsis_classification(
    chromoth_result,
    show_all_chromosomes = FALSE  # Set TRUE to see all chromosomes
)

# Print detailed report
print_diagnostic_report(diagnostic)
```

### Step 3: Interpret Results

The diagnostic report shows:
- Which criteria each chromosome passes/fails
- Specific p-values for statistical tests
- Blocking reasons for "Not chromothripsis" classifications

## Common Issues and Solutions

### Issue 1: Fragment Joins Test Failing (Most Common)

**Symptom**: Diagnostic shows "Fragment joins FAILED (p≤0.05)"

**What this means**:
- Your SV types (DEL, DUP, INV) are NOT randomly distributed
- Chi-squared test detects bias toward specific SV types
- This suggests a replication-based mechanism (not chromothripsis)

**Why this is NEW**:
- Previous versions ignored this test
- v2.0.1 enforces it per original publication

**Solution options**:

1. **Accept the classification** - Your event may genuinely not be chromothripsis:
   - Consider chromosynthesis if you see tandem duplications
   - Check for CN gradients (incremental increases)
   - Look for FoSTeS/MMBIR signatures

2. **Check your SV type encoding**:
   ```R
   # Correct encoding:
   # DEL: strand1="+", strand2="-"
   # DUP: strand1="-", strand2="+"
   # h2hINV: strand1="+", strand2="+"
   # t2tINV: strand1="-", strand2="-"

   # Verify your data:
   table(SV_data@SVtype)
   table(paste(SV_data@strand1, SV_data@strand2))
   ```

3. **Understand the biology**:
   - Chromothripsis involves random NHEJ repair → equal SV type proportions
   - If your sample is biased (e.g., mostly DUPs), it's likely NOT chromothripsis

### Issue 2: CN Oscillations Insufficient

**Symptom**: "CN oscillations insufficient (X, need ≥4)" or "CN oscillations in low conf range 4-6, not high ≥7"

**What this means**:
- Not enough copy number changes
- Or borderline between low/high confidence

**Solutions**:

1. **Check your CNV data quality**:
   ```R
   # Are adjacent segments different?
   cn_df <- as.data.frame(CN_data)
   head(cn_df, 20)

   # Look for merged segments with same CN
   ```

2. **Verify CN calling**:
   - Ensure CNV caller detects fine-grained segments
   - Some callers merge small segments (reduces oscillation count)
   - Consider using a more sensitive caller

### Issue 3: Too Few Intrachromosomal SVs

**Symptom**: "Too few intrachromosomal SVs (X, need ≥6)"

**What this means**:
- Not enough clustered SVs on the chromosome
- Below the published threshold

**Solutions**:

1. **Check SV filtering**:
   - Are you filtering by size?
   - min_sv_size parameter may exclude valid SVs

2. **Verify SV caller sensitivity**:
   - Some callers miss smaller events
   - Consider using multiple callers and merging

### Issue 4: Clustering Tests Failed

**Symptom**: "Clustering tests FAILED (both p≥0.05)"

**What this means**:
- SVs not significantly clustered in space
- Chromosome not enriched for breakpoints

**This is expected for**:
- Scattered SVs across the genome
- Low SV counts
- Not true chromothripsis

## Advanced Diagnostic: Compare Old vs New

To understand what changed for your specific sample:

```R
# 1. Check the chromSummary for raw values
chromSummary <- chromoth_result@chromSummary

# Focus on chromosomes with clusters
chromSummary[chromSummary$clusterSize > 0, c(
    "chrom",
    "clusterSize",
    "max_number_oscillating_CN_segments_2_states",
    "pval_fragment_joins",
    "pval_exp_cluster",
    "chr_breakpoint_enrichment"
)]

# 2. Identify borderline cases
# Fragment joins close to 0.05?
borderline_frag <- chromSummary[
    !is.na(chromSummary$pval_fragment_joins) &
    chromSummary$pval_fragment_joins > 0.03 &
    chromSummary$pval_fragment_joins <= 0.05,
]

if (nrow(borderline_frag) > 0) {
    cat("BORDERLINE fragment joins p-values:\n")
    print(borderline_frag[, c("chrom", "pval_fragment_joins")])
    cat("\nThese would have been classified as chromothripsis in v2.0.0\n")
    cat("but fail in v2.0.1 due to fragment joins test enforcement.\n")
}
```

## Understanding Fragment Joins Test

### What It Tests

The fragment joins test (chi-squared test) checks if SV types are randomly distributed:

- **Expected for chromothripsis**: ~Equal proportions of DEL, DUP, INV
  - Random NHEJ repair produces all types equally

- **NOT expected for other mechanisms**:
  - Replication-based: Biased toward DUPs
  - BFB (breakage-fusion-bridge): Mostly INVs and DELs
  - FoSTeS/MMBIR: Tandem duplications

### How to Check Your SV Distribution

```R
# Extract SVs for a specific chromosome
chrom_svs <- SV_data[SV_data@chrom1 == "21" & SV_data@chrom2 == "21", ]

# Count SV types
sv_types <- table(chrom_svs@SVtype)
print(sv_types)

# Calculate proportions
prop.table(sv_types)

# Random distribution would be ~33% each for DEL, DUP, INV
# If you see 70% DUP, 20% DEL, 10% INV → NOT chromothripsis
```

## Validation: Test with Example Data

To verify the diagnostic tools work correctly:

```bash
cd tests/validation
Rscript diagnose_example.R
```

Expected output:
- 4 high confidence chromothripsis (chr 2, 3, 21, X)
- All pass fragment joins test (p > 0.05)
- All have ≥7 CN oscillations

## What If My Sample Is Genuinely Not Chromothripsis?

If the diagnostic shows your sample fails criteria, consider:

1. **Chromosynthesis**: Serial replication-based rearrangements
   - Biased toward tandem duplications
   - CN gradients (incremental increases)
   - Detection: `detect_chromosynthesis()`

2. **Chromoplexy**: Translocation chains
   - Multi-chromosome involvement
   - Stable copy number
   - Detection: `detect_chromoplexy()`

3. **Other complex rearrangements**:
   - BFB (breakage-fusion-bridge)
   - Localized hypermutation
   - Replication stress signatures

Run comprehensive analysis:
```R
results <- detect_chromoanagenesis(SV_data, CN_data)
```

## Getting Help

If you believe there's an error in the classification:

1. **Collect diagnostic output**:
   ```R
   # Save diagnostic report
   diagnostic <- diagnose_chromothripsis_classification(chromoth_result)
   write.csv(diagnostic, "my_sample_diagnostic.csv")
   ```

2. **Check raw data**:
   ```R
   # SV type distribution
   table(SV_data@SVtype)

   # Fragment joins p-values
   chromoth_result@chromSummary[, c("chrom", "pval_fragment_joins")]
   ```

3. **Report issue**:
   - GitHub: https://github.com/godkin1211/ShatterSeek/issues
   - Include diagnostic output
   - Mention your SV caller and CNV caller

## Summary Decision Tree

```
Sample previously showed chromothripsis, now doesn't
    |
    ├── Run diagnostic → Fragment joins p ≤ 0.05?
    |   ├── YES → Check SV type distribution
    |   |   ├── Biased toward DUP → Likely chromosynthesis
    |   |   └── Equal proportions → Check data encoding
    |   |
    |   ├── CN oscillations < 4?
    |   |   └── Check CNV data quality/caller sensitivity
    |   |
    |   ├── Too few SVs (<6)?
    |   |   └── Check SV filtering/caller sensitivity
    |   |
    |   └── Clustering failed?
    |       └── SVs not spatially clustered → Not chromothripsis
```

## References

1. **Cortes-Ciriano, I. et al.** (2020). Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat. Genet.* 52, 331-341.

2. **Korbel, J. O. & Campbell, P. J.** (2013). Criteria for inference of chromothripsis in cancer genomes. *Cell* 152, 1226-1236.

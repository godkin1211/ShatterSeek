# P-value Thresholds Configuration Guide

## Overview

Starting from v2.0.2, ShatterSeek allows you to configure the p-value thresholds used in chromothripsis classification. This provides flexibility to adjust the stringency based on your research goals.

## Default Thresholds (Tutorial Standard)

```R
classify_chromothripsis(chromoth_output)
# Equivalent to:
classify_chromothripsis(chromoth_output,
                        pval_fragment_joins_threshold = 0.05,
                        pval_clustering_threshold = 0.05)
```

**These are the standard significance levels from the original publication.**

## Parameters

### 1. `pval_fragment_joins_threshold` (default: 0.05)

**What it controls**: Fragment joins test threshold

**Test logic**: Fragment joins p-value must be **GREATER** than this threshold to pass
- Pass: p > threshold → SV types are randomly distributed (consistent with chromothripsis)
- Fail: p ≤ threshold → SV types are NOT random (suggests other mechanism)

**Example**:
```R
# Default (0.05): Accepts p-values like 0.06, 0.10, 0.50
# Stringent (0.01): Rejects p-values between 0.01-0.05 (e.g., 0.02, 0.03)
```

### 2. `pval_clustering_threshold` (default: 0.05)

**What it controls**: Clustering tests threshold (chr enrichment & exp distribution)

**Test logic**: Clustering p-values must be **LESS** than this threshold to pass
- Pass: p < threshold → Significant clustering/enrichment
- Fail: p ≥ threshold → Not significantly clustered

**Example**:
```R
# Default (0.05): Accepts p-values like 0.04, 0.01, 0.001
# Stringent (0.01): Rejects p-values between 0.01-0.05 (e.g., 0.02, 0.03)
```

## Use Cases

### Scenario 1: Discovery Phase (Default: 0.05)

**Goal**: Identify all potential chromothripsis events

```R
results <- classify_chromothripsis(chromoth_output)
```

**Characteristics**:
- Balanced sensitivity and specificity
- Standard statistical significance
- Recommended for initial screening
- May include borderline cases

### Scenario 2: Validation Phase (Stringent: 0.01)

**Goal**: Only report high-confidence events

```R
results <- classify_chromothripsis(chromoth_output,
                                   pval_fragment_joins_threshold = 0.01,
                                   pval_clustering_threshold = 0.01)
```

**Characteristics**:
- Higher specificity (fewer false positives)
- Lower sensitivity (may miss borderline cases)
- Suitable for clinical reporting
- Reduces risk of over-calling

### Scenario 3: Custom Thresholds

**Goal**: Different stringency for different tests

```R
# More stringent fragment joins, standard clustering
results <- classify_chromothripsis(chromoth_output,
                                   pval_fragment_joins_threshold = 0.01,
                                   pval_clustering_threshold = 0.05)
```

## Risk Analysis

### Changing to 0.01: What You Lose and Gain

#### Fragment Joins Threshold: 0.05 → 0.01

**Events Lost**:
- Chromosomes with fragment joins p-values between 0.01-0.05
- Example: p = 0.02 (weak randomness, but still possible chromothripsis)

**Impact**:
- ❌ Reduced sensitivity for "weakly random" events
- ✅ Increased confidence that detected events are truly random

**Risk Level**: **Low** for most real chromothripsis
- True chromothripsis usually has p >> 0.05 (often > 0.5)
- Borderline cases (0.01 < p < 0.05) may be ambiguous anyway

#### Clustering Threshold: 0.05 → 0.01

**Events Lost**:
- Chromosomes with clustering p-values between 0.01-0.05
- Example: p = 0.02 (weakly significant clustering)

**Impact**:
- ❌ May miss events with modest clustering
- ✅ Ensures highly significant spatial clustering

**Risk Level**: **Moderate**
- Some real chromothripsis may have moderate clustering (0.01 < p < 0.05)
- Consider using 0.05 for clustering if discovery is the goal

### Recommended Thresholds by Use Case

| Use Case | Fragment Joins | Clustering | Rationale |
|----------|---------------|------------|-----------|
| **Discovery** | 0.05 | 0.05 | Standard; balanced |
| **Validation** | 0.01 | 0.01 | High confidence; clinical |
| **Mixed** | 0.01 | 0.05 | Strict randomness, allow moderate clustering |
| **Conservative** | 0.01 | 0.001 | Maximum specificity |

## Testing Your Thresholds

### Check Impact on Your Data

```R
library(ShatterSeek)

# Load your data
SV_data <- ...
CN_data <- ...
result <- detect_chromothripsis(SV_data, CN_data, genome="hg19")

# Compare thresholds
class_default <- classify_chromothripsis(result)
class_stringent <- classify_chromothripsis(result,
                                          pval_fragment_joins_threshold = 0.01,
                                          pval_clustering_threshold = 0.01)

# Count changes
cat("Default (0.05):\n")
table(class_default$classification)

cat("\nStringent (0.01):\n")
table(class_stringent$classification)

# Identify borderline cases
library(dplyr)
comparison <- data.frame(
    chrom = class_default$chrom,
    default = class_default$classification,
    stringent = class_stringent$classification,
    pval_frag = class_default$pval_fragment_joins,
    pval_enrich = class_default$chr_breakpoint_enrichment
)

# Find cases that changed
changed <- comparison[comparison$default != comparison$stringent, ]
if (nrow(changed) > 0) {
    cat("\nCases that changed classification:\n")
    print(changed)
} else {
    cat("\nNo cases changed classification ✓\n")
}
```

### Diagnostic Analysis

```R
# Run diagnostic with custom thresholds
diagnostic <- diagnose_chromothripsis_classification(
    result,
    pval_fragment_joins_threshold = 0.01,
    pval_clustering_threshold = 0.01
)

# View detailed report
print_diagnostic_report(diagnostic)
```

The report will show:
- Thresholds used at the top
- Which criteria each chromosome passes/fails
- Specific p-values compared to your thresholds

## Real-World Examples

### Example 1: DO17373 (Sample from Original Paper)

```R
data(DO17373)
# ... prepare data ...
result <- detect_chromothripsis(SV_data, CN_data, genome="hg19")

# Test both thresholds
class_005 <- classify_chromothripsis(result)
class_001 <- classify_chromothripsis(result,
                                     pval_fragment_joins_threshold = 0.01,
                                     pval_clustering_threshold = 0.01)
```

**Result**: No difference
- All detected chromothripsis have very extreme p-values
- Fragment joins: p > 0.5 (far above both thresholds)
- Clustering: p < 1e-08 (far below both thresholds)

**Conclusion**: Well-defined chromothripsis is robust to threshold choice

### Example 2: Borderline Case (Hypothetical)

```
Chromosome 5:
  - Intrachromosomal SVs: 8 ✓
  - CN oscillations: 9 ✓
  - Fragment joins: p = 0.03 (borderline!)
  - Chr enrichment: p = 0.02 (borderline!)
```

**With threshold = 0.05**: **High confidence** (both pass)
**With threshold = 0.01**: **Not chromothripsis** (both fail)

**Interpretation**: This is an ambiguous case
- Weak randomness (p=0.03 not very high)
- Weak enrichment (p=0.02 not very low)
- May be mixed mechanism or early-stage event
- Consider additional evidence (e.g., breakpoint sequences)

## Best Practices

### 1. Report Your Thresholds

**Always document which thresholds you used:**

```R
# Good
cat("Classification performed with p-value thresholds:\n")
cat("  - Fragment joins: p > 0.01\n")
cat("  - Clustering: p < 0.01\n")
```

### 2. Analyze Borderline Cases

```R
# Identify borderline cases
summary <- result@chromSummary[result@chromSummary$clusterSize > 0, ]

borderline_frag <- summary[!is.na(summary$pval_fragment_joins) &
                           summary$pval_fragment_joins > 0.01 &
                           summary$pval_fragment_joins <= 0.05, ]

if (nrow(borderline_frag) > 0) {
    cat("⚠ Borderline fragment joins cases:\n")
    print(borderline_frag[, c("chrom", "pval_fragment_joins")])
}
```

### 3. Consider Cohort Analysis

For large cohorts:
- Use default (0.05) for discovery
- Apply stringent (0.01) for subset validation
- Report both analyses for transparency

### 4. Biological Context

Remember:
- Thresholds are statistical conventions, not biological truth
- A p-value of 0.03 vs 0.07 may not be biologically meaningful
- Consider other evidence (CN patterns, SV types, clinical context)

## References

1. **Fisher, R. A.** (1925). Statistical Methods for Research Workers.
   - Origin of p < 0.05 convention

2. **Cortes-Ciriano, I. et al.** (2020). Nat. Genet. 52, 331-341.
   - Original ShatterSeek publication using 0.05 threshold

3. **Benjamin, D. J. et al.** (2018). Nat. Hum. Behav. 2, 6-10.
   - Arguments for p < 0.005 in some contexts

## Support

For questions about threshold selection:
- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues
- Email: n28111021@gs.ncku.edu.tw

# Fragment Joins Test Bug Analysis

## Summary

**Bug discovered**: The bash script used **inverted logic** for the fragment joins test, causing false positive chromothripsis detections.

## The Bug

### Wrong Logic (Original Script)
```bash
if (numIntra>=6 && $16>=7 && $12<0.05 && ($13<0.05 || $14<0.05))
#                              ^^^^^^^^ WRONG!
```

### Correct Logic
```bash
if (numIntra>=6 && $16>=7 && $12>0.05 && ($13<0.05 || $14<0.05))
#                              ^^^^^^^^ CORRECT!
```

## Why This Matters

The fragment joins test (chi-squared test) evaluates if SV types (DEL, DUP, INV) are **randomly distributed**:

| P-value | Interpretation | Chromothripsis? |
|---------|---------------|----------------|
| **p > 0.05** | SV types are **random** (no bias) | **PASS** ✓ - Consistent with chromothripsis |
| **p ≤ 0.05** | SV types are **non-random** (biased) | **FAIL** ✗ - Suggests other mechanism |

**Chromothripsis** is characterized by **random reassembly** of chromosome fragments, so we expect:
- Approximately equal proportions of DEL, DUP, and INV (~33% each)
- High p-value (>>0.05) indicating random distribution

**The bug**: Script accepted events with **p ≤ 0.05**, meaning it detected events with **biased** SV type distributions as chromothripsis!

## Impact on chromothripsis.tsv Analysis

### With WRONG Logic (p < 0.05)

**Detected**: 4 "low confidence" events

| Chr | IntraSVs | CN_Osc | Fragment Joins p | Chr Enrich p | Classification |
|-----|----------|--------|-----------------|--------------|----------------|
| 6   | 52       | 5      | 0.00001         | 1.09e-63     | Low            |
| 15  | 23       | 4      | 0.00026         | 2.60e-30     | Low            |
| 17  | 23       | 6      | 0.00135         | 1.53e-52     | Low            |
| 20  | 24       | 5      | 0.02999         | 4.55e-24     | Low            |

**Problem**: All have fragment joins p < 0.05, meaning SV types are **non-random**!

**Example - Chr 6**:
- Fragment joins p = 0.00001 (highly non-random)
- SV distribution: DEL=39 (75%), h2hINV=0, t2tINV=0, DUP=13 (25%)
- **Interpretation**: Strongly biased toward DEL; NOT chromothripsis-like
- **Likely mechanism**: BFB (Bridge-Fusion-Breakage) or other deletion-based process

### With CORRECT Logic (p > 0.05)

**Detected**: Only 1 chromosome (Chr 1)

| Chr | IntraSVs | CN_Osc | Fragment Joins p | Chr Enrich p | Exp Chr p | Classification |
|-----|----------|--------|-----------------|--------------|-----------|----------------|
| 1   | 21       | 7      | 0.910           | 5.89e-76     | 2.22e-16  | High confidence |

**Why Chr 1 passes**:
- Fragment joins p = 0.910 (very high, highly random)
- SV distribution: DEL=7 (33%), h2hINV=0, t2tINV=0, DUP=14 (67%)
- **Interpretation**: Reasonably balanced distribution
- **Note**: Still no inversions (h2hINV=0, t2tINV=0), which is somewhat unusual

**Why the 4 previous events now FAIL**:
- All have fragment joins p < 0.05
- All show biased SV type distributions
- These are NOT chromothripsis by current standards

## SV Type Distribution Analysis

### Chromothripsis (Expected)
- Random reassembly → Random SV types
- Fragment joins p >> 0.05 (often > 0.5)
- Approximately equal DEL, DUP, INV proportions

### Chr 6 (Failed with correct logic)
```
DEL: 75% |████████████████████████████████████████|
DUP: 25% |█████████████|
INV:  0% |
Fragment joins p = 0.00001 (FAIL)
```
**Mechanism**: BFB or deletion-dominant process

### Chr 15 (Failed with correct logic)
```
DEL: 39% |████████████████████|
DUP: 61% |██████████████████████████████|
INV:  0% |
Fragment joins p = 0.00026 (FAIL)
```
**Mechanism**: Duplication-dominant, possibly chromosynthesis

### Chr 17 (Failed with correct logic)
```
DEL: 52% |██████████████████████████|
DUP: 48% |████████████████████████|
INV:  0% |
Fragment joins p = 0.00135 (FAIL)
```
**Mechanism**: Despite balanced DEL/DUP, p-value indicates non-random pattern

### Chr 1 (PASSED with correct logic)
```
DEL: 33% |█████████████████|
DUP: 67% |█████████████████████████████████|
INV:  0% |
Fragment joins p = 0.910 (PASS)
```
**Mechanism**: High p-value indicates random joining despite DUP bias

## Biological Interpretation

### The 4 False Positives (from wrong logic)

These events have **non-random SV type distributions** (p < 0.05), suggesting:

1. **Chr 6 (75% DEL)**: Likely **BFB** (Bridge-Fusion-Breakage)
   - Characterized by repeated deletion cycles
   - Creates clustered deletions on chromosome arms

2. **Chr 15 (61% DUP)**: Possibly **chromosynthesis**
   - Serial template-switching during replication
   - Generates tandem duplications

3. **Chr 17 & Chr 20**: Mixed mechanisms or **chromoplexy**
   - Even though DEL/DUP proportions appear balanced
   - Statistical test detected non-random patterns

### True Chromothripsis (Chr 1)

- High fragment joins p-value (0.910) indicates random reassembly
- Passes all other criteria (≥7 CN oscillations, clustering tests)
- Consistent with catastrophic chromosome shattering

## Recommendation

**For chromothripsis.tsv data**:

1. **Do NOT use the 4 events** detected with wrong logic
   - They are not chromothripsis by current standards
   - They represent other chromoanagenesis mechanisms

2. **Only Chr 1** meets chromothripsis criteria
   - High confidence classification
   - Random fragment joining pattern

3. **Consider comprehensive analysis**:
   ```R
   # Detect all chromoanagenesis types
   results <- detect_chromoanagenesis(SV_data, CN_data, genome="hg19")

   # This will identify:
   # - Chromothripsis (random joining, Chr 1)
   # - Chromoplexy (translocation chains)
   # - Chromosynthesis (replication-based, Chr 15?)
   ```

## Testing the Fix

To verify the corrected logic matches our R implementation:

```R
library(ShatterSeek)

# Load chromothripsis.tsv data
chrom_summary <- read.table("chromothripsis.tsv", header=TRUE, sep="\t")

# Apply correct classification logic
classify_chromothripsis(chrom_output,
                        pval_fragment_joins_threshold = 0.05,
                        pval_clustering_threshold = 0.05)

# Expected result: Only Chr 1 detected as high confidence
```

## Key Takeaway

**Fragment joins test logic**:
- **Correct**: `p > 0.05` (random joining → chromothripsis)
- **Wrong**: `p < 0.05` (non-random joining → other mechanisms)

This single character difference (`>` vs `<`) caused the detection of 4 false positive chromothripsis events that are actually other types of chromoanagenesis.

# Understanding v2.0.1 Classification Changes

## Quick Answer: Why Are My Samples Not Showing Chromothripsis Anymore?

The most likely reason is the **fragment joins test**, which was previously calculated but **not enforced**. This test is now mandatory in v2.0.1 to comply with the original publication criteria.

## What Changed?

### 1. Fragment Joins Test Now Enforced (CRITICAL)

**Before v2.0.1:**
- The test was calculated (`pval_fragment_joins`)
- But **NOT checked** during classification
- Samples with non-random SV distributions could still be classified as chromothripsis

**After v2.0.1:**
- Fragment joins test p-value > 0.05 is **REQUIRED**
- Samples with biased SV type distributions are now correctly excluded

**Impact:** This is the most common reason for missing detections.

### 2. What Does the Fragment Joins Test Check?

The test asks: **Are SV types (DEL, DUP, INV) randomly distributed?**

**Expected for chromothripsis:**
- Random NHEJ repair → approximately equal proportions
- DEL ~33%, DUP ~33%, INV ~33%
- Chi-squared test p-value > 0.05

**NOT chromothripsis if:**
- Biased toward specific SV types
- Example: 70% DUP, 20% DEL, 10% INV
- Chi-squared test p-value ≤ 0.05
- Suggests replication-based mechanism (chromosynthesis, FoSTeS/MMBIR)

### 3. Other Changes

- **CN oscillation upper bound** for low confidence (now 4-6, not ≥4)
- **Classification names** changed:
  - "Likely" → "High confidence"
  - "Possible" → "Low confidence"
  - "Unlikely" → "Not chromothripsis"

## How to Diagnose Your Samples

### Step 1: Install Updated Package

```bash
cd /Users/godkin/Projects/ShatterSeek
R CMD INSTALL .
```

### Step 2: Run Quick Diagnostic

```bash
# Use the diagnostic script
Rscript diagnose_my_sample.R
```

**Edit the script first** to load your data (see "USER INPUT" section).

### Step 3: Interpret Results

The diagnostic report will show:

1. **Overall classification** (High/Low/Not chromothripsis)
2. **Detailed criteria breakdown** for each chromosome
3. **Fragment joins p-values** with interpretation
4. **Blocking reasons** for "Not chromothripsis" classifications

### Step 4: Understand Common Patterns

#### Pattern 1: Fragment Joins Test Failing

```
Chromosome 5: Not chromothripsis
  - Fragment joins: p=0.03 ✗ (≤0.05 = non-random)
  - SV types: DEL=2, h2hINV=1, t2tINV=0, DUP=15 (DUP=83%, DEL=11%, INV=6%)
  → HIGH DUP proportion suggests replication-based mechanism (chromosynthesis?)
```

**What this means:**
- Your sample has too many duplications
- This is NOT chromothripsis (which has random joining)
- Consider running `detect_chromosynthesis()` instead

**Action:**
```R
# Check for other mechanisms
results <- detect_chromoanagenesis(SV_data, CN_data)
```

#### Pattern 2: Insufficient CN Oscillations

```
Chromosome 8: Not chromothripsis
  - CN oscillations: 3 ✗ (<4)
  - Intrachromosomal SVs: 8 ✓ (≥6)
  - Fragment joins: p=0.45 ✓ (>0.05)
```

**What this means:**
- Not enough copy number changes
- May be CNV caller sensitivity issue

**Action:**
- Check CNV data quality
- Verify adjacent segments are properly merged
- Consider more sensitive CNV caller

#### Pattern 3: Too Few SVs

```
Chromosome 12: Not chromothripsis
  - Intrachromosomal SVs: 4 ✗ (<6)
```

**What this means:**
- Below published threshold
- May be SV caller sensitivity or filtering

**Action:**
- Check SV size filtering
- Review SV caller sensitivity
- Consider using multiple callers

## Is This a Bug or Expected Behavior?

**This is EXPECTED behavior** and represents a bug fix, not a new bug.

### Why Previous Versions Were Wrong

The original ShatterSeek publication (Cortes-Ciriano et al. 2020) clearly states:

> "Chromothripsis involves random joining of DNA fragments through NHEJ,
> resulting in random SV type distribution (p > 0.05 by chi-squared test)"

Previous versions calculated this test but didn't enforce it, leading to **false positive chromothripsis classifications** for events that were actually:
- Chromosynthesis (serial replication)
- FoSTeS/MMBIR (template switching)
- BFB (breakage-fusion-bridge)
- Other non-random mechanisms

### Why Your Samples May Be Different Mechanisms

If your samples consistently fail the fragment joins test, they may genuinely be:

1. **Chromosynthesis** (replication-based)
   - Serial template switching
   - Tandem duplications predominate
   - Gradual CN increases
   - Detection: `detect_chromosynthesis()`

2. **Chromoplexy** (translocation chains)
   - Multi-chromosome chains
   - Balanced translocations
   - Stable CN
   - Detection: `detect_chromoplexy()`

3. **BFB** (breakage-fusion-bridge)
   - Inverted duplications
   - CN amplifications
   - Classic pattern

## What Should You Do?

### Option 1: Accept the New Classification (Recommended)

If your samples fail the fragment joins test:

1. **They are likely NOT chromothripsis** (by strict definition)
2. Run comprehensive analysis:
   ```R
   results <- detect_chromoanagenesis(SV_data, CN_data)
   mixed <- classify_mixed_mechanisms(results)
   ```
3. Check for chromosynthesis or chromoplexy
4. Report the correct mechanism

### Option 2: Investigate Data Quality

If you believe there's a data issue:

1. **Check SV type encoding:**
   ```R
   table(SV_data@SVtype)
   table(paste(SV_data@strand1, SV_data@strand2))
   ```

2. **Verify caller output:**
   - Are SV types correctly assigned?
   - Is there caller-specific bias?

3. **Review filtering:**
   - Size filters
   - Quality filters
   - Mappability filters

### Option 3: Compare with Published Data

Run the diagnostic on the example data:

```bash
cd tests/validation
Rscript diagnose_example.R
```

The DO17373 example should show:
- 4 high confidence chromothripsis (chr 2, 3, 21, X)
- All pass fragment joins test (p > 0.05)

If the example works but your data doesn't, the issue is with your data/caller, not ShatterSeek.

## Real-World Example

Let's say your sample previously classified as "Likely chromothripsis" on chromosome 5:

**Old v2.0.0 classification:**
```
Chr 5: Likely chromothripsis
  - 8 intrachromosomal SVs ✓
  - 5 CN oscillations ✓
  - Clustering p=0.01 ✓
  - Fragment joins p=0.03 (IGNORED)
```

**New v2.0.1 classification:**
```
Chr 5: Not chromothripsis
  - 8 intrachromosomal SVs ✓
  - 5 CN oscillations ✓
  - Clustering p=0.01 ✓
  - Fragment joins p=0.03 ✗ (FAIL - now enforced!)

SV types: DEL=1, DUP=6, INV=1 (DUP=75%)
→ High DUP proportion → Likely chromosynthesis
```

**The truth:** This was never chromothripsis. It was always a replication-based event, but v2.0.0 incorrectly classified it.

## Summary Decision Tree

```
My sample no longer shows chromothripsis
    |
    ├── Run diagnostic script (diagnose_my_sample.R)
    |
    ├── Fragment joins p ≤ 0.05?
    |   ├── YES → Check SV type distribution
    |   |   ├── Biased toward DUP (>50%) → Likely chromosynthesis
    |   |   ├── Biased toward INV (>50%) → Likely BFB
    |   |   ├── Random but p≤0.05 → Small sample size, borderline
    |   |   └── Verify SV type encoding correctness
    |   |
    |   ├── CN oscillations < 4?
    |   |   └── Check CNV caller sensitivity / segment merging
    |   |
    |   ├── Too few SVs (<6)?
    |   |   └── Check SV caller sensitivity / filtering
    |   |
    |   └── Clustering failed?
    |       └── SVs not spatially clustered → Not chromothripsis
    |
    └── All tests pass but still "Not chromothripsis"?
        └── Review all criteria together - may need ALL to pass
```

## Files to Help You

1. **diagnose_my_sample.R** - Quick diagnostic script (edit and run)
2. **DIAGNOSTIC_GUIDE.md** - Comprehensive troubleshooting guide
3. **CHROMOTHRIPSIS_CLASSIFICATION.md** - Detailed classification criteria
4. **CHANGES_v2.0.1.md** - Complete changelog

## Getting Help

If you're still confused or believe there's an error:

1. **Run diagnostic and save output:**
   ```R
   diagnostic <- diagnose_chromothripsis_classification(chromoth_result)
   write.csv(diagnostic, "my_diagnostic.csv")
   ```

2. **Check your SV type distribution:**
   ```R
   chromoth_result@chromSummary[, c("chrom", "number_DEL", "number_DUP",
                                     "number_h2hINV", "number_t2tINV",
                                     "pval_fragment_joins")]
   ```

3. **Report issue with details:**
   - GitHub: https://github.com/godkin1211/ShatterSeek/issues
   - Email: n28111021@gs.ncku.edu.tw
   - Include: SV caller, CNV caller, diagnostic output

## References

1. **Cortes-Ciriano, I. et al.** (2020). Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat. Genet.* 52, 331-341.
   - **Key quote (p. 333):** "We applied a chi-squared test to assess whether the proportions of deletion, tandem duplication, and inversion breakpoints were consistent with random DNA fragment joining."

2. **Korbel, J. O. & Campbell, P. J.** (2013). Criteria for inference of chromothripsis in cancer genomes. *Cell* 152, 1226-1236.
   - **Key quote (Box 1):** "Random joining of fragments... implies that deletion-type, tandem duplication-type, and inversion-type rearrangements occur in similar proportions."

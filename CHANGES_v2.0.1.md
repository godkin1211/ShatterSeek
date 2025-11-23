# ShatterSeek Extended Edition v2.0.1 - Classification Update

## Summary

This update aligns ShatterSeek's chromothripsis classification system with the original publication criteria from Cortes-Ciriano et al. (Nature Genetics 2020) and Korbel & Campbell (Cell 2013).

## Major Changes

### 1. Classification Categories Renamed

**Old System (v2.0.0)**:
- Likely chromothripsis
- Possible chromothripsis
- Unlikely

**New System (v2.0.1)**:
- High confidence
- Low confidence
- Not chromothripsis

### 2. Fragment Joins Test Now Enforced

**Critical Bug Fix**: The fragment joins test (chi-squared test for random SV type distribution) was previously calculated but **not enforced** in the classification logic.

**Impact**:
- Previous versions could incorrectly classify events with non-random fragment joining as chromothripsis
- Events biased toward specific SV types (e.g., replication-based mechanisms) may have been misclassified

**Fix**: Fragment joins test (p > 0.05) is now **mandatory** for both high and low confidence classification.

### 3. Classification Criteria Fully Implemented

#### High Confidence Chromothripsis (Type 1)
```
✓ ≥6 interleaved intrachromosomal SVs
✓ ≥7 contiguous CN segments oscillating between 2 states
✓ Fragment joins test passed (p > 0.05)
✓ Clustering test passed (chr enrichment p < 0.05 OR exp. distribution p < 0.05)
```

#### High Confidence Chromothripsis (Type 2)
```
✓ ≥3 interleaved intrachromosomal SVs
✓ ≥4 interchromosomal SVs (translocations)
✓ ≥7 contiguous CN segments oscillating between 2 states
✓ Fragment joins test passed (p > 0.05)
```

#### Low Confidence Chromothripsis
```
✓ ≥6 interleaved intrachromosomal SVs
✓ 4-6 contiguous CN segments oscillating between 2 states (NOTE: upper bound!)
✓ Fragment joins test passed (p > 0.05)
✓ Clustering test passed (chr enrichment p < 0.05 OR exp. distribution p < 0.05)
```

### 4. CN Oscillation Upper Bound for Low Confidence

**Old System**: Low confidence could have ≥4 CN oscillations (no upper limit)

**New System**: Low confidence requires 4-6 CN oscillations (upper bound enforced)

**Rationale**: Events with ≥7 oscillations should be classified as high confidence, not low confidence.

## Files Modified

### Core Classification Logic
- ✅ `R/scoring_functions.R`
  - `classify_chromothripsis()`: Complete rewrite with tutorial-compliant criteria
  - `summarize_chromothripsis()`: Updated output format and terminology
  - Added detailed documentation with references

### Integration and Output
- ✅ `R/chromoanagenesis_integrated.R`
  - Updated `detect_chromoanagenesis()` to use new classification levels
  - Updated `create_integrated_summary()` for high/low confidence counts
  - Updated `print.chromoanagenesis()` output format

### Visualization Functions
- ✅ `R/visualization_utils.R` (2 locations)
- ✅ `R/ggnome_style_visualization.R`
- ✅ `R/genome_dashboard_visualization.R`
- ✅ `R/circos_visualization.R`
- ✅ `R/integrated_classifier.R`

All updated to recognize "High confidence" and "Low confidence" instead of "Likely" and "Possible".

### Documentation
- ✅ `README.md`
  - Added classification criteria summary
  - Updated examples to show new classification workflow
  - Added reference to detailed classification guide

- ✅ `CHROMOTHRIPSIS_CLASSIFICATION.md` (NEW)
  - Comprehensive guide to classification criteria
  - Detailed explanation of statistical tests
  - CN oscillation counting methodology
  - Usage examples with expected output
  - Complete references

- ✅ `CLAUDE.md`
  - Updated project documentation references

### Testing
- ✅ `tests/validation/test_new_classification.R` (NEW)
  - Validates new classification system
  - Tests fragment joins enforcement
  - Verifies CN oscillation requirements
  - Checks classification category validity

## Migration Guide

### For Users

**Code changes required**: Minimal

**Old code**:
```R
result <- classify_chromothripsis(chromoth_output,
                                 min_cluster_size = 6,
                                 min_cn_oscillations = 4,
                                 max_pval_clustering = 0.05)
```

**New code**:
```R
# Parameters removed - now uses tutorial criteria
result <- classify_chromothripsis(chromoth_output)
```

**Output changes**:
```R
# Old
result$classification  # "Likely chromothripsis", "Possible chromothripsis", "Unlikely"

# New
result$classification  # "High confidence", "Low confidence", "Not chromothripsis"
```

### For Developers

**Breaking Changes**:
- `classify_chromothripsis()` function signature changed (parameters removed)
- Classification levels renamed throughout codebase
- Fragment joins test now mandatory

**Backward Compatibility**:
- `confidence_score` field still present for ranking
- `chromSummary` slot structure unchanged
- All statistical tests still calculated and available

## Testing

Run the new validation script:

```bash
cd tests/validation
Rscript test_new_classification.R
```

Expected output should show:
- High confidence events with ≥7 CN oscillations
- Low confidence events with 4-6 CN oscillations
- All chromothripsis events passing fragment joins test (p > 0.05)

## References

1. **Cortes-Ciriano, I. et al.** Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. *Nat. Genet.* 52, 331-341 (2020).

2. **Korbel, J. O. & Campbell, P. J.** Criteria for inference of chromothripsis in cancer genomes. *Cell* 152, 1226-1236 (2013).

3. **Stephens, P. J. et al.** Massive genomic rearrangement acquired in a single catastrophic event during cancer development. *Cell* 144, 27-40 (2011).

## Version History

- **v2.0.1** (Current): Tutorial-compliant classification with fragment joins test enforcement
- **v2.0.0**: Extended Edition with chromoplexy, chromosynthesis, and visualization enhancements
- **v1.x**: Original chromothripsis-only framework

## Contact

For questions about the classification update:
- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues
- Email: n28111021@gs.ncku.edu.tw

## Acknowledgments

Classification criteria based on the original ShatterSeek framework developed by:
- Isidro Cortes-Ciriano (University of Cambridge)
- Peter J. Park (Harvard Medical School)

Tutorial implementation verified against:
- Tutorial.pdf from original ShatterSeek distribution
- Nature Genetics 2020 publication supplementary methods

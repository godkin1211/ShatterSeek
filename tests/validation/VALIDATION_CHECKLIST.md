# Phase 1 Validation Checklist

Use this checklist to track your validation progress.

## Pre-Validation Setup

- [ ] R (>= 3.5.0) installed
- [ ] ShatterSeek Extended Edition installed
  ```bash
  cd ShatterSeek
  R -e "devtools::install()"
  ```
- [ ] Navigate to validation directory
  ```bash
  cd tests/validation
  ```

## Phase 1 Validation Tests

### Test 1: Built-in Dataset (DO17373)

**Status**: [ ] PASS / [ ] FAIL / [ ] NOT RUN

**Commands**:
```bash
# Run in validation directory
Rscript phase1_validation.R
```

**What to verify**:
- [ ] SV and CNV data loads successfully
- [ ] Chromothripsis detected on chromosome(s)
- [ ] Classification produces results
- [ ] No errors in console output

**Expected results**:
```
✓ Detected chromothripsis on X chromosomes
✓ Classification complete: X chromosomes classified
✓ Original and extended results are consistent
✓ Built-in data test PASSED
```

---

### Test 2: Synthetic Data Generation & Detection

**Status**: [ ] PASS / [ ] FAIL / [ ] NOT RUN

**What to verify**:
- [ ] 10 synthetic chromothripsis events generated
- [ ] At least 7 out of 10 detected (≥70%)
- [ ] No errors during generation or detection

**Expected results**:
```
Test case 1/10: DETECTED
Test case 2/10: DETECTED
...
Sensitivity on synthetic data: XX.X% (≥7/10)
✓ Synthetic data test PASSED
```

---

### Test 3: BND/TRA ALT Field Parsing

**Status**: [ ] PASS / [ ] FAIL / [ ] NOT RUN

**What to verify**:
- [ ] Test VCF created successfully
- [ ] VCF parsing completes without errors
- [ ] TRA records identified
- [ ] All TRA records are interchromosomal
- [ ] All bracket formats parse correctly

**Expected results**:
```
✓ Parsed X SV records
✓ Found X TRA records
✓ All TRA records are interchromosomal
'[chr13:49291490[TT' -> chr=13, pos=49291490 ✓
'T]chr1:245088964]' -> chr=1, pos=245088964 ✓
...
✓ BND parsing test PASSED
```

---

## Overall Validation Summary

**Status**: [ ] ALL PASSED / [ ] SOME FAILED / [ ] NOT RUN

**Final metrics**:
- Total tests: 3
- Passed: ___
- Failed: ___
- Success rate: ___%

**Required for success**: 100% (3/3 tests passed)

---

## Output Files Checklist

After running validation, verify these files exist:

- [ ] `phase1_validation_results.rds` - Detailed test results
- [ ] `validation_report.html` - HTML report
- [ ] `validation_summary.txt` - Text summary

**To view HTML report**:
```bash
open validation_report.html  # Mac/Linux
# Or open manually in web browser
```

---

## Troubleshooting Log

Use this space to note any issues encountered:

### Issue 1:
**Problem**:
**Solution**:
**Status**: [ ] RESOLVED / [ ] PENDING

### Issue 2:
**Problem**:
**Solution**:
**Status**: [ ] RESOLVED / [ ] PENDING

---

## Next Steps

Once all Phase 1 tests pass:

- [ ] Review HTML validation report
- [ ] Test with your own sample data
- [ ] Generate visualizations (circos, dashboard, regional plots)
- [ ] Plan Phase 2 validation (PCAWG datasets)
- [ ] Consider publishing validation results

---

## Sign-off

**Validation completed by**: ___________________

**Date**: ___________________

**R version**: ___________________

**ShatterSeek version**: 2.0.0

**Overall result**: [ ] PASS / [ ] FAIL

**Notes**:
_____________________________________________________________
_____________________________________________________________
_____________________________________________________________

---

## Quick Commands Reference

```bash
# Run everything
./run_all_validation.sh

# Individual tests
Rscript quick_check.R              # Quick check (~30s)
Rscript phase1_validation.R        # Full validation (~2-3min)
Rscript generate_report.R          # Generate reports
Rscript example_usage.R            # See examples

# View results
cat validation_summary.txt         # Text summary
open validation_report.html        # HTML report
```

## Support

If you encounter issues:
- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues
- Email: n28111021@gs.ncku.edu.tw
- Documentation: `RUNNING_VALIDATION.md`

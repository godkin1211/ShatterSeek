# Running Phase 1 Validation - Quick Start Guide

## Prerequisites

You need R (>= 3.5.0) installed on your system with the following packages:
- ShatterSeek Extended Edition (this package)
- All its dependencies (will be installed automatically)

## Step-by-Step Instructions

### 1. Clone and Install ShatterSeek

```bash
# Clone the repository
git clone https://github.com/godkin1211/ShatterSeek.git
cd ShatterSeek

# Install the package (choose one method)

# Method A: Using devtools (recommended)
R -e "devtools::install()"

# Method B: Using R CMD INSTALL
R CMD INSTALL .
```

### 2. Navigate to Validation Directory

```bash
cd tests/validation
```

### 3. Run Validation Tests

#### Option A: Run All Tests (Recommended)

```bash
# Run the complete validation suite
./run_all_validation.sh
```

This will:
1. Run quick check (~30 seconds)
2. Run comprehensive validation (~2-3 minutes)
3. Generate HTML and text reports

**Expected output:**
```
==========================================
ShatterSeek Extended Edition
Phase 1 Validation Suite
==========================================

Step 1/3: Running quick check...
----------------------------------------
=== Quick Validation Check ===

1. Loading built-in dataset... ✓
2. Testing chromothripsis detection... ✓
3. Testing comprehensive detection... ✓
4. Testing BND parsing... ✓

✓ Quick check PASSED - All core functions working
✓ Quick check passed

Step 2/3: Running comprehensive validation...
----------------------------------------
===========================================
ShatterSeek Extended Edition - Phase 1 Validation
===========================================

[Test 1] Built-in Dataset (DO17373)
------------------------------------------------------------
Loading built-in dataset DO17373...
  SV records: XXX
  CNV segments: XXX
  ...
✓ PASSED

[Test 2] Synthetic Data Generation & Detection
------------------------------------------------------------
Generating synthetic chromothripsis data...
  Test case 1/10: DETECTED
  Test case 2/10: DETECTED
  ...
✓ PASSED

[Test 3] BND/TRA ALT Field Parsing
------------------------------------------------------------
Testing BND ALT field parsing...
  ...
✓ PASSED

===========================================
VALIDATION SUMMARY
===========================================
Total tests: 3
Passed: 3
Failed: 0
Success rate: 100.0%
===========================================

✓ ALL TESTS PASSED

Step 3/3: Generating validation report...
----------------------------------------
✓ Report generated

==========================================
Validation Complete
==========================================

Output files:
  - phase1_validation_results.rds (detailed results)
  - validation_report.html (HTML report)
  - validation_summary.txt (text summary)

✓ All Phase 1 validation tests passed!
```

#### Option B: Run Individual Tests

```bash
# Quick sanity check (30 seconds)
Rscript quick_check.R

# Full validation (2-3 minutes)
Rscript phase1_validation.R

# Generate reports
Rscript generate_report.R

# See usage examples
Rscript example_usage.R
```

#### Option C: Run from R Console

```R
# Start R in the validation directory
setwd("tests/validation")

# Run quick check
source("quick_check.R")

# Run full validation
source("phase1_validation.R")

# Generate report
source("generate_report.R")
```

## Output Files

After running validation, you will get:

1. **phase1_validation_results.rds** - Detailed R object with all test results
2. **validation_report.html** - Beautiful HTML report (open in browser)
3. **validation_summary.txt** - Plain text summary

### Viewing the HTML Report

```bash
# Linux/Mac
open validation_report.html

# Or just open it manually in any web browser
```

## What Each Test Validates

### Test 1: Built-in Dataset (DO17373)
- ✅ Loads example data from original ShatterSeek paper
- ✅ Detects chromothripsis (original implementation)
- ✅ Detects all chromoanagenesis types (extended edition)
- ✅ Verifies classification works
- ✅ Confirms consistency between original and extended

### Test 2: Synthetic Data
- ✅ Generates 10 synthetic chromothripsis events
- ✅ Tests detection sensitivity (target: ≥70%)
- ✅ Validates different SV patterns
- ✅ Ensures robustness across parameters

### Test 3: BND Parsing
- ✅ Tests VCF breakend record parsing
- ✅ Validates bracket notation: `[chr:pos[`, `]chr:pos]`
- ✅ Verifies mate coordinate extraction
- ✅ Confirms translocation detection

## Success Criteria

✅ **All tests pass**: Success rate = 100%
- Built-in dataset: Chromothripsis detected
- Synthetic data: ≥70% detection rate (≥7 out of 10)
- BND parsing: All formats correctly parsed

⚠️ **Some tests fail**: Success rate < 100%
- Review error messages in console output
- Check `phase1_validation_results.rds` for details
- Report issues at: https://github.com/godkin1211/ShatterSeek/issues

## Troubleshooting

### Issue: "package 'ShatterSeek' is not installed"

**Solution:**
```R
# Install from local directory
devtools::install()

# Or install from GitHub
devtools::install_github("godkin1211/ShatterSeek")
```

### Issue: "DO17373 dataset not found"

**Solution:** The dataset should be included in the package. Try:
```R
data(package = "ShatterSeek")  # List all datasets
data(DO17373)                  # Load the dataset
```

### Issue: VCF reading errors

**Solution:** Install VCF parsing dependencies:
```R
# Option 1: VariantAnnotation (recommended)
BiocManager::install("VariantAnnotation")

# Option 2: vcfR (lightweight)
install.packages("vcfR")
```

### Issue: Visualization errors

**Solution:** Install optional visualization packages:
```R
# For circos plots
BiocManager::install("circlize")

# For breakpoint analysis
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
```

## Next Steps After Validation

Once Phase 1 passes:

1. **Review the HTML report**
   - Open `validation_report.html` in browser
   - Check all tests show "PASSED"

2. **Test with your own data**
   ```R
   # Load your VCF files
   sv_data <- read_sv_vcf("your_sample.sv.vcf.gz", caller = "dragen")
   cnv_data <- read_cnv_vcf("your_sample.cnv.vcf.gz", caller = "cnvkit")

   # Run detection
   results <- detect_chromoanagenesis(sv_data, cnv_data)

   # Visualize
   plot_chromoanagenesis_circos(results, sv_data, cnv_data, "YourSample")
   ```

3. **Proceed to Phase 2 validation**
   - PCAWG dataset comparison
   - Statistical properties validation
   - Cross-tool comparison

## Getting Help

- **GitHub Issues**: https://github.com/godkin1211/ShatterSeek/issues
- **Email**: n28111021@gs.ncku.edu.tw
- **Documentation**: See `README.md` in validation directory

## Quick Reference

```bash
# Complete validation workflow
cd ShatterSeek/tests/validation
./run_all_validation.sh          # Run everything
open validation_report.html      # View results
```

That's it! Happy validating! 🚀

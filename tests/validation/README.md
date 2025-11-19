# ShatterSeek Extended Edition - Validation Framework

This directory contains validation tests for ShatterSeek Extended Edition v2.0+.

## Phase 1: Core Functionality Validation

### Tests Included

1. **Built-in Dataset Test (DO17373)**
   - Tests chromothripsis detection on the example dataset from the original paper
   - Verifies both original and extended edition detection
   - Checks consistency between implementations
   - Validates classification functionality

2. **Synthetic Data Test**
   - Generates synthetic chromothripsis events with known characteristics
   - Tests detection sensitivity (target: ≥70%)
   - Validates handling of different SV types and CN patterns
   - Runs 10 independent test cases

3. **BND/TRA Parsing Test**
   - Tests VCF breakend (BND) ALT field parsing
   - Validates all bracket notation formats: `[chr:pos[`, `]chr:pos]`, etc.
   - Verifies interchromosomal event detection
   - Ensures TRA records are correctly identified

## Running the Tests

### Quick Start

```bash
# From ShatterSeek root directory
cd tests/validation
Rscript phase1_validation.R
```

### Expected Output

```
===========================================
ShatterSeek Extended Edition - Phase 1 Validation
===========================================

Starting Phase 1 Validation Tests...
====================================

[Test 1] Built-in Dataset (DO17373)
------------------------------------------------------------
Loading built-in dataset DO17373...
  SV records: XXX
  CNV segments: XXX
  Testing chromothripsis detection...
  ✓ Detected chromothripsis on X chromosomes
  ...
✓ PASSED

[Test 2] Synthetic Data Generation & Detection
------------------------------------------------------------
Generating synthetic chromothripsis data...
  Test case 1/10: DETECTED
  Test case 2/10: DETECTED
  ...
  Sensitivity on synthetic data: XX.X% (X/10)
✓ PASSED

[Test 3] BND/TRA ALT Field Parsing
------------------------------------------------------------
Testing BND ALT field parsing...
  Created test VCF: /tmp/RtmpXXXXXX/fileXXXXXX.vcf
  Parsing VCF...
  ✓ Parsed X SV records
  ✓ Found X TRA records
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
```

## Interpreting Results

### Success Criteria

- **Built-in Dataset**: Must detect chromothripsis on known positive sample
- **Synthetic Data**: Must achieve ≥70% detection rate
- **BND Parsing**: Must correctly parse all BND bracket formats

### Output Files

- `phase1_validation_results.rds`: Detailed test results (R object)
- Console output: Human-readable summary

## Troubleshooting

### Common Issues

**Error: "package 'ShatterSeek' is not installed"**
```R
# Install ShatterSeek first
devtools::install()  # from ShatterSeek root directory
```

**Error: "DO17373 dataset not found"**
```R
# The dataset should be included in the package
# Check installation and try:
data(package = "ShatterSeek")
```

**BND parsing test fails**
- Ensure VCF reading dependencies are installed:
```R
# Option 1: VariantAnnotation
BiocManager::install("VariantAnnotation")

# Option 2: vcfR
install.packages("vcfR")
```

## Validation Phases

- **Phase 1** (Current): Core functionality, synthetic data, BND parsing
- **Phase 2** (Planned): PCAWG dataset comparison, statistical properties
- **Phase 3** (Planned): Cross-tool comparison, literature case studies
- **Phase 4** (Planned): Clinical data, biological validation

## Adding Custom Tests

To add your own validation tests:

```R
# Add to phase1_validation.R

test_my_custom_validation <- function() {
    cat("Running custom test...\n")

    # Your test code here
    # Use stopifnot() for assertions

    cat("✓ Custom test PASSED\n")
    return(TRUE)
}

# Add to test execution
run_test("My Custom Validation", test_my_custom_validation)
```

## Reporting Issues

If validation tests fail:

1. Check the error messages in console output
2. Review `phase1_validation_results.rds` for details
3. Report issues at: https://github.com/godkin1211/ShatterSeek/issues

Include:
- Full console output
- R session info: `sessionInfo()`
- Test environment (OS, R version)
- Sample data (if applicable)

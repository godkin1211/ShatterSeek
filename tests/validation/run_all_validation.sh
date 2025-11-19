#!/bin/bash
# Run complete Phase 1 validation suite for ShatterSeek Extended Edition

set -e  # Exit on error

echo "=========================================="
echo "ShatterSeek Extended Edition"
echo "Phase 1 Validation Suite"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Change to validation directory
cd "$(dirname "$0")"

echo "Working directory: $(pwd)"
echo ""

# Step 1: Quick check
echo "Step 1/3: Running quick check..."
echo "----------------------------------------"
if Rscript quick_check.R; then
    echo -e "${GREEN}✓ Quick check passed${NC}"
else
    echo -e "${RED}✗ Quick check failed${NC}"
    echo "Aborting validation suite"
    exit 1
fi
echo ""

# Step 2: Full validation
echo "Step 2/3: Running comprehensive validation..."
echo "----------------------------------------"
if Rscript phase1_validation.R; then
    echo -e "${GREEN}✓ Comprehensive validation passed${NC}"
else
    echo -e "${RED}✗ Comprehensive validation failed${NC}"
    echo "Check output above for details"
    exit 1
fi
echo ""

# Step 3: Generate report
echo "Step 3/3: Generating validation report..."
echo "----------------------------------------"
if Rscript generate_report.R; then
    echo -e "${GREEN}✓ Report generated${NC}"
else
    echo -e "${YELLOW}⚠ Report generation failed (non-critical)${NC}"
fi
echo ""

# Summary
echo "=========================================="
echo "Validation Complete"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - phase1_validation_results.rds (detailed results)"
echo "  - validation_report.html (HTML report)"
echo "  - validation_summary.txt (text summary)"
echo ""
echo -e "${GREEN}✓ All Phase 1 validation tests passed!${NC}"
echo ""
echo "Next steps:"
echo "  1. Open validation_report.html in a web browser"
echo "  2. Review validation_summary.txt"
echo "  3. Proceed to Phase 2 validation (PCAWG data)"
echo ""

#!/usr/bin/env Rscript
# Generate HTML validation report from test results

# Load results
if (!file.exists("phase1_validation_results.rds")) {
    cat("Error: phase1_validation_results.rds not found\n")
    cat("Run 'Rscript phase1_validation.R' first\n")
    quit(status = 1)
}

results <- readRDS("phase1_validation_results.rds")

# Generate HTML report
html <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<title>ShatterSeek Extended Edition - Validation Report</title>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }",
    ".container { max-width: 900px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }",
    "h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }",
    "h2 { color: #34495e; margin-top: 30px; }",
    ".summary { background: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }",
    ".pass { color: #27ae60; font-weight: bold; }",
    ".fail { color: #e74c3c; font-weight: bold; }",
    ".metric { display: inline-block; margin: 10px 20px; }",
    ".metric-label { font-weight: bold; color: #7f8c8d; }",
    ".metric-value { font-size: 24px; margin-left: 10px; }",
    ".test-detail { background: #fff; border-left: 4px solid #3498db; padding: 15px; margin: 15px 0; }",
    ".test-pass { border-left-color: #27ae60; }",
    ".test-fail { border-left-color: #e74c3c; }",
    "table { width: 100%; border-collapse: collapse; margin: 20px 0; }",
    "th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }",
    "th { background-color: #34495e; color: white; }",
    "tr:hover { background-color: #f5f5f5; }",
    ".timestamp { color: #95a5a6; font-size: 14px; }",
    "</style>",
    "</head>",
    "<body>",
    "<div class='container'>",
    sprintf("<h1>ShatterSeek Extended Edition v%s</h1>", results$version),
    "<h2>Phase 1 Validation Report</h2>",
    sprintf("<p class='timestamp'>Generated: %s</p>", format(results$timestamp, "%Y-%m-%d %H:%M:%S")),
    "",
    "<!-- Summary -->",
    "<div class='summary'>",
    "<h3>Summary</h3>",
    sprintf("<div class='metric'><span class='metric-label'>Total Tests:</span><span class='metric-value'>%d</span></div>",
            results$tests_total),
    sprintf("<div class='metric'><span class='metric-label'>Passed:</span><span class='metric-value pass'>%d</span></div>",
            results$tests_passed),
    sprintf("<div class='metric'><span class='metric-label'>Failed:</span><span class='metric-value fail'>%d</span></div>",
            results$tests_failed),
    sprintf("<div class='metric'><span class='metric-label'>Success Rate:</span><span class='metric-value'>%.1f%%</span></div>",
            100 * results$tests_passed / results$tests_total),
    "</div>",
    ""
)

# Test details
html <- c(html,
    "<h2>Test Results</h2>",
    "<table>",
    "<tr><th>Test Name</th><th>Status</th><th>Details</th></tr>"
)

for (test_name in names(results$details)) {
    test <- results$details[[test_name]]
    status_class <- ifelse(test$status == "PASS", "pass", "fail")
    status_symbol <- ifelse(test$status == "PASS", "✓", "✗")

    error_msg <- if (!is.null(test$error)) {
        sprintf("<br><small style='color: #e74c3c;'>Error: %s</small>", test$error)
    } else {
        ""
    }

    html <- c(html,
        "<tr>",
        sprintf("<td>%s</td>", test_name),
        sprintf("<td class='%s'>%s %s</td>", status_class, status_symbol, test$status),
        sprintf("<td>%s%s</td>", ifelse(test$status == "PASS", "All checks passed", "Test failed"), error_msg),
        "</tr>"
    )
}

html <- c(html,
    "</table>",
    ""
)

# Test descriptions
html <- c(html,
    "<h2>Test Descriptions</h2>",
    "",
    "<div class='test-detail test-pass'>",
    "<h3>1. Built-in Dataset (DO17373)</h3>",
    "<p><strong>Purpose:</strong> Validate chromothripsis detection on the example dataset from Cortes-Ciriano et al. (2020)</p>",
    "<p><strong>Checks:</strong></p>",
    "<ul>",
    "<li>SV and CNV data loading</li>",
    "<li>Chromothripsis detection (original implementation)</li>",
    "<li>Comprehensive chromoanagenesis detection (extended edition)</li>",
    "<li>Classification functionality</li>",
    "<li>Consistency between original and extended implementations</li>",
    "</ul>",
    "</div>",
    "",
    "<div class='test-detail test-pass'>",
    "<h3>2. Synthetic Data Generation & Detection</h3>",
    "<p><strong>Purpose:</strong> Test detection sensitivity on synthetic chromothripsis events with known characteristics</p>",
    "<p><strong>Checks:</strong></p>",
    "<ul>",
    "<li>Synthetic data generation with realistic properties</li>",
    "<li>Detection sensitivity (target: ≥70%)</li>",
    "<li>Handling of different SV types and CN patterns</li>",
    "<li>Robustness across 10 independent test cases</li>",
    "</ul>",
    "</div>",
    "",
    "<div class='test-detail test-pass'>",
    "<h3>3. BND/TRA ALT Field Parsing</h3>",
    "<p><strong>Purpose:</strong> Validate VCF breakend (BND) record parsing for translocation detection</p>",
    "<p><strong>Checks:</strong></p>",
    "<ul>",
    "<li>Parsing of VCF 4.2+ bracket notation</li>",
    "<li>All bracket formats: [chr:pos[, ]chr:pos], etc.</li>",
    "<li>Correct mate chromosome and position extraction</li>",
    "<li>Interchromosomal event identification</li>",
    "<li>TRA record creation from BND records</li>",
    "</ul>",
    "</div>",
    ""
)

# System information
html <- c(html,
    "<h2>System Information</h2>",
    "<table>",
    sprintf("<tr><td><strong>R Version</strong></td><td>%s</td></tr>", R.version.string),
    sprintf("<tr><td><strong>Platform</strong></td><td>%s</td></tr>", R.version$platform),
    sprintf("<tr><td><strong>ShatterSeek Version</strong></td><td>%s</td></tr>", results$version),
    "</table>",
    ""
)

# Conclusion
if (results$tests_failed == 0) {
    html <- c(html,
        "<div class='summary' style='background: #d4edda; border: 1px solid #c3e6cb;'>",
        "<h3 style='color: #155724;'>✓ All Tests Passed</h3>",
        "<p>ShatterSeek Extended Edition passed all Phase 1 validation tests. Core functionality is working correctly.</p>",
        "<p><strong>Next Steps:</strong></p>",
        "<ul>",
        "<li>Proceed to Phase 2: PCAWG dataset validation</li>",
        "<li>Test on real clinical data</li>",
        "<li>Compare with other tools</li>",
        "</ul>",
        "</div>"
    )
} else {
    html <- c(html,
        "<div class='summary' style='background: #f8d7da; border: 1px solid #f5c6cb;'>",
        "<h3 style='color: #721c24;'>✗ Some Tests Failed</h3>",
        "<p>Please review the failed tests above and address the issues before proceeding.</p>",
        "</div>"
    )
}

html <- c(html,
    "</div>",
    "</body>",
    "</html>"
)

# Write HTML file
writeLines(html, "validation_report.html")

cat("\n===========================================\n")
cat("Validation Report Generated\n")
cat("===========================================\n")
cat("File: validation_report.html\n")
cat("Open this file in a web browser to view the full report\n\n")

# Also create a simple text summary
summary_file <- "validation_summary.txt"
sink(summary_file)
cat("ShatterSeek Extended Edition - Phase 1 Validation Summary\n")
cat(strrep("=", 60), "\n\n")
cat(sprintf("Timestamp: %s\n", format(results$timestamp, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Version: %s\n\n", results$version))
cat("Test Results:\n")
cat(sprintf("  Total: %d\n", results$tests_total))
cat(sprintf("  Passed: %d\n", results$tests_passed))
cat(sprintf("  Failed: %d\n", results$tests_failed))
cat(sprintf("  Success Rate: %.1f%%\n\n", 100 * results$tests_passed / results$tests_total))

if (results$tests_failed > 0) {
    cat("Failed Tests:\n")
    for (test_name in names(results$details)) {
        if (results$details[[test_name]]$status == "FAIL") {
            cat(sprintf("  - %s\n", test_name))
        }
    }
} else {
    cat("✓ All tests passed!\n")
}
sink()

cat(sprintf("Text summary: %s\n\n", summary_file))

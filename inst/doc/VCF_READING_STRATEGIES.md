# VCF Reading Strategies for ShatterSeek

## Overview

ShatterSeek Extended Edition (v2.0.0+) provides **three different ways** to read SV VCF files, each with different trade-offs. This guide helps you choose the right approach.

## Strategy Comparison

| Feature | Built-in VCF Readers | StructuralVariantAnnotation | Manual BEDPE |
|---------|---------------------|----------------------------|--------------|
| **Complexity** | Low (one function call) | Medium (two steps) | High (manual conversion) |
| **Robustness** | Good for standard callers | Excellent for complex VCFs | Depends on implementation |
| **Flexibility** | Caller-specific optimizations | Highly flexible | Full control |
| **Dependencies** | VariantAnnotation OR vcfR | StructuralVariantAnnotation | Custom code |
| **Speed** | Fast | Moderate | Depends |
| **Metadata preservation** | High | Medium | Low |

---

## Strategy 1: Built-in VCF Readers (Recommended for Most Users)

### When to Use
- Standard SV callers (Manta, DRAGEN, Delly, GRIDSS, LUMPY, Sniffles)
- You want the simplest workflow
- You trust caller-specific parsing optimizations

### Workflow

```R
library(ShatterSeek)

# One-step VCF reading
sv_data <- read_sv_vcf("sample.manta.vcf.gz", caller = "manta")
cn_data <- read_cnv_vcf("sample.cnvkit.vcf.gz", caller = "cnvkit")

# Immediate analysis
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")

# Generate reports
plot_genome_dashboard(results, sv_data, cn_data, "Sample_01")
```

### Supported Callers

**SV callers:**
- Manta (Illumina)
- DRAGEN (Illumina)
- GRIDSS (complex event caller)
- Delly
- LUMPY
- Sniffles (long-read)

**CNV callers:**
- CNVkit
- GATK
- Control-FREEC
- Canvas
- DRAGEN

### Pros
- ✅ Simplest workflow (one function call)
- ✅ Caller-specific optimizations
- ✅ Automatic BND parsing
- ✅ Preserves VCF metadata
- ✅ Integrated error handling

### Cons
- ⚠️ Limited to supported callers (though auto-detection fallback exists)
- ⚠️ Newer implementation (less battle-tested than StructuralVariantAnnotation)

---

## Strategy 2: StructuralVariantAnnotation (Recommended for Complex Cases)

### When to Use
- Complex or unusual VCF formats
- Unsupported SV callers
- You need additional SV QC and annotation
- You're already using Bioconductor ecosystem
- Maximum robustness for BND parsing

### Workflow

```R
library(ShatterSeek)
library(StructuralVariantAnnotation)
library(VariantAnnotation)

# Step 1: Read VCF with StructuralVariantAnnotation
vcf <- readVcf("sample.vcf.gz", "hg38")

# Step 2a: Extract breakpoints
gr <- breakpointRanges(vcf)

# Step 2b: Optionally include unpaired breakends
gr_all <- c(breakpointRanges(vcf), breakendRanges(vcf))

# Step 3: Convert to ShatterSeek format
sv_data <- granges_to_svs(gr_all)

# Step 4: Continue with ShatterSeek analysis
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
```

### Alternative: One-Step Wrapper

```R
# ShatterSeek v2.0.2+ provides a convenience wrapper
sv_data <- read_sv_vcf_structuralvariant("sample.vcf.gz", genome = "hg38")

# Then analyze as usual
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
```

### Pros
- ✅ Most robust BND parsing (well-tested Bioconductor package)
- ✅ Handles complex SV types automatically
- ✅ Excellent for unusual VCF formats
- ✅ Access to additional StructuralVariantAnnotation features
- ✅ GRanges format enables downstream Bioconductor analysis

### Cons
- ⚠️ Extra dependency (StructuralVariantAnnotation)
- ⚠️ Two-step process (unless using wrapper)
- ⚠️ May lose some VCF metadata in conversion

---

## Strategy 3: Manual BEDPE Conversion (Advanced Users)

### When to Use
- You already have BEDPE files
- Custom SV caller with non-standard VCF
- You need full control over data processing
- Integrating with existing pipeline that outputs BEDPE

### Workflow (Example from GitHub)

```R
library(StructuralVariantAnnotation)
library(VariantAnnotation)

# Step 1: Read VCF
vcf <- readVcf("PASS_somaticSV.vcf", "hg38")

# Step 2: Extract breakpoints
gr <- c(breakpointRanges(vcf), breakendRanges(vcf))

# Step 3: Convert to BEDPE (custom function needed)
bedpe <- breakpointgr2bedpe(gr)  # User-defined function

# Step 4: Read BEDPE into ShatterSeek
sv_data <- read_bedpe_to_svs(bedpe)  # Custom conversion

# Step 5: Analysis
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
```

### Pros
- ✅ Full control over conversion process
- ✅ Can handle any custom format
- ✅ BEDPE is widely supported format

### Cons
- ⚠️ Most complex workflow
- ⚠️ Need to write custom conversion functions
- ⚠️ Risk of errors in manual conversion
- ⚠️ May lose important VCF information
- ⚠️ Strand inference can be error-prone

### BEDPE to SVs Conversion Template

If you have BEDPE files, here's a template for conversion:

```R
#' Convert BEDPE to ShatterSeek SVs object
#'
#' @param bedpe_file Path to BEDPE file
#' @return SVs object
read_bedpe_to_svs <- function(bedpe_file) {

    # Read BEDPE (standard format: chrom1, start1, end1, chrom2, start2, end2, ...)
    bedpe <- read.table(bedpe_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

    # Extract required fields
    chrom1 <- bedpe$chrom1
    pos1 <- bedpe$start1  # or end1 depending on convention
    chrom2 <- bedpe$chrom2
    pos2 <- bedpe$start2  # or end2

    # Infer SV type from strands or use provided column
    if ("svtype" %in% names(bedpe)) {
        SVtype <- bedpe$svtype
    } else {
        # Infer from chromosome and strands
        SVtype <- ifelse(chrom1 != chrom2, "TRA",
                        ifelse(bedpe$strand1 == "+" & bedpe$strand2 == "-", "DEL",
                        ifelse(bedpe$strand1 == "-" & bedpe$strand2 == "+", "DUP",
                        ifelse(bedpe$strand1 == "+" & bedpe$strand2 == "+", "h2hINV",
                        ifelse(bedpe$strand1 == "-" & bedpe$strand2 == "-", "t2tINV",
                               "UNKNOWN")))))
    }

    # Get strands
    strand1 <- bedpe$strand1
    strand2 <- bedpe$strand2

    # Create SVs object
    sv_data <- SVs(
        chrom1 = chrom1,
        pos1 = pos1,
        chrom2 = chrom2,
        pos2 = pos2,
        SVtype = SVtype,
        strand1 = strand1,
        strand2 = strand2
    )

    return(sv_data)
}
```

---

## Decision Tree

```
Do you have a VCF file?
│
├─ YES → Is your caller supported by read_sv_vcf()?
│         │
│         ├─ YES → Use Strategy 1 (Built-in VCF Readers) ← EASIEST
│         │
│         └─ NO → Is the VCF format complex or unusual?
│                 │
│                 ├─ YES → Use Strategy 2 (StructuralVariantAnnotation) ← MOST ROBUST
│                 │
│                 └─ NO → Try Strategy 1 with caller="auto" first
│
└─ NO → Do you have BEDPE files?
         │
         ├─ YES → Use Strategy 3 (Manual BEDPE) ← ADVANCED
         │
         └─ NO → You need VCF or BEDPE format
                  Contact your SV caller documentation
```

---

## Real-World Examples

### Example 1: Standard Manta Analysis (Strategy 1)

```R
library(ShatterSeek)

# Simplest approach
sv_data <- read_sv_vcf("tumor.manta.vcf.gz")
cn_data <- read_cnv_vcf("tumor.cnvkit.vcf.gz")

results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
classify_mixed_mechanisms(results)
```

### Example 2: Complex GRIDSS VCF (Strategy 2)

```R
library(ShatterSeek)
library(StructuralVariantAnnotation)
library(VariantAnnotation)

# GRIDSS produces complex BND records - use StructuralVariantAnnotation
vcf <- readVcf("tumor.gridss.vcf.gz", "hg38")

# StructuralVariantAnnotation handles complex pairing
gr <- c(breakpointRanges(vcf), breakendRanges(vcf))

# Convert to ShatterSeek
sv_data <- granges_to_svs(gr)

# Or use convenience wrapper (v2.0.2+)
sv_data <- read_sv_vcf_structuralvariant("tumor.gridss.vcf.gz", genome = "hg38")

# Continue analysis
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
```

### Example 3: Custom Caller with BEDPE Output (Strategy 3)

```R
library(ShatterSeek)

# You have custom BEDPE from your pipeline
bedpe <- read.table("tumor.custom_caller.bedpe", header = TRUE, sep = "\t")

# Convert to ShatterSeek format (use template above)
sv_data <- read_bedpe_to_svs("tumor.custom_caller.bedpe")

# Analyze
results <- detect_chromoanagenesis(sv_data, cn_data, genome = "hg38")
```

---

## Troubleshooting

### Problem: "No SVs found" with Strategy 1

**Solution**: Try Strategy 2 (StructuralVariantAnnotation)
- Your VCF may have unusual format
- BND records may not be standard

### Problem: "Partner not found" with Strategy 2

**Solution**: Check if breakends are unpaired
```R
# Separate paired and unpaired
gr_paired <- breakpointRanges(vcf)
gr_unpaired <- breakendRanges(vcf)

# Analyze paired only
sv_data <- granges_to_svs(gr_paired)
```

### Problem: Strand errors with Strategy 3

**Solution**: Verify BEDPE strand encoding
- ShatterSeek expects: DEL (+/-), DUP (-/+), h2hINV (+/+), t2tINV (-/-)
- Some BEDPE formats use different conventions

### Problem: Performance issues

**Solution**: Use Strategy 1 for large VCFs
- Built-in readers are optimized for speed
- StructuralVariantAnnotation is more thorough but slower

---

## Best Practices

1. **Start with Strategy 1** for standard callers
   - Fastest and simplest
   - Good for most use cases

2. **Use Strategy 2** when Strategy 1 fails or for complex VCFs
   - More robust BND parsing
   - Better for unusual formats

3. **Reserve Strategy 3** for custom pipelines
   - Only when you already have BEDPE
   - Verify strand encoding carefully

4. **Always validate your data**
   ```R
   # Check your SVs object
   summary(sv_data)

   # Verify SV types distribution
   table(sv_data@SVtype)

   # Check for strand consistency
   data.frame(
       SVtype = sv_data@SVtype,
       strand1 = sv_data@strand1,
       strand2 = sv_data@strand2
   ) %>% table()
   ```

5. **Use diagnostic tools**
   ```R
   # Run quality check
   qc <- check_data_quality(sv_data, cn_data)
   print(qc)
   ```

---

## Recommendation Summary

| Your Situation | Recommended Strategy | Why |
|----------------|---------------------|-----|
| Manta, Delly, LUMPY VCF | **Strategy 1** | Built-in support, simplest |
| GRIDSS, complex VCF | **Strategy 2** | More robust BND parsing |
| Unknown/custom caller | **Try 1 → 2** | Auto-detect, fallback to robust |
| Already have BEDPE | **Strategy 3** | Direct conversion |
| Research/exploration | **Strategy 2** | Access to full Bioconductor tools |
| Production pipeline | **Strategy 1** | Fast, integrated, caller-optimized |

---

## References

1. **StructuralVariantAnnotation**
   - Bioconductor package: https://bioconductor.org/packages/StructuralVariantAnnotation/
   - Paper: Cameron et al. (2017) doi:10.1371/journal.pcbi.1005595

2. **ShatterSeek Extended Edition**
   - GitHub: https://github.com/godkin1211/ShatterSeek
   - Original paper: Cortes-Ciriano et al. (2020) Nat. Genet. 52, 331-341

3. **VCF Specification**
   - VCF 4.2: https://samtools.github.io/hts-specs/VCFv4.2.pdf
   - BND notation for structural variants

---

## Support

For questions about VCF reading strategies:
- GitHub Issues: https://github.com/godkin1211/ShatterSeek/issues
- Email: n28111021@gs.ncku.edu.tw

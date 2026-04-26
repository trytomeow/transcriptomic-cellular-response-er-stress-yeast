#!/bin/bash
# ==============================================================================
# Upstream Step 5: Variant Calling
# ==============================================================================
# Script: 05_variant_calling.sh
# Output: VCF files in 05_variants/
# Description: SNV/InDel calling using BCFtools.
#              Filter: QUAL >= 20, DP >= 10.
# ==============================================================================

set -euo pipefail

# Source configuration (only if not already loaded by master)
if [[ -z "${PROJECT_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config.sh"
fi

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

run_variant_calling() {
    log_message "=== STEP 5: Variant Calling (BCFtools) ==="
    log_message "  Filter: QUAL >= ${VCF_MIN_QUAL}, DP >= ${VCF_MIN_DEPTH}"
    
    local bam_files
    bam_files=$(find "${SORTED_DIR}" -name "*_sorted.bam" | sort)
    
    if [[ -z "${bam_files}" ]]; then
        log_message "ERROR: No sorted BAM files found in ${SORTED_DIR}"
        exit 1
    fi
    
    while IFS= read -r bam_file; do
        local sample_name
        sample_name=$(basename "${bam_file}" "_sorted.bam")
        local raw_vcf="${VCF_DIR}/${sample_name}_raw.vcf.gz"
        local filtered_vcf="${VCF_DIR}/${sample_name}_filtered.vcf.gz"
        
        log_message "Calling variants for ${sample_name}..."
        
        # Step 1: mpileup -> call
        bcftools mpileup \
            -Ou \
            -f "${REFERENCE_FASTA}" \
            --min-BQ "${MIN_BASE_QUALITY}" \
            --min-MQ "${MIN_MAPPING_QUALITY}" \
            --threads "${THREADS}" \
            "${bam_file}" | \
        bcftools call \
            -mv \
            -Oz \
            -o "${raw_vcf}" \
            --threads "${THREADS}"
        
        # Index raw VCF
        bcftools index "${raw_vcf}"
        
        # Step 2: Filter variants
        bcftools filter \
            -i "${VCF_FILTER_EXPR}" \
            -Oz \
            -o "${filtered_vcf}" \
            "${raw_vcf}"
        
        # Index filtered VCF
        bcftools index "${filtered_vcf}"
        
        # Count variants
        local raw_count
        raw_count=$(bcftools view -H "${raw_vcf}" | wc -l)
        local filtered_count
        filtered_count=$(bcftools view -H "${filtered_vcf}" | wc -l)
        
        log_message "  ✓ ${sample_name}: ${raw_count} raw -> ${filtered_count} filtered variants"
        
    done <<< "${bam_files}"
    
    log_message "✓ Variant calling complete"
    log_message "  Output: ${VCF_DIR}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${VCF_DIR}"
    run_variant_calling
fi

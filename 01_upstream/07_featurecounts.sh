#!/bin/bash
# ==============================================================================
# Upstream Step 7: Gene Quantification
# ==============================================================================
# Script: 07_featurecounts.sh
# Output: Count matrices in 07_counts/
# Description: Gene and exon quantification using FeatureCounts.
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

run_featurecounts() {
    log_message "=== STEP 7: FeatureCounts (Gene Quantification) ==="
    
    local bam_files
    bam_files=$(find "${SORTED_DIR}" -name "*_sorted.bam" | sort | tr '\n' ' ')
    
    if [[ -z "${bam_files}" ]]; then
        log_message "ERROR: No sorted BAM files found in ${SORTED_DIR}"
        exit 1
    fi
    
    log_message "Counting reads for all samples..."
    
    # Gene-level counts
    featureCounts \
        -T "${THREADS}" \
        -a "${GTF_FILE}" \
        -o "${COUNT_MATRIX_FILE}" \
        -g "${GENE_ID_ATTR}" \
        -t "${FEATURE_TYPE}" \
        -s "${STRAND_SPECIFIC}" \
        ${FEATURECOUNTS_OPTS} \
        ${bam_files} \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ Gene-level counts complete"
    log_message "  Count matrix: ${COUNT_MATRIX_FILE}"
    
    # Exon-level counts for DEXSeq
    log_message "Counting reads at exon level for DEXSeq..."
    
    featureCounts \
        -T "${THREADS}" \
        -a "${GTF_FILE}" \
        -o "${EXON_COUNT_MATRIX_FILE}" \
        -g "${GENE_ID_ATTR}" \
        -t "${EXON_FEATURE_TYPE}" \
        -s "${STRAND_SPECIFIC}" \
        -f \
        -O \
        ${FEATURECOUNTS_OPTS} \
        ${bam_files} \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ Exon-level counts complete"
    log_message "  Exon count matrix: ${EXON_COUNT_MATRIX_FILE}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${COUNTS_DIR}"
    run_featurecounts
fi

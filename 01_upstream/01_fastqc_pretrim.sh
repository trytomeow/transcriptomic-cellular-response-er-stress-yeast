#!/bin/bash
# ==============================================================================
# Upstream Step 1: Pre-trimming FastQC
# ==============================================================================
# Script: 01_fastqc_pretrim.sh
# Output: FASTQC reports in 01_fastqc_pretrim/
# Description: Pre-trimming quality control using FastQC.
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

run_fastqc_pretrim() {
    log_message "=== STEP 1: FastQC (Pre-trimming) ==="
    
    local fastq_files
    fastq_files=$(find "${RAW_DATA_DIR}" -type f \( -name "*.fq" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fastq.gz" \) | sort)
    
    if [[ -z "${fastq_files}" ]]; then
        log_message "ERROR: No FASTQ files found in ${RAW_DATA_DIR}"
        exit 1
    fi
    
    local file_count
    file_count=$(echo "${fastq_files}" | wc -l)
    log_message "Found ${file_count} FASTQ files"
    
    echo "${fastq_files}" | xargs fastqc \
        --threads "${THREADS}" \
        --outdir "${FASTQC_PRE_DIR}" \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ FastQC (pre-trim) complete"
    log_message "  Output: ${FASTQC_PRE_DIR}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${FASTQC_PRE_DIR}"
    run_fastqc_pretrim
fi

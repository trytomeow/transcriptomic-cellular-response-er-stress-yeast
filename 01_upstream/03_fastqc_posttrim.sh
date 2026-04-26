#!/bin/bash
# ==============================================================================
# Upstream Step 3: Post-trimming FastQC
# ==============================================================================
# Script: 03_fastqc_posttrim.sh
# Output: FASTQC reports in 03_fastqc_posttrim/
# Description: Post-trimming quality control using FastQC.
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

run_fastqc_posttrim() {
    log_message "=== STEP 3: FastQC (Post-trimming) ==="
    
    local trimmed_files
    trimmed_files=$(find "${TRIMMED_DIR}" -name "*_paired.fq.gz" | sort)
    
    if [[ -z "${trimmed_files}" ]]; then
        log_message "ERROR: No trimmed files found in ${TRIMMED_DIR}"
        exit 1
    fi
    
    local file_count
    file_count=$(echo "${trimmed_files}" | wc -l)
    log_message "Found ${file_count} trimmed files"
    
    echo "${trimmed_files}" | xargs fastqc \
        --threads "${THREADS}" \
        --outdir "${FASTQC_POST_DIR}" \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ FastQC (post-trim) complete"
    log_message "  Output: ${FASTQC_POST_DIR}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${FASTQC_POST_DIR}"
    run_fastqc_posttrim
fi

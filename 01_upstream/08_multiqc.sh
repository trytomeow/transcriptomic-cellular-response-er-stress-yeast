#!/bin/bash
# ==============================================================================
# Upstream Step 8: Quality Report Aggregation
# ==============================================================================
# Script: 08_multiqc.sh
# Output: MultiQC report in 08_multiqc/
# Description: Aggregate QC reports using MultiQC.
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

run_multiqc() {
    log_message "=== STEP 8: MultiQC (Aggregate Reports) ==="
    
    multiqc \
        "${RESULTS_DIR}" \
        -o "${MULTIQC_DIR}" \
        -n "multiqc_report" \
        --force \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ MultiQC complete"
    log_message "  Report: ${MULTIQC_DIR}/multiqc_report.html"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${MULTIQC_DIR}"
    run_multiqc
fi

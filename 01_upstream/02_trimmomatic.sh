#!/bin/bash
# ==============================================================================
# Upstream Step 2: Read Trimming
# ==============================================================================
# Script: 02_trimmomatic.sh
# Output: Trimmed FASTQ files in 02_trimmed/
# Description: Read trimming using Trimmomatic (PE mode, Phred 20, MinLen 75).
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

run_trimmomatic() {
    log_message "=== STEP 2: Trimmomatic (PE Mode, Phred 20, MinLen 75) ==="
    
    local r1_files
    r1_files=$(find "${RAW_DATA_DIR}" -type f \( -name "*${READ1_PATTERN}*.fq" -o -name "*${READ1_PATTERN}*.fq.gz" \) | sort)
    
    if [[ -z "${r1_files}" ]]; then
        log_message "ERROR: No R1 files found matching pattern ${READ1_PATTERN}"
        exit 1
    fi
    
    while IFS= read -r r1_file; do
        local sample_name
        sample_name=$(get_sample_name "${r1_file}")
        log_message "Processing sample: ${sample_name}"
        
        local r2_file="${r1_file/${READ1_PATTERN}/${READ2_PATTERN}}"
        
        if [[ ! -f "${r2_file}" ]]; then
            log_message "WARNING: R2 file not found for ${sample_name}, skipping..."
            continue
        fi
        
        local out_r1_paired="${TRIMMED_DIR}/${sample_name}_R1_paired.fq.gz"
        local out_r1_unpaired="${TRIMMED_DIR}/${sample_name}_R1_unpaired.fq.gz"
        local out_r2_paired="${TRIMMED_DIR}/${sample_name}_R2_paired.fq.gz"
        local out_r2_unpaired="${TRIMMED_DIR}/${sample_name}_R2_unpaired.fq.gz"
        
        trimmomatic PE \
            -threads "${THREADS}" \
            -phred33 \
            "${r1_file}" "${r2_file}" \
            "${out_r1_paired}" "${out_r1_unpaired}" \
            "${out_r2_paired}" "${out_r2_unpaired}" \
            "${ILLUMINACLIP}" \
            LEADING:"${LEADING}" \
            TRAILING:"${TRAILING}" \
            SLIDINGWINDOW:"${SLIDINGWINDOW}" \
            MINLEN:"${MINLEN}" \
            >> "${LOG_FILE}" 2>&1
        
        log_message "  ✓ ${sample_name} trimmed"
        
    done <<< "${r1_files}"
    
    log_message "✓ Trimmomatic complete"
    log_message "  Output: ${TRIMMED_DIR}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${TRIMMED_DIR}"
    run_trimmomatic
fi

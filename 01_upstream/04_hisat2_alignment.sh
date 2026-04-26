#!/bin/bash
# ==============================================================================
# Upstream Step 4: RNA-Seq Alignment
# ==============================================================================
# Script: 04_hisat2_alignment.sh
# Output: Sorted BAM files in 04_sorted_bam/
# Description: HISAT2 alignment piped directly to sorted BAM (no intermediate SAM).
# ==============================================================================

# Source configuration (only if not already loaded by master)
if [[ -z "${PROJECT_DIR:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config.sh"
fi

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

run_hisat2() {
    log_message "=== STEP 4: HISAT2 Alignment -> Sorted BAM (Piped, No SAM) ==="
    
    local r1_files
    r1_files=$(find "${TRIMMED_DIR}" -name "*_R1_paired.fq.gz" | sort)
    
    if [[ -z "${r1_files}" ]]; then
        log_message "ERROR: No trimmed R1 files found in ${TRIMMED_DIR}"
        exit 1
    fi
    
    local total_samples
    total_samples=$(echo "${r1_files}" | wc -l)
    local current=0
    
    while IFS= read -r r1_file; do
        current=$((current + 1))
        local sample_name
        sample_name=$(basename "${r1_file}" "_R1_paired.fq.gz")
        local r2_file="${TRIMMED_DIR}/${sample_name}_R2_paired.fq.gz"
        local sorted_bam="${SORTED_DIR}/${sample_name}_sorted.bam"
        local summary_file="${SORTED_DIR}/${sample_name}_hisat2_summary.txt"
        
        if [[ ! -f "${r2_file}" ]]; then
            log_message "WARNING: R2 file not found for ${sample_name}, skipping..."
            continue
        fi
        
        log_message "[${current}/${total_samples}] Aligning ${sample_name} (piped to sorted BAM)..."
        
        # OPTIMIZED PIPELINE: HISAT2 -> samtools view -> samtools sort -> BAM
        # Use set +e temporarily so PIPESTATUS is captured before set -e can exit
        set +e
        hisat2 \
            -p "${THREADS}" \
            -x "${HISAT2_INDEX}" \
            -1 "${r1_file}" \
            -2 "${r2_file}" \
            ${HISAT2_OPTS} \
            --new-summary \
            --summary-file "${summary_file}" \
            2>>"${LOG_FILE}" \
        | samtools view \
            -@ "${THREADS}" \
            -bS \
            - \
        | samtools sort \
            -@ "${THREADS}" \
            -m "${SORT_MEMORY}" \
            -o "${sorted_bam}" \
            -
        local pipe_status=("${PIPESTATUS[@]}")
        set -e

        if [[ ${pipe_status[0]} -ne 0 ]]; then
            log_message "  ✗ ERROR: HISAT2 alignment failed for ${sample_name} (exit: ${pipe_status[0]})"
            exit 1
        elif [[ ${pipe_status[1]} -ne 0 ]]; then
            log_message "  ✗ ERROR: samtools view failed for ${sample_name} (exit: ${pipe_status[1]})"
            exit 1
        elif [[ ${pipe_status[2]} -ne 0 ]]; then
            log_message "  ✗ ERROR: samtools sort failed for ${sample_name} (exit: ${pipe_status[2]})"
            exit 1
        elif [[ ! -f "${sorted_bam}" ]]; then
            log_message "  ✗ ERROR: Sorted BAM not created for ${sample_name}"
            exit 1
        fi

        samtools index -@ "${THREADS}" "${sorted_bam}"
        log_message "  ✓ ${sample_name}: aligned, sorted, and indexed"
        
    done <<< "${r1_files}"
    
    log_message "✓ HISAT2 alignment complete"
    log_message "  Output: ${SORTED_DIR}"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${SORTED_DIR}"
    run_hisat2
fi

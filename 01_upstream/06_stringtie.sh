#!/bin/bash
# ==============================================================================
# Upstream Step 6: Transcript Assembly
# ==============================================================================
# Script: 06_stringtie.sh
# Output: GTF files and annotations in 06_stringtie/
# Description: Novel transcript assembly using StringTie.
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

run_stringtie() {
    log_message "=== STEP 6: StringTie (Novel Transcript Assembly) ==="
    
    local bam_files
    bam_files=$(find "${SORTED_DIR}" -name "*_sorted.bam" | sort)
    local gtf_list="${STRINGTIE_DIR}/gtf_list.txt"
    
    if [[ -z "${bam_files}" ]]; then
        log_message "ERROR: No sorted BAM files found in ${SORTED_DIR}"
        exit 1
    fi
    
    # Clear previous list
    > "${gtf_list}"
    
    while IFS= read -r bam_file; do
        local sample_name
        sample_name=$(basename "${bam_file}" "_sorted.bam")
        local sample_gtf="${STRINGTIE_DIR}/${sample_name}.gtf"
        
        log_message "Assembling transcripts for ${sample_name}..."
        
        stringtie \
            -p "${THREADS}" \
            -G "${GTF_FILE}" \
            -o "${sample_gtf}" \
            ${STRINGTIE_OPTS} \
            "${bam_file}" \
            >> "${LOG_FILE}" 2>&1
        
        echo "${sample_gtf}" >> "${gtf_list}"
        log_message "  ✓ ${sample_name} assembled"
        
    done <<< "${bam_files}"
    
    # Merge all GTFs
    log_message "Merging all assembled GTFs..."
    stringtie --merge \
        -p "${THREADS}" \
        -G "${GTF_FILE}" \
        -o "${MERGED_GTF}" \
        "${gtf_list}" \
        >> "${LOG_FILE}" 2>&1
    
    # Compare with reference
    log_message "Comparing with reference annotation..."
    gffcompare \
        -r "${GTF_FILE}" \
        -o "${STRINGTIE_DIR}/gffcompare" \
        "${MERGED_GTF}" \
        >> "${LOG_FILE}" 2>&1
    
    log_message "✓ StringTie complete"
    log_message "  Merged GTF: ${MERGED_GTF}"
    log_message "  Novel transcripts: ${STRINGTIE_DIR}/gffcompare.annotated.gtf"
}

# ==============================================================================
# EXECUTION
# ==============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    mkdir -p "${STRINGTIE_DIR}"
    run_stringtie
fi

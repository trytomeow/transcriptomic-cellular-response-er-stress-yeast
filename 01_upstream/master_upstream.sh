#!/bin/bash
# ==============================================================================
# Upstream: Master Pipeline Script
# ==============================================================================
# Script: master_upstream.sh
# Output: (Executes everything)
# Description: Master script to orchestrate all upstream pipeline modules.
#              Use flags to skip specific steps if needed.
# ==============================================================================

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source configuration
source "${SCRIPT_DIR}/00_config.sh"

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================

SKIP_QC=false
SKIP_VARIANTS=false
ONLY_COUNTS=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-qc)
            SKIP_QC=true
            shift
            ;;
        --skip-variants)
            SKIP_VARIANTS=true
            shift
            ;;
        --only-counts)
            ONLY_COUNTS=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main() {
    log_message "=============================================="
    log_message "  RNA-Seq Upstream Pipeline (Modular)"
    log_message "  Organism: S. cerevisiae (BY4742)"
    log_message "  Threads: ${THREADS}"
    log_message "=============================================="
    
    # Create directories
    create_directories
    
    # Check dependencies and print versions
    check_dependencies
    print_tool_versions
    check_reference_files
    
    # ==================================================================
    # STEP 1: FastQC (Pre-trimming)
    # ==================================================================
    if [[ "${SKIP_QC}" == "false" ]]; then
        source "${SCRIPT_DIR}/01_fastqc_pretrim.sh"
        run_fastqc_pretrim
    else
        log_message "⏭ Skipping pre-trim FastQC (--skip-qc)"
    fi
    
    # ==================================================================
    # STEP 2: Trimmomatic
    # ==================================================================
    source "${SCRIPT_DIR}/02_trimmomatic.sh"
    run_trimmomatic
    
    # ==================================================================
    # STEP 3: FastQC (Post-trimming)
    # ==================================================================
    if [[ "${SKIP_QC}" == "false" ]]; then
        source "${SCRIPT_DIR}/03_fastqc_posttrim.sh"
        run_fastqc_posttrim
    else
        log_message "⏭ Skipping post-trim FastQC (--skip-qc)"
    fi
    
    # ==================================================================
    # STEP 4: HISAT2 Alignment -> Sorted BAM
    # ==================================================================
    source "${SCRIPT_DIR}/04_hisat2_alignment.sh"
    run_hisat2
    
    # ==================================================================
    # STEP 5: Variant Calling (Optional)
    # ==================================================================
    if [[ "${SKIP_VARIANTS}" == "false" ]] && [[ "${ONLY_COUNTS}" == "false" ]]; then
        source "${SCRIPT_DIR}/05_variant_calling.sh"
        run_variant_calling
    else
        log_message "⏭ Skipping variant calling"
    fi
    
    # ==================================================================
    # STEP 6: StringTie (Optional)
    # ==================================================================
    if [[ "${ONLY_COUNTS}" == "false" ]]; then
        source "${SCRIPT_DIR}/06_stringtie.sh"
        run_stringtie
    else
        log_message "⏭ Skipping StringTie"
    fi
    
    # ==================================================================
    # STEP 7: FeatureCounts
    # ==================================================================
    source "${SCRIPT_DIR}/07_featurecounts.sh"
    run_featurecounts
    
    # ==================================================================
    # STEP 8: MultiQC
    # ==================================================================
    source "${SCRIPT_DIR}/08_multiqc.sh"
    run_multiqc
    
    # ==================================================================
    # COMPLETION
    # ==================================================================
    log_message "=============================================="
    log_message "  PIPELINE COMPLETE!"
    log_message "=============================================="
    log_message "Outputs:"
    log_message "  - Pre-trim FastQC:   ${FASTQC_PRE_DIR}"
    log_message "  - Trimmed:           ${TRIMMED_DIR}"
    log_message "  - Post-trim FastQC:  ${FASTQC_POST_DIR}"
    log_message "  - HISAT2 Alignments: ${SORTED_DIR}"
    log_message "  - Variants:          ${VCF_DIR}"
    log_message "  - StringTie GTFs:    ${STRINGTIE_DIR}"
    log_message "  - Counts:            ${COUNTS_DIR}"
    log_message "  - MultiQC:           ${MULTIQC_DIR}"
    log_message "  - Logs:              ${LOGS_DIR}"
}

# Run main
main "$@"

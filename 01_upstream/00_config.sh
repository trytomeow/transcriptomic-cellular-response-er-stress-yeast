#!/bin/bash
# ==============================================================================
# Upstream: Master Configuration File
# ==============================================================================
# Script: 00_config.sh
# Output: (none — sourced by all upstream scripts)
# Description: Shared configuration file for upstream pipeline modules.
#              All paths, parameters, and helper functions defined here.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# █▀▀ █▀█ █▄ █ █▀▀ █ █▀▀ █ █ █▀█ ▄▀█ ▀█▀ █ █▀█ █▄ █   █▀ █▀▀ █▀▀ ▀█▀ █ █▀█ █▄ █
# █▄▄ █▄█ █ ▀█ █▀  █ █▄█ █▄█ █▀▄ █▀█  █  █ █▄█ █ ▀█   ▄█ ██▄ █▄▄  █  █ █▄█ █ ▀█
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: SYSTEM RESOURCES
# ------------------------------------------------------------------------------
THREADS=16

# ------------------------------------------------------------------------------
# SECTION 2: PROJECT DIRECTORIES
# ------------------------------------------------------------------------------
# Detect repo root from this file's location:
# 00_config.sh is at <repo_root>/01_upstream/ -> one dirname() up = repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Input directories
RAW_DATA_DIR="${PROJECT_DIR}/data/raw"    # Raw FASTQ files
REFS_DIR="${PROJECT_DIR}/refs"            # Reference genome, HISAT2 index, GTF

# Output base directory
RESULTS_DIR="${PROJECT_DIR}/results/upstream"

# Output subdirectories
FASTQC_PRE_DIR="${RESULTS_DIR}/01_fastqc_pretrim"
TRIMMED_DIR="${RESULTS_DIR}/02_trimmed"
FASTQC_POST_DIR="${RESULTS_DIR}/03_fastqc_posttrim"
SORTED_DIR="${RESULTS_DIR}/04_sorted_bam"
VCF_DIR="${RESULTS_DIR}/05_variants"
STRINGTIE_DIR="${RESULTS_DIR}/06_stringtie"
COUNTS_DIR="${RESULTS_DIR}/07_counts"
MULTIQC_DIR="${RESULTS_DIR}/08_multiqc"
LOGS_DIR="${RESULTS_DIR}/logs"

# ------------------------------------------------------------------------------
# SECTION 3: REFERENCE GENOME FILES
# ------------------------------------------------------------------------------
HISAT2_INDEX="${REFS_DIR}/hisat2_index/r64_tran/genome_tran"
REFERENCE_FASTA="${REFS_DIR}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
GTF_FILE="${REFS_DIR}/Saccharomyces_cerevisiae.R64.gtf"
GENE_ID_ATTR="gene_id"
FEATURE_TYPE="exon"

# ------------------------------------------------------------------------------
# SECTION 4: TRIMMOMATIC SETTINGS
# ------------------------------------------------------------------------------
TRIMMOMATIC_JAR="${CONDA_PREFIX}/share/trimmomatic/trimmomatic.jar"
ADAPTER_FILE="${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
ILLUMINACLIP="ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:keepBothReads"
LEADING=20
TRAILING=20
SLIDINGWINDOW="4:20"
MINLEN=75

# ------------------------------------------------------------------------------
# SECTION 5: ALIGNMENT SETTINGS (HISAT2)
# ------------------------------------------------------------------------------
HISAT2_OPTS="--dta"

# ------------------------------------------------------------------------------
# SECTION 6: SAMTOOLS SETTINGS
# ------------------------------------------------------------------------------
SORT_MEMORY="1500M"

# ------------------------------------------------------------------------------
# SECTION 7: VARIANT CALLING SETTINGS (BCFtools)
# ------------------------------------------------------------------------------
MIN_BASE_QUALITY=20
MIN_MAPPING_QUALITY=20
VCF_MIN_QUAL=20
VCF_MIN_DEPTH=10
VCF_FILTER_EXPR="QUAL>=${VCF_MIN_QUAL} && INFO/DP>=${VCF_MIN_DEPTH}"

# ------------------------------------------------------------------------------
# SECTION 8: STRINGTIE SETTINGS
# ------------------------------------------------------------------------------
STRINGTIE_MIN_COV=2.5
STRINGTIE_MIN_LEN=200
STRINGTIE_OPTS="-m ${STRINGTIE_MIN_LEN} -c ${STRINGTIE_MIN_COV}"
MERGED_GTF="${STRINGTIE_DIR}/merged_transcripts.gtf"

# ------------------------------------------------------------------------------
# SECTION 9: FEATURECOUNTS SETTINGS
# ------------------------------------------------------------------------------
FEATURECOUNTS_OPTS="-p --countReadPairs"
STRAND_SPECIFIC=0
COUNT_MATRIX_FILE="${COUNTS_DIR}/gene_counts.txt"
EXON_COUNT_MATRIX_FILE="${COUNTS_DIR}/exon_counts.txt"
EXON_FEATURE_TYPE="exon"

# ------------------------------------------------------------------------------
# SECTION 10: INPUT FILE PATTERNS
# ------------------------------------------------------------------------------
READ1_PATTERN="_L1_1"
READ2_PATTERN="_L1_2"

# ------------------------------------------------------------------------------
# SECTION 11: LOGGING
# ------------------------------------------------------------------------------
LOG_FILE="${LOGS_DIR}/upstream_pipeline_$(date +%Y%m%d_%H%M%S).log"

# ==============================================================================
# █▀▀ █▄ █ █▀▄   █▀█ █▀▀   █▀▀ █▀█ █▄ █ █▀▀ █ █▀▀ █ █ █▀█ ▄▀█ ▀█▀ █ █▀█ █▄ █
# ██▄ █ ▀█ █▄▀   █▄█ █▀    █▄▄ █▄█ █ ▀█ █▀  █ █▄█ █▄█ █▀▄ █▀█  █  █ █▄█ █ ▀█
# ==============================================================================

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

log_message() {
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] $1" | tee -a "${LOG_FILE}"
}

create_directories() {
    log_message "Creating output directories..."
    mkdir -p "${FASTQC_PRE_DIR}" "${TRIMMED_DIR}" "${FASTQC_POST_DIR}"
    mkdir -p "${SORTED_DIR}" "${VCF_DIR}"
    mkdir -p "${STRINGTIE_DIR}" "${COUNTS_DIR}" "${LOGS_DIR}" "${MULTIQC_DIR}"
    log_message "✓ Output directories created"
}

check_dependencies() {
    log_message "Checking dependencies..."
    local missing=()
    
    for cmd in fastqc trimmomatic hisat2 samtools bcftools stringtie gffcompare featureCounts multiqc; do
        if ! command -v "${cmd}" &> /dev/null; then
            missing+=("${cmd}")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log_message "ERROR: Missing tools: ${missing[*]}"
        exit 1
    fi
    log_message "✓ All dependencies available"
}

print_tool_versions() {
    log_message "=============================================="
    log_message "  SESSION INFO: Tool Versions"
    log_message "=============================================="
    log_message "  Date/Time: $(date '+%Y-%m-%d %H:%M:%S')"
    log_message "  Hostname:  $(hostname)"
    log_message "  User:      $(whoami)"
    log_message "----------------------------------------------"
    log_message "Tool Versions:"
    
    local fastqc_ver=$(fastqc --version 2>&1 | head -n1 || echo "N/A")
    log_message "  FastQC:       ${fastqc_ver}"
    
    local trimmo_ver=$(trimmomatic -version 2>&1 | head -n1 || echo "N/A")
    log_message "  Trimmomatic:  ${trimmo_ver}"
    
    local hisat2_ver=$(hisat2 --version 2>&1 | head -n1 | awk '{print $NF}' || echo "N/A")
    log_message "  HISAT2:       ${hisat2_ver}"
    
    local samtools_ver=$(samtools --version 2>&1 | head -n1 || echo "N/A")
    log_message "  Samtools:     ${samtools_ver}"
    
    local bcftools_ver=$(bcftools --version 2>&1 | head -n1 || echo "N/A")
    log_message "  BCFtools:     ${bcftools_ver}"
    
    local stringtie_ver=$(stringtie --version 2>&1 || echo "N/A")
    log_message "  StringTie:    ${stringtie_ver}"
    
    local fc_ver=$(featureCounts -v 2>&1 | head -n1 || echo "N/A")
    log_message "  FeatureCounts: ${fc_ver}"
    
    local multiqc_ver=$(multiqc --version 2>&1 | awk '{print $NF}' || echo "N/A")
    log_message "  MultiQC:      ${multiqc_ver}"
    
    log_message "=============================================="
}

check_reference_files() {
    log_message "Checking reference files..."
    
    if [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
        log_message "ERROR: HISAT2 index not found at ${HISAT2_INDEX}"
        exit 1
    fi
    
    if [[ ! -f "${REFERENCE_FASTA}" ]]; then
        log_message "ERROR: Reference FASTA not found at ${REFERENCE_FASTA}"
        exit 1
    fi
    
    if [[ ! -f "${GTF_FILE}" ]]; then
        log_message "ERROR: GTF file not found at ${GTF_FILE}"
        exit 1
    fi
    
    log_message "✓ Reference files found"
}

get_sample_name() {
    local filename
    filename=$(basename "$1")
    echo "${filename}" | sed -E 's/_L1_[12]\.(fq|fastq)(\.gz)?$//'
}

# Initialize logs directory
mkdir -p "${LOGS_DIR}"

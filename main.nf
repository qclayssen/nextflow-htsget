nextflow.enable.dsl=2

// Default parameters
params.outdir = params.outdir ?: 'results'
params.samplesheet = params.samplesheet ?: 'samplesheet.csv'
params.publish_dir_mode = params.publish_dir_mode ?: 'copy'

// Read the samplesheet, validate input, and prepare metadata for processing
def prepare_input() {
    // Supported file types with their metadata
    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]
    // Validate required parameters
    if (!params.samplesheet) {
        log.error "Parameter --samplesheet is required"
        System.exit(1)
    }

    if (!file(params.samplesheet).exists()) {
        log.error "Samplesheet file '${params.samplesheet}' does not exist"
        System.exit(1)
    }

    return Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sampleId = row.id?.trim()
            def sampleName = row.name?.trim()
            def sampleLabel = sampleId ?: sampleName
            if (!sampleLabel) {
                log.error "Row is missing a sample identifier (name or id)."
                System.exit(1)
            }

            def uri = row.uri?.trim()
            if (!uri) {
                log.error "Row for '${sampleLabel}' is missing a uri"
                System.exit(1)
            }

            def filetype = row.filetype?.trim()?.toLowerCase()
            if (!filetype || !FILETYPE_METADATA.containsKey(filetype)) {
                log.error "Row for '${sampleLabel}' has unsupported filetype '${row.filetype}'. Supported: ${FILETYPE_METADATA.keySet().join(', ')}"
                System.exit(1)
            }

            def referenceName = row.containsKey('reference_name') ? row.reference_name?.trim() : null
            def startRaw = row.containsKey('start') ? row.start?.trim() : null
            def endRaw = row.containsKey('end') ? row.end?.trim() : null

            Long startPos = null
            if (startRaw && startRaw != '') {
                try {
                    startPos = startRaw.toLong()
                    if (startPos < 0) {
                        log.error "Row for '${sampleLabel}' has invalid start position '${startRaw}' (must be >= 0)"
                        System.exit(1)
                    }
                } catch (NumberFormatException _e) {
                    log.error "Row for '${sampleLabel}' has invalid start '${startRaw}' (must be a number)"
                    System.exit(1)
                }
            }

            Long endPos = null
            if (endRaw && endRaw != '') {
                try {
                    endPos = endRaw.toLong()
                    if (endPos < 0) {
                        log.error "Row for '${sampleLabel}' has invalid end position '${endRaw}' (must be >= 0)"
                        System.exit(1)
                    }
                    if (startPos != null && endPos <= startPos) {
                        log.error "Row for '${sampleLabel}' has end position '${endRaw}' <= start position '${startRaw}'"
                        System.exit(1)
                    }
                } catch (NumberFormatException _e) {
                    log.error "Row for '${sampleLabel}' has invalid end '${endRaw}' (must be a number)"
                    System.exit(1)
                }
            }

            // Create metadata with attached file type information
            def meta = [
                id             : sampleLabel,
                sample_id      : sampleId ?: sampleLabel,
                sample_name    : sampleName ?: sampleLabel,
                filetype       : filetype,
                uri            : uri,
                reference_name : referenceName,
                start          : startPos,
                end            : endPos
            ] + FILETYPE_METADATA[filetype]

            return tuple(meta, meta.uri)
        }
}

process FETCH_FILE_PYTHON {
    tag { meta.id }
    publishDir "${params.outdir}/downloads/python", mode: params.publish_dir_mode

    input:
        tuple val(meta), val(uri)

    output:
        tuple val(meta), path("${meta.id}.${meta.extension}")

    script:
    def args = [
        "--url '${uri}'",
        "--output '${meta.id}.${meta.extension}'"
    ]
    if (meta.format)          args << "--format '${meta.format}'"
    if (meta.reference_name)  args << "--reference-name '${meta.reference_name}'"
    if (meta.start != null)   args << "--start ${meta.start}"
    if (meta.end != null)     args << "--end ${meta.end}"
    """
    set -euo pipefail
    htsget_fetch.py ${args.join(' ')}
    """
}

process FETCH_FILE_CLI {
    tag { meta.id }
    publishDir "${params.outdir}/downloads/cli", mode: params.publish_dir_mode

    input:
        tuple val(meta), val(uri)

    output:
        tuple val(meta), path("${meta.id}.cli.${meta.extension}")

    script:
    def args = ['htsget']
    if (meta.format)          args << "--format=${meta.format}"
    if (meta.reference_name)  args << "--reference-name=${meta.reference_name}"
    if (meta.start != null)   args << "--start=${meta.start}"
    if (meta.end != null)     args << "--end=${meta.end}"
    args << "--output=${meta.id}.cli.${meta.extension}"
    args << "'${uri}'"
    """
    set -euo pipefail
    
    # Check if htsget CLI tool is available
    if ! command -v htsget &> /dev/null; then
        echo "Error: htsget CLI tool not found. Please install it first." >&2
        exit 1
    fi
    
    # Execute the htsget command
    ${args.join(' ')}
    
    # Validate output file was created
    if [[ ! -f "${meta.id}.cli.${meta.extension}" ]]; then
        echo "Error: Output file '${meta.id}.cli.${meta.extension}' was not created" >&2
        exit 1
    fi
    """
}

process RUN_FASTQC {
    tag { meta.id }
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"),  emit: fastqc_reports
    tuple val(meta), path("*_fastqc.html"), emit: fastqc_html

    when:
    meta.filetype == 'fastq'

    script:
    """
    set -euo pipefail
    
    # Check if fastqc is available
    if ! command -v fastqc &> /dev/null; then
        echo "Error: FastQC not found. Please install it first." >&2
        exit 1
    fi
    
    fastqc --quiet --outdir . ${reads}
    """
}

process RUN_SAMTOOLS_STATS {
    tag { meta.id }
    publishDir "${params.outdir}/qc/samtools", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.samtools_stats.txt"),   emit: samtools_stats
    tuple val(meta), path("${meta.id}.samtools_flagstat.txt"), emit: samtools_flagstat

    when:
    meta.filetype == 'bam'

    script:
    """
    set -euo pipefail
    
    # Check if samtools is available
    if ! command -v samtools &> /dev/null; then
        echo "Error: Samtools not found. Please install it first." >&2
        exit 1
    fi
    
    samtools stats ${bam} > ${meta.id}.samtools_stats.txt
    samtools flagstat ${bam} > ${meta.id}.samtools_flagstat.txt
    """
}

process RUN_BCFTOOLS_STATS {
    tag { meta.id }
    publishDir "${params.outdir}/qc/bcftools", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.bcftools_stats.txt"), emit: bcftools_stats

    when:
    meta.filetype == 'vcf'

    script:
    """
    set -euo pipefail
    
    # Check if bcftools is available
    if ! command -v bcftools &> /dev/null; then
        echo "Error: BCFtools not found. Please install it first." >&2
        exit 1
    fi
    
    bcftools stats ${vcf} > ${meta.id}.bcftools_stats.txt
    """
}

process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/qc/multiqc", mode: params.publish_dir_mode

    input:
    path(qc_files)

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    when:
    qc_files

    script:
    """
    set -euo pipefail
    
    # Check if multiqc is available
    if ! command -v multiqc &> /dev/null; then
        echo "Error: MultiQC not found. Please install it first." >&2
        exit 1
    fi
    
    # Run MultiQC
    multiqc . -o . --force
    
    # Validate output was created
    if [[ ! -f "multiqc_report.html" ]]; then
        echo "Warning: MultiQC report was not generated" >&2
    fi
    """
}

workflow {
    meta_uri_ch = prepare_input()

    fetched_python_ch = FETCH_FILE_PYTHON(meta_uri_ch)
    // CLI fetching is available but not used in the main workflow
    // fetched_cli_ch = FETCH_FILE_CLI(meta_uri_ch)

    fastqc_out    = RUN_FASTQC(fetched_python_ch)
    samtools_out  = RUN_SAMTOOLS_STATS(fetched_python_ch)
    bcftools_out  = RUN_BCFTOOLS_STATS(fetched_python_ch)

    // Collect QC files for MultiQC
    multiqc_inputs = Channel.empty()
        .mix(
            fastqc_out.fastqc_reports.map { _meta, file -> file },
            samtools_out.samtools_stats.map { _meta, file -> file },
            samtools_out.samtools_flagstat.map { _meta, file -> file },
            bcftools_out.bcftools_stats.map { _meta, file -> file }
        )
        .collect()
    
    MULTIQC(multiqc_inputs)
}

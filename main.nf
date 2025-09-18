nextflow.enable.dsl=2

params.outdir = params.outdir ?: 'results'
params.samplesheet = params.samplesheet ?: 'samplesheet.csv'
params.publish_dir_mode = params.publish_dir_mode ?: 'copy'

// Helper that returns the rows from samplesheet.
def sample_rows() {
    // Supported file types with their corresponding output extensions and HTSGET formats.
    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.name?.trim() ?: row.id?.trim()
            if (!sample) {
                log.error "Row is missing a sample identifier (name or id)."
                System.exit(1)
            }

            def uri = row.uri?.trim()
            if (!uri) {
                log.error "Row for '${sample}' is missing a uri"
                System.exit(1)
            }

            def filetype = row.filetype?.trim()?.toLowerCase()
            if (!filetype || !FILETYPE_METADATA.containsKey(filetype as String)) {
                log.error "Row for '${sample}' has unsupported filetype '${row.filetype}'. Supported: ${FILETYPE_METADATA.keySet().join(', ')}"
                System.exit(1)
            }

            tuple(sample, filetype as String, uri)
        }
}

process FETCH_FILE {
    tag { sample }
    publishDir "${params.outdir}/downloads", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(filetype), val(extension), val(format), val(uri)

    output:
        tuple val(sample), val(filetype), path("${sample}.${extension}")

    script:
    """
    set -euo pipefail
    fetch_htsget.py '${uri}' '${sample}.${extension}'
    """
}

process RUN_FASTQC {
    tag { sample }
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(filetype), path(reads)

    output:
        path "${sample}_fastqc.zip", emit: fastqc_reports
        path "${sample}_fastqc.html"

    script:
    """
    set -euo pipefail
    fastqc --quiet --outdir . ${reads}
    """
}

process RUN_SAMTOOLS_STATS {
    tag { sample }
    publishDir "${params.outdir}/qc/samtools", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(filetype), path(bam)

    output:
        path "${sample}.samtools_stats.txt", emit: samtools_stats
        path "${sample}.samtools_flagstat.txt", emit: samtools_flagstat

    script:
    """
    set -euo pipefail
    samtools stats ${bam} > ${sample}.samtools_stats.txt
    samtools flagstat ${bam} > ${sample}.samtools_flagstat.txt
    """
}

process RUN_BCFTOOLS_STATS {
    tag { sample }
    publishDir "${params.outdir}/qc/bcftools", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(filetype), path(vcf)

    output:
        path "${sample}.bcftools_stats.txt", emit: bcftools_stats

    script:
    """
    set -euo pipefail
    bcftools stats ${vcf} > ${sample}.bcftools_stats.txt
    """
}

process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/qc/multiqc", mode: params.publish_dir_mode

    input:
        path qc_files

    when:
        qc_files && qc_files.size() > 0

    output:
        path "multiqc_report.html"
        path "multiqc_data"

    script:
    """
    set -euo pipefail
    multiqc . -o .
    """
}

workflow {
    // Supported file types with their corresponding output extensions and HTSGET formats.
    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]

    sample_ch = sample_rows()

    prepared_ch = sample_ch.map { sample, filetype, uri ->
        def meta = FILETYPE_METADATA[filetype]
        tuple(sample, filetype, meta.extension, meta.format, uri)
    }

    fetched_ch = FETCH_FILE(prepared_ch)

    fastq_ch = fetched_ch.filter { sample, filetype, _file -> filetype == 'fastq' }
    fastqc_out = RUN_FASTQC(fastq_ch)

    bam_ch = fetched_ch.filter { sample, filetype, _file -> filetype == 'bam' }
    samtools_out = RUN_SAMTOOLS_STATS(bam_ch)

    vcf_ch = fetched_ch.filter { sample, filetype, _file -> filetype == 'vcf' }
    bcftools_out = RUN_BCFTOOLS_STATS(vcf_ch)

    // Collect QC files conditionally
    multiqc_inputs = Channel.empty()

    // Add FASTQC reports if any FASTQ files were processed
    fastqc_out.fastqc_reports
        .ifEmpty([])
        .set { fastqc_files }
    multiqc_inputs = multiqc_inputs.mix(fastqc_files)

    // Add SAMTOOLS reports if any BAM files were processed
    samtools_out.samtools_stats
        .ifEmpty([])
        .set { samtools_stats_files }
    multiqc_inputs = multiqc_inputs.mix(samtools_stats_files)

    samtools_out.samtools_flagstat
        .ifEmpty([])
        .set { samtools_flagstat_files }
    multiqc_inputs = multiqc_inputs.mix(samtools_flagstat_files)

    // Add BCFTOOLS reports if any VCF files were processed
    bcftools_out.bcftools_stats
        .ifEmpty([])
        .set { bcftools_files }
    multiqc_inputs = multiqc_inputs.mix(bcftools_files)

    multiqc_inputs = multiqc_inputs.flatten().collect()

    MULTIQC(multiqc_inputs)
}

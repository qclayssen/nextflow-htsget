nextflow.enable.dsl=2

// Default parameters
params.outdir = params.outdir ?: 'results'
params.samplesheet = params.samplesheet ?: 'samplesheet.csv'
params.publish_dir_mode = params.publish_dir_mode ?: 'copy'

include { PREPARE_INPUT } from '../subworkflows/prepare_input'

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

    python htsget_fetch.py ${args.join(' ')}

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
    if (meta.format)          args << "--format=${meta.format}"
    if (meta.reference_name)  args << "--reference-name=${meta.reference_name}"
    if (meta.start != null)   args << "--start=${meta.start}"
    if (meta.end != null)     args << "--end=${meta.end}"

    """

    htsget "'${uri}'" "--output=${meta.id}.cli.${meta.extension}" ${args.join(' ')}

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

    multiqc . -o . --force

    """
}

workflow {
    meta_uri_ch = PREPARE_INPUT()

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

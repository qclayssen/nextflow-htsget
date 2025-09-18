process FETCH_WITH_CLI {
    tag "$sample"
    label 'process_single'

    publishDir "${params.outdir}/cli", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(uri)

    output:
    path "${sample}_cli.bam", emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    htsget get ${uri} -O ${sample}_cli.bam ${args}
    """

    stub:
    """
    touch ${sample}_cli.bam
    """
}
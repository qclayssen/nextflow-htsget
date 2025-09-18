process FETCH_WITH_CURL {
    tag "$sample"
    label 'process_single'

    publishDir "${params.outdir}/curl", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(uri)

    output:
    path "${sample}_discovery.json", emit: json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Convert htsget:// to https:// for curl
    def curl_uri = uri.replace('htsget://', 'https://')
    """
    curl -sS '${curl_uri}' -o '${sample}_discovery.json' ${args}
    """

    stub:
    """
    echo '{"format": "BAM", "urls": []}' > ${sample}_discovery.json
    """
}
process FETCH_WITH_PYTHON {
    tag "$sample"
    label 'process_single'

    publishDir "${params.outdir}/python", mode: params.publish_dir_mode

    input:
    tuple val(sample), val(uri)

    output:
    path "${sample}_python.bam", emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python3 -c "
    import htsget
    import sys
    url = '${uri}'
    try:
        with open('${sample}_python.bam', 'wb') as output:
            htsget.get(url, output)
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)
    "
    """

    stub:
    """
    touch ${sample}_python.bam
    """
}
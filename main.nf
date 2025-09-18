nextflow.enable.dsl=2

params.outdir = params.outdir ?: 'results'
params.samplesheet = params.samplesheet ?: 'samplesheet.csv'
params.publish_dir_mode = params.publish_dir_mode ?: 'copy'

// Helper that returns the rows from samplesheet.
def sample_rows() {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def name = row.name?.trim() ?: row.id?.trim()
            def uri = row.uri?.trim()
            if (!uri) {
                log.error "Row for '${name}' is missing a uri"
                System.exit(1)
            }
            tuple(name, uri)
        }
}

process FETCH_WITH_CLI {
    tag { sample }
    publishDir "${params.outdir}/cli", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(uri)

    output:
        path "${sample}_cli.bam"

    script:
    """
    htsget ${uri} -O ${sample}_cli.bam
    """
}

process FETCH_WITH_PYTHON {
    tag { sample }
    publishDir "${params.outdir}/python", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(uri)

    output:
        path "${sample}_python.bam"

    script:
    """
    set -euo pipefail
    python3 -c "
    import htsget
    url = '${uri}'
    with open('${sample}_python.bam', 'wb') as output:
        htsget.get(url, output)
    "
    """
}

process FETCH_WITH_CURL {
    tag { sample }
    publishDir "${params.outdir}/curl", mode: params.publish_dir_mode

    input:
        tuple val(sample), val(uri)

    output:
        path "${sample}_discovery.json"

    script:
    """
    set -euo pipefail
    curl -sS '${uri.replace('htsget://', 'https://')}' -o '${sample}_discovery.json'
    """
}

workflow {
    sample_ch = sample_rows()
    FETCH_WITH_CLI(sample_ch)
    FETCH_WITH_PYTHON(sample_ch)
    FETCH_WITH_CURL(sample_ch)
}

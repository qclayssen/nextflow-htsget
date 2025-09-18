nextflow.enable.dsl=2

params.outdir = params.outdir ?: 'results'
params.samplesheet = params.samplesheet ?: 'samplesheet.csv'

// Helper that returns the rows for a given method.
def method_rows(String method) {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .filter { row ->
            row.filetype?.trim()?.toLowerCase() == method
        }
        .map { row ->
            def name = row.name?.trim() ?: row.id?.trim() ?: method
            def uri = row.uri?.trim()
            if (!uri) {
                log.error "Row for '${name}' (${method}) is missing a uri"
                System.exit(1)
            }
            def filename = row.filename?.trim()
            if (!filename) {
                filename = "${name}_${method}"
            }
            tuple(name, uri, filename)
        }
}

process FETCH_WITH_CLI {
    tag { sample }
    publishDir "${params.outdir}/cli", mode: 'copy'

    input:
        tuple val(sample), val(uri), val(filename)

    output:
        path filename

    script:
    """
    set -euo pipefail
    htsget get '${uri}' --output '${filename}'
    """
}

process FETCH_WITH_PYTHON {
    tag { sample }
    publishDir "${params.outdir}/python", mode: 'copy'

    input:
        tuple val(sample), val(uri), val(filename)

    output:
        path filename

    script:
    """
    set -euo pipefail
    bin/fetch_htsget.py '${uri}' '${filename}'
    """
}

process FETCH_WITH_CURL {
    tag { sample }
    publishDir "${params.outdir}/curl", mode: 'copy'

    input:
        tuple val(sample), val(uri), val(filename)

    output:
        path filename

    script:
    """
    set -euo pipefail
    curl -sS '${uri}' -o '${filename}'
    """
}

workflow {
    FETCH_WITH_CLI(method_rows('cli'))
    FETCH_WITH_PYTHON(method_rows('python'))
    FETCH_WITH_CURL(method_rows('curl'))
}

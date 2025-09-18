/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    qclayssen/nextflow-htsget WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_WITH_CLI    } from '../modules/local/fetch_cli'
include { FETCH_WITH_PYTHON } from '../modules/local/fetch_python'  
include { FETCH_WITH_CURL   } from '../modules/local/fetch_curl'

workflow NEXTFLOWHTSGET {

    //
    // Read samplesheet and create input channel
    //
    ch_samplesheet = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.id?.trim()
            def sample_name = row.name?.trim() ?: sample_id
            def uri = row.uri?.trim()
            
            if (!sample_id) {
                log.error "Missing sample ID in samplesheet row: ${row}"
                System.exit(1)
            }
            if (!uri) {
                log.error "Missing URI for sample '${sample_id}'"
                System.exit(1)
            }
            
            tuple(sample_name, uri)
        }

    //
    // Fetch data using different methods
    //
    FETCH_WITH_CLI    ( ch_samplesheet )
    FETCH_WITH_PYTHON ( ch_samplesheet )
    FETCH_WITH_CURL   ( ch_samplesheet )

}
// Subworkflow to prepare input channel from samplesheet
workflow PREPARE_INPUT {
    params.outdir = params.outdir ?: 'results'
    params.samplesheet = params.samplesheet ?: 'samplesheet.csv'

    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]

    ch_meta_uri = Channel
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

            if ((startPos != null || endPos != null) && !referenceName) {
                log.error "Row for '${sampleLabel}' specifies start/end positions but is missing reference_name. Reference name is required when using genomic ranges."
                System.exit(1)
            }

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

    emit:
    ch_meta_uri
}
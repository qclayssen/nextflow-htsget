// Subworkflow to prepare input channel from samplesheet or direct URL list
workflow PREPARE_INPUT {
    main:
    params.outdir        = params.outdir ?: 'results'
    params.samplesheet   = params.samplesheet ?: 'samplesheet.csv'

    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]

    // Normalise manual URL inputs from single-field UI or CLI list flags
    def manualUrls = []
    def addUrls = { value ->
        if (!value) {
            return
        }
        if (value instanceof Collection) {
            value.each { addUrls(it) }
            return
        }

        def text = value.toString()
        if (!text) {
            return
        }

        text
            .split(/[\n,]/)
            .collect { it?.trim() }
            .findAll { it }
            .each { manualUrls << it }
    }

    addUrls(params.htsget_urls)
    addUrls(params.htsget_url)

    manualUrls = manualUrls
        .collect { it?.toString()?.trim() }
        .findAll { it }

    def configuredFiletype = params.htsget_filetype?.toString()?.trim()

    def manualTuples = []
    manualUrls.eachWithIndex { url, idx ->
        def suffix        = idx + 1
        def baseComponent = url.tokenize('/')?.last()?.split('\\?')[0]
        def safeLabel     = baseComponent?.replaceAll(/[^A-Za-z0-9._\-]/, '_')?.take(40)
        def sampleLabel   = safeLabel ?: "url_${suffix}"

        def detectedFiletype = configuredFiletype?.toLowerCase()
        if (!detectedFiletype) {
            // Auto-detect filetype from URL path when not specified
            if (url.contains('/variants/')) {
                detectedFiletype = 'vcf'
            } else if (url.contains('/reads/')) {
                detectedFiletype = 'bam'
            } else {
                detectedFiletype = 'bam'  // default fallback
            }
        }
        if (!FILETYPE_METADATA.containsKey(detectedFiletype)) {
            def invalidValue = configuredFiletype ?: params.htsget_filetype
            log.error "Invalid htsget_filetype '${invalidValue}'. Supported: ${FILETYPE_METADATA.keySet().join(', ')}"
            System.exit(1)
        }

        def meta = [
            id             : sampleLabel,
            sample_id      : sampleLabel,
            sample_name    : sampleLabel,
            filetype       : detectedFiletype,
            uri            : url,
            reference_name : null,
            start          : null,
            end            : null
        ] + FILETYPE_METADATA[detectedFiletype]

        manualTuples << tuple(meta, meta.uri)
    }

    def manualChannel = manualTuples ? Channel.of(*manualTuples) : Channel.empty()

    def sheetPathStr = params.samplesheet?.toString()?.trim()
    if (sheetPathStr == '') {
        sheetPathStr = null
    }
    def sheetPath    = sheetPathStr ? file(sheetPathStr, checkIfExists: false) : null
    def sheetExists  = sheetPath?.exists()

    def rowToTuple = { row ->
        def sampleId    = row.id?.trim()
        def sampleName  = row.name?.trim()
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
        def startRaw      = row.containsKey('start') ? row.start?.trim() : null
        def endRaw        = row.containsKey('end') ? row.end?.trim() : null

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

        tuple(meta, meta.uri)
    }

    if (sheetExists && !manualTuples.isEmpty()) {
        log.error "Both samplesheet '${sheetPathStr}' and manual HTSGET URLs were provided. Choose only one input method."
        System.exit(1)
    }

    def sheetChannel = Channel.empty()
    if (sheetExists) {
        def delimiter = sheetPathStr.toLowerCase().endsWith('.tsv') ? "\t" : ','
        sheetChannel = Channel
            .fromPath(sheetPath)
            .splitCsv(header: true, sep: delimiter)
            .map(rowToTuple)
    }

    if (!sheetExists && manualTuples.isEmpty()) {
        log.error "No records to process: provide a samplesheet or at least one manual HTSGET URL."
        System.exit(1)
    }

    ch_meta_uri = sheetExists ? sheetChannel : manualChannel

    emit:
    ch_meta_uri
}

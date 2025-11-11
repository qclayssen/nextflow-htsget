// Normalise any mix of collections or comma/newline separated strings into a flat list of URLs
def normaliseUrlInput(value) {
    if (!value) {
        return []
    }

    def candidates = value instanceof Collection ? value : [value]

    candidates
        .collect { it?.toString() }
        .findAll { it }
        .collectMany { item ->
            item
                .split(/[\n,]/)
                .collect { it?.trim() }
                .findAll { it }
        }
}

// Subworkflow to prepare input channel from samplesheet or direct URL list
workflow PREPARE_INPUT {
    main:
    // assume `params.outdir` validated at workflow entry (handled in main.nf)

    def FILETYPE_METADATA = [
        bam  : [extension: 'bam',  format: 'BAM'],
        fastq: [extension: 'fastq', format: 'FASTQ'],
        vcf  : [extension: 'vcf',  format: 'VCF']
    ]

    // Normalise manual URL inputs from form or CLI overrides
    def manualUrls = normaliseUrlInput(params.htsget_urls)

    def configuredFiletype = params.htsget_filetype?.toString()?.trim()

    def manualTuples = []
    def usedLabels = [:] // Track used labels to ensure uniqueness
    manualUrls.eachWithIndex { url, idx ->
        def suffix        = idx + 1
        def baseComponent = url.tokenize('/')?.last()?.split('\\?')[0]
        def safeLabel     = baseComponent?.replaceAll(/[^A-Za-z0-9._\-]/, '_')?.take(40)
        def sampleLabel   = safeLabel ?: "url_${suffix}"

        // Ensure unique sample labels
        if (usedLabels.containsKey(sampleLabel)) {
            def counter = usedLabels[sampleLabel] + 1
            usedLabels[sampleLabel] = counter
            sampleLabel = "${sampleLabel}_${counter}"
        } else {
            usedLabels[sampleLabel] = 1
        }

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

    def manualChannel = manualTuples ? Channel.from(manualTuples) : Channel.empty()

    def sheetPathStr = params.samplesheet?.toString()?.trim()
    if (sheetPathStr == '') {
        sheetPathStr = null
    }
    def sheetProvided = sheetPathStr != null

    def rowToTuple = { row ->
        def caseId    = row.CASEID?.trim()
        def patientId  = row.PATIENTID?.trim()
        def specimenId  = row.SPECIMENID?.trim()
        def sampleLabel = caseId + "_" + patientId + "_" + specimenId
        if (!sampleLabel) {
            log.error "Row is missing a unique identifier (must have at least one of CASEID, PATIENTID, SPECIMENID)."
            System.exit(1)
        }

        def uri = row.OBJECTSTOREURL?.trim()
        if (!uri) {
            log.error "Row for '${sampleLabel}' is missing a OBJECTSTOREURL"
            System.exit(1)
        }

        def filetype = row.OBJECTTYPE?.trim()?.toLowerCase()
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
            sample_id      : specimenId ?: sampleLabel,
            sample_name    : specimenId ?: sampleLabel,
            filetype       : filetype,
            uri            : uri,
            reference_name : referenceName,
            start          : startPos,
            end            : endPos
        ] + FILETYPE_METADATA[filetype]

        tuple(meta, meta.uri)
    }

    if (sheetProvided && !manualTuples.isEmpty()) {
        log.error "Both samplesheet '${sheetPathStr}' and manual HTSGET URLs were provided. Choose only one input method."
        System.exit(1)
    }

    def sheetChannel = Channel.empty()
    if (sheetProvided) {
        def sheetPath = file(sheetPathStr, checkIfExists: false)
        def delimiter = sheetPathStr.toLowerCase().endsWith('.tsv') ? "\t" : ','
        sheetChannel = Channel
            .fromPath(sheetPath)
            .splitCsv(header: true, sep: delimiter)
            .map(rowToTuple)
    }

    if (!sheetProvided && manualTuples.isEmpty()) {
        log.error "No records to process: provide a samplesheet or at least one manual HTSGET URL."
        System.exit(1)
    }

    ch_meta_uri = sheetProvided ? sheetChannel : manualChannel

    emit:
    ch_meta_uri
}

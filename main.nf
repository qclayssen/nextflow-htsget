#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    qclayssen/nextflow-htsget
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GitHub : https://github.com/qclayssen/nextflow-htsget
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Print help message if needed
if (params.help) {
   log.info """
   ========================================================================================
                      q c l a y s s e n / n e x t f l o w - h t s g e t
   ========================================================================================
   Demonstration pipeline that downloads GA4GH HTSGET data via three methods.
   
   Usage:
   nextflow run qclayssen/nextflow-htsget --samplesheet samplesheet.csv --outdir results -profile docker
   
   Mandatory arguments:
     --samplesheet         Path to CSV samplesheet containing sample information
     --outdir              Output directory for results
   
   Options:
     -profile              Configuration profile (docker, singularity, etc.)
     --help                Show this help message and exit
     --version             Show version and exit
   
   ========================================================================================
   """.stripIndent()
   exit 0
}

// Print version and exit
if (params.version) {
   log.info "qclayssen/nextflow-htsget version: ${workflow.manifest.version}"
   exit 0
}

// Check mandatory parameters
if (!params.samplesheet) {
   log.error "Please provide a samplesheet with --samplesheet"
   exit 1
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INCLUDE WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NEXTFLOWHTSGET } from './workflows/nextflowhtsget'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    NEXTFLOWHTSGET()
}

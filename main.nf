#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GdTBioinfo-nf/snvs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/GdTBioinfo-nf/snvs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ojo, hay que instalar el modulo para que funcione : nf-core modules install pangolin
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GRCh37 -profile docker"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log, args)

// Include the SNVS workflow
include { SNVS } from './workflows/snvs'

// Include the pangolin module
include { PANGOLIN } from 'nf-core/modules/pangolin'

// Define your input for pangolin (Este moduloparte del fastq, pero podriamos modificarlo si deseais otra cosa)
params.input_file = file('ruta/al/archivo/input.fasta')

// Define el workflow para pangolin
workflow PANGOLIN_WORKFLOW {
    PANGOLIN(input_file: params.input_file)
}

// MAIN WORKFLOW
workflow GDTBIOINFONF_SNVS {
    SNVS()
    PANGOLIN_WORKFLOW()
}

// Execute all workflows
workflow {
    GDTBIOINFONF_SNVS()
}


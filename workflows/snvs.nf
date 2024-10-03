/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowSnvs.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters

if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta file not specified!' }
if (params.fai) { ch_fai = file(params.fai) } else { exit 1, 'Fai file not specified!' }
if (params.dict) { ch_dict = file(params.dict) } else { exit 1, 'Dict file not specified!' }
if (params.index) { ch_index = file(params.index) } else { exit 1, 'Index file not specified!' }
if (params.bed) { ch_bed = Channel.fromPath(params.bed) } else { ch_bed = [] }
if (params.dgn_model) { ch_dgn_model = Channel.fromPath(params.dgn_model) } else { ch_dgn_model = [] }
if (params.dbsnp) { ch_dbsnp = Channel.fromPath(params.dbsnp) } else { ch_dbsnp = [] }
if (params.dbsnp_tbi) { ch_dbsnp_tbi = Channel.fromPath(params.dbsnp_tbi) } else { ch_dbsnp_tbi = [] }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_ALIGN_BWA } from '../subworkflows/nf-core/fastq_align_bwa/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { BWA_MEM               } from '../modules/nf-core/bwa/mem/main'
include { GATK4_FASTQTOSAM  } from '../modules/nf-core/gatk4/fastqtosam/main'
include { SAMTOOLS_SORT     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT  as  SAMTOOLS_SORT_MAPPED  } from '../modules/nf-core/samtools/sort/main'
include { GATK4_MERGEBAMALIGNMENT  } from '../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4_ADDORREPLACEREADGROUPS } from '../modules/nf-core/gatk4/addorreplacereadgroups/main' 

include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4SPARK_MARKDUPLICATES  } from '../modules/nf-core/gatk4spark/markduplicates/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../modules/nf-core/picard/addorreplacereadgroups/main'                                     


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //FASTQ_ALIGN_BWA (
    //    INPUT_CHECK.out.reads, tuple([], ch_index), false, tuple([], ch_fasta)
    //)

    BWA_MEM(
        INPUT_CHECK.out.reads,
        tuple([], ch_index),
        tuple([], ch_fasta),
        true
    )

    PICARD_ADDORREPLACEREADGROUPS (
        BWA_MEM.out.bam,
        tuple([], ch_fasta),
        tuple([], ch_fai)
    )

    
    GATK4_FASTQTOSAM (
        INPUT_CHECK.out.reads
    )

    SAMTOOLS_SORT (
        GATK4_FASTQTOSAM.out.bam,
        tuple([], ch_fasta)

    )

    SAMTOOLS_SORT_MAPPED (
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        tuple([], ch_fasta)

    )

    //SAMTOOLS_SORT_MAPPED.out.bam.join(SAMTOOLS_SORT.out.bam).view()
    //GATK4_FASTQTOSAM.out.bam.view()

    //GATK4_MERGEBAMALIGNMENT (
//
  //      SAMTOOLS_SORT_MAPPED.out.bam.join(SAMTOOLS_SORT.out.bam),
    //    tuple([], ch_fasta), 
    //    tuple([], ch_dict)
    //)

    //GATK4_MERGEBAMALIGNMENT.out.bam.view()

    GATK4SPARK_MARKDUPLICATES (
        PICARD_ADDORREPLACEREADGROUPS.out.bam,
        ch_fasta, 
        ch_fai,
        ch_dict
    )

    //GATK4SPARK_MARKDUPLICATES.out.output.view()

    //GATK4SPARK_BASERECALIBRATOR ()

    //ch_bed.view()
    //ch_bambai = FASTQ_ALIGN_BWA.out.bam.join(FASTQ_ALIGN_BWA.out.bai).concat(ch_bed, ch_dgn_model ).collect()//.view()
    if (params.bed) {
        ch_bambai = PICARD_ADDORREPLACEREADGROUPS.out.bam.combine(PICARD_ADDORREPLACEREADGROUPS.out.bai).combine(ch_bed).map{ meta, bam, meta2, bai, intervals -> [ meta, bam, bai, intervals, [] ] }
        } else {
        ch_bambai = PICARD_ADDORREPLACEREADGROUPS.out.bam.combine(PICARD_ADDORREPLACEREADGROUPS.out.bai).map{ meta, bam, meta2, bai -> [ meta, bam, bai, [], [] ] }
    }
    //PICARD_ADDORREPLACEREADGROUPS.out.bai.view()
    //ch_bambai.view()


    GATK4_HAPLOTYPECALLER (
        ch_bambai, 
        tuple([], ch_fasta), 
        tuple([], ch_fai), 
        tuple([], ch_dict), 
        tuple([], ch_dbsnp), 
        tuple([], ch_dbsnp_tbi)
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSnvs.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSnvs.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

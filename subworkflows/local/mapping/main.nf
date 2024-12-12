include { BWA_MEM                 } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_FASTQTOSAM } from '../../../modules/nf-core/gatk4/fastqtosam/main'
include { GATK4_MERGEBAMALIGNMENT } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4SPARK_MARKDUPLICATES } from '../../../modules/nf-core/gatk4spark/markduplicates/main'
include { PICARD_SORTSAM } from '../../../modules/nf-core/picard/sortsam/main'
// faltaría SetNmMdAndUqTags que no está en nf-core
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR } from '../../../modules/nf-core/gatk4/applybqsr/main'


workflow MAPPING {

    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    //val_sort_bam    // boolean (mandatory): true or false
    ch_fasta        // channel (optional) : [ val(meta3), path(fasta) ]
    ch_fai
    ch_refdict
    ch_intervals
    ch_known_sites
    ch_known_sites_tbi

    main:

    BWA_MEM (
        ch_reads,
        ch_index,
        ch_fasta,
        true
    )

    GATK4_FASTQTOSAM (
        ch_reads
    )
    
    GATK4_MERGEBAMALIGNMENT (
        BWA_MEM.out.bam.join(GATK4_FASTQTOSAM.out.bam),
        ch_fasta,
        ch_refdict
    )

    GATK4_MERGEBAMALIGNMENT.out.bam.view()


    GATK4SPARK_MARKDUPLICATES (
        GATK4_MERGEBAMALIGNMENT.out.bam,
        ch_fasta.map {meta , fasta -> [fasta]},
        ch_fai.map {meta, fai -> [fai]}, 
        ch_refdict.map {meta, dict -> [dict]}
    )


    // ch_known_sites.view()

    //GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index.combine(ch_intervals)).view()
    //GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index.join(ch_intervals)).view()

    GATK4_BASERECALIBRATOR (   
         GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index.join(ch_intervals)),
         ch_fasta.map {meta , fasta -> [fasta]},
         ch_fai.map {meta, fai -> [fai]}, 
         ch_refdict.map {meta, dict -> [dict]},
         ch_known_sites,
         ch_known_sites_tbi
   )

   GATK4_APPLYBQSR (
        GATK4SPARK_MARKDUPLICATES.out.output.join(GATK4SPARK_MARKDUPLICATES.out.bam_index).join(GATK4_BASERECALIBRATOR.out.table).join(ch_intervals),
        ch_fasta.map {meta , fasta -> [fasta]},
        ch_fai.map {meta, fai -> [fai]}, 
        ch_refdict.map {meta, dict -> [dict]}
   )
    
    emit:
    bam = GATK4_APPLYBQSR.out.bam // channel: [ val(meta), path(bam)]

}
// Run GATK4 Base Recalibration

include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR_1ST_PASS } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR_2ND_PASS } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR         as APPLYBQSR                 } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_ANALYZECOVARIATES as ANALYZECOVARIATES         } from '../../../modules/local/gatk4/analyzecovariates/main'

workflow BAM_BASERECALIBRATOR_APPLYBQSR_GATK {
    take:
    bam                // [mandatory] channel: tuple(meta, bam)
    bai                // [mandatory] channel: tuple(meta, bai)
    intervals          // [mandatory] channel: path(intervals)
    fasta              // [mandatory] channel: tuple(meta, fa)
    fasta_fai          // [mandatory] channel: tuple(meta, fai)
    dict               // [mandatory] channel: tuple(meta, dict)
    known_indels       // [mandatory] channel: tuple(meta, known_indels)
    known_indels_tbi   // [mandatory] channel: tuple(meta, known_indels_tbi)
    known_snps         // [mandatory] channel: tuple(meta, known_snps)
    known_snps_tbi     // [mandatory] channel: tuple(meta, known_snps_tbi)
    get_metrics        // [optional]  value: boolean

    main:
    versions   = Channel.empty()
    bam_recall = Channel.empty()
    bai_recall = Channel.empty()
    plots      = Channel.empty()

    known_sites = known_indels
        .mix(known_snps)
        .flatMap { meta, vcfs ->
            if (vcfs instanceof List) {
                vcfs.collect { vcf ->
                    tuple([id:"known_sites"], vcf)
                }
            } else {
                [ tuple([id:"known_sites"], vcfs) ]
            }
        }
        .groupTuple()

    known_sites_tbi = known_indels_tbi
        .mix(known_snps_tbi)
        .flatMap { meta, tbis ->
            if (tbis instanceof List) {
                tbis.collect { tbi ->
                    tuple([id:"known_sites"], tbi)
                }
            } else {
                [ tuple([id:"known_sites"], tbis) ]
            }
        }
        .groupTuple()

    // First pass of Base Quality Score Recalibration (BQSR)
    BASERECALIBRATOR_1ST_PASS(
        bam
            .join(bai)
            .combine(intervals)
            .map { meta, bam, bai, bed ->
                tuple(meta, bam, bai, bed)
            },
        fasta,
        fasta_fai,
        dict,
        known_sites,
        known_sites_tbi
    )

    // Apply Base Quality Score Recalibration (BQSR)
    APPLYBQSR(
        bam
            .join(bai)
            .join(BASERECALIBRATOR_1ST_PASS.out.table)
            .combine(intervals)
            .map { meta, bam, bai, table, bed ->
                tuple(meta, bam, bai, table, bed)
            },
        fasta.map     { it[1] },
        fasta_fai.map { it[1] },
        dict.map      { it[1] }
    )

    bam_recall = APPLYBQSR.out.bam
    bai_recall = APPLYBQSR.out.bai
    versions   = versions.mix(BASERECALIBRATOR_1ST_PASS.out.versions)
    versions   = versions.mix(APPLYBQSR.out.versions)

    if (get_metrics) {
        // Second pass of Base Quality Score Recalibration (BQSR)
        // Generate error model based on corrected BAM file
        // Only used to Analyze Covariates
        BASERECALIBRATOR_2ND_PASS(
            bam_recall
                .join(bai_recall)
                .combine(intervals)
                .map { meta, bam, bai, bed ->
                    tuple(meta, bam, bai, bed)
                },
            fasta,
            fasta_fai,
            dict,
            known_sites,
            known_sites_tbi
        )

        // Analyze covariates (visually see the effects of recalibration)
        ANALYZECOVARIATES(
            BASERECALIBRATOR_1ST_PASS.out.table
                .join(BASERECALIBRATOR_2ND_PASS.out.table)
        )

        plots    = ANALYZECOVARIATES.out.plots
        versions = versions.mix(BASERECALIBRATOR_2ND_PASS.out.versions)
        versions = versions.mix(ANALYZECOVARIATES.out.versions)
    }

    emit:
    bam      = bam_recall  // channel: tuple(meta, bam)
    bai      = bai_recall  // channel: tuple(meta, bai)
    plots    = plots       // channel: tuple(meta, pdf)
    versions = versions    // channel: path(versions.yml)
}

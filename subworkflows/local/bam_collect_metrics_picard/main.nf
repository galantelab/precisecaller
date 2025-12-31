// Collect alignment and insert size metrics

include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../../../modules/nf-core/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS       } from '../../../modules/nf-core/picard/collectinsertsizemetrics/main'

workflow BAM_COLLECT_METRICS_PICARD {
    take:
    bam   // [mandatory] channel: tuple(meta, bam)
    fasta // [mandatory] channel: tuple(meta, fasta)

    main:
    versions = Channel.empty()

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(bam, fasta)
    PICARD_COLLECTINSERTSIZEMETRICS(bam)

    versions = versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions)
    versions = versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)

    emit:
    alignment_metrics = PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.metrics  // channel: tuple(meta, txt)
    insert_metrics    = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics        // channel: tuple(meta, txt)
    insert_histogram  = PICARD_COLLECTINSERTSIZEMETRICS.out.histogram      // channel: tuple(meta, pdf)
    versions          = versions                                           // channel: path(versions.yml)
}

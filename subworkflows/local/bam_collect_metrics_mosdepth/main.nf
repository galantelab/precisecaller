// Collect coverage metrics

include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'

workflow BAM_COLLECT_METRICS_MOSDEPTH {
    take:
    bam         // [mandatory] channel: tuple(meta, bam)
    bai         // [mandatory] channel: tuple(meta, bai)
    intervals   // [mandatory] channel: path(bed)

    main:
    versions = Channel.empty()

    MOSDEPTH(
        bam.join(bai)
           .combine(intervals)
           .map { meta, bam, bai, bed ->
               tuple(meta, bam, bai, bed)
           },
        [[id:"none"], []]
    )

    versions = versions.mix(MOSDEPTH.out.versions)

    emit:
    global_txt     = MOSDEPTH.out.global_txt      // channel: tuple(meta, txt)
    summary_txt    = MOSDEPTH.out.summary_txt     // channel: tuple(meta, txt)
    regions_txt    = MOSDEPTH.out.regions_txt     // channel: tuple(meta, txt)
    regions_bed    = MOSDEPTH.out.regions_bed     // channel: tuple(meta, bed)
    quantized_bed  = MOSDEPTH.out.quantized_bed   // channel: tuple(meta, bed)
    thresholds_bed = MOSDEPTH.out.thresholds_bed  // channel: tuple(meta, bed)
    versions       = versions                     // channel: path(versions.yml)
}

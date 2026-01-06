// Prepare intervals BED file
//
// This subworkflow normalizes user-provided or reference-derived interval
// definitions into a sorted BED file and its corresponding tabix index
//
// Supported input formats:
//   - .fai (converted to whole-genome BED; typically for WGS)
//   - .bed
//   - .bed.gz
//   - .list / .intervals (GATK-style chr:start-end)
//   - .interval_list (Picard-style)
//

include { CREATE_INTERVALS                      } from '../../../modules/local/create_intervals/main'
include { GUNZIP            as GUNZIP_INTERVALS } from '../../../modules/nf-core/gunzip/main'
include { TABIX_BGZIPTABIX  as TABIX_INTERVALS  } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow PREPARE_INTERVALS {
    take:
    fasta_fai   // [mandatory] channel: tuple(meta, fai)
    intervals   // [optional]  value: string - list/intervals, interval_list, bed

    main:
    versions             = Channel.empty()
    intervals_bed        = Channel.empty()
    intervals_bed_gz_tbi = Channel.empty()

    // Determine the source of the intervals
    //
    // If the user provides an intervals file, it is wrapped into a single-element
    // channel with a dummy meta object, since downstream nf-core modules expect
    // tuple(meta, path) inputs
    //
    // If no intervals are provided, we fall back to the FASTA index (FAI),
    // which is expected to already be a channel in tuple(meta, path) form
    intervals_path = Channel.empty()
    if (intervals != null) {
        intervals_path = Channel.value(
            [[id:"intervals"], file(intervals, checkIfExists:true)]
        )

        if (intervals.endsWith('.gz')) {
            GUNZIP_INTERVALS(intervals_path)
            intervals_path = GUNZIP_INTERVALS.out.gunzip.first()
            versions       = versions.mix(GUNZIP_INTERVALS.out.versions)
        }
    } else {
        intervals_path = fasta_fai
    }

    // Normalize all supported input formats into a sorted BED file
    CREATE_INTERVALS(intervals_path.map { it[1] })
    intervals_bed = CREATE_INTERVALS.out.bed.first()
    versions      = versions.mix(CREATE_INTERVALS.out.versions)

    // Compress and index the BED file using bgzip + tabix
    TABIX_INTERVALS(intervals_bed.map { [[id:"intervals"], it] })
    intervals_bed_gz_tbi = TABIX_INTERVALS.out.gz_index.map { meta, bed, tbi -> [bed, tbi] }.first()
    versions             = versions.mix(TABIX_INTERVALS.out.versions)

    emit:
    intervals        = intervals_bed         // channel: path(bed)
    intervals_gz_tbi = intervals_bed_gz_tbi  // channel: tuple(bed.gz, bed.gz.tbi)
    versions         = versions              // channel: path(versions.yml)
}

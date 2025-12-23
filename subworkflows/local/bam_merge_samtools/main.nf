// Merge multiple BAM files belonging to the same sample into a single BAM
// Metadata is normalized to sample-level after the merge; lane- or chunk-
// specific attributes are intentionally discarded

include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'

workflow BAM_MERGE_SAMTOOLS {
    take:
    bam       // channel: [mandatory] [ val(meta), bam ]

    main:
    versions = Channel.empty()

    grouped = bam
        .map { meta, bam_file ->
            [ meta.sample, bam_file ]
        }
        .groupTuple()
        .map { sample, bam_files ->
            def sorted_bams = bam_files.sort { a, b ->
                a.getName() <=> b.getName()
            }

            // Normalize metadata to sample-level after merging lanes/chunks
            def meta = [
                id:     sample,
                sample: sample,
                n_bams: sorted_bams.size()
            ]

            [ meta, sorted_bams ]
        }

    SAMTOOLS_MERGE(
        grouped,
        [[id:'no_fasta'], []],
        [[id:'no_fai'], []],
        [[id:'no_gzi'], []]
    )

    versions = versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    bam      = SAMTOOLS_MERGE.out.bam   // channel: [ val(meta), [ bam ] ]
    versions = versions                 // channel: [ versions.yml ]
}

// Subworkflow responsible for splitting FASTQ files using seqkit split2
// and normalizing the output into paired FASTQ chunks suitable for downstream
// alignment. Each chunk represents a consistent R1/R2 pair with updated metadata

include { SEQKIT_SPLIT2 } from '../../../modules/nf-core/seqkit/split2/main'

workflow FASTQ_SPLIT_SEQKIT {
    take:
    reads  // channel: [mandatory] [ val(meta), [ reads ] ]

    main:
    versions = Channel.empty()

    SEQKIT_SPLIT2(reads)

    // Regroup split FASTQs into consistent chunks
    // Each chunk corresponds to a single split unit (e.g. part_001)
    regrouped = SEQKIT_SPLIT2.out.reads
        .flatMap { meta, files ->
            files.collect { f ->
                // Expected filename pattern: *_R1.part_001.fq.gz
                // Regex must capture chunk id as group(1)
                def matcher = f.getName() =~ /[rR]?[12]\.(part_\d+)\./

                assert matcher.find() :
                    "File ${f.getName()} does not match split FASTQ pattern"

                [ meta + [ chunk: matcher.group(1) ], f ]
            }
        }
        .groupTuple(sort:true)
        .map { meta, pair ->
            if (meta.single_end)
                assert pair.size() == 1 :
                    "Expected 1 FASTQ for ${meta.id}:${meta.chunk}, got ${pair.size()}"
            else
                assert pair.size() == 2 :
                    "Expected 2 FASTQs for ${meta.id}:${meta.chunk}, got ${pair.size()}"
            [ meta, pair ]
        }

    versions = versions.mix(SEQKIT_SPLIT2.out.versions)

    emit:
    reads    = regrouped  // [ val(meta), [ reads ] ]
    versions = versions  // channel: [ versions.yml ]
}

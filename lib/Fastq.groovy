import nextflow.Channel
import java.util.regex.Pattern

class Fastq {

    public static Channel regroupByPattern(Channel input_ch, Pattern pattern) {
        input_ch
            .flatMap { meta, reads ->
                reads.collect { read ->
                    def matcher = read.getName() =~ pattern

                    assert matcher.find() :
                        "File ${read.getName()} does not match split FASTQ pattern"

                    assert matcher.groupCount() >= 1 :
                        "Pattern must capture chunk id as group(1)"

                    def chunk = matcher.group(1)
                    [ meta + [ chunk: chunk ], read ]
                }
            }
            .groupTuple(sort: true)
            .map { meta, pair ->
                if (meta.single_end) {
                    assert pair.size() == 1 :
                        "Expected 1 FASTQ for ${meta.id}:${meta.chunk}, got ${pair.size()}"
                } else {
                    assert pair.size() == 2 :
                        "Expected 2 FASTQs for ${meta.id}:${meta.chunk}, got ${pair.size()}"
                }
                [ meta, pair ]
            }
    }

}

import java.nio.file.Path
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

class ReadGroup {

    public static Map extractFromFastq(
        Map meta,
        List<Path> files,
        String seq_center,
        String seq_platform
    ) {
        assert files.size() >= 1 :
            "No FASTQ files provided for sample ${meta.id}"

        assert seq_platform :
            "--seq_platform is a mandatory option"

        def attrs1 = this.readNameAttrs(files[0])

        if (!meta.single_end) {
            def attrs2 = this.readNameAttrs(files[1])
            assert attrs1.flowcell == attrs2.flowcell :
                "Flowcell ID does not match for paired reads of sample ${meta.id}: ${files*.getName()}"
        }

        def ID = attrs1.flowcell ?
            "${attrs1.flowcell}.${meta.lane}" :
            "${meta.sample}.${meta.lane}"

        def PU = attrs1.barcode ?
            "${ID}.${attrs1.barcode}" :
            "${ID}"

        return [
            ID: ID,
            PU: PU,
            CN: seq_center,
            PL: seq_platform,
            SM: meta.sample,
            LB: "${meta.sample}_lib"
        ]
    }

    public static String toSam(Map rg) {
    def order = ['ID', 'CN', 'PU', 'SM', 'LB', 'PL']
    
    return '@RG\t' + order
        .findAll { rg[it] != null }
        .collect { "${it}:${rg[it]}" }
        .join('\t')
    }

    private static Map readNameAttrs(Path path) {
        def firstLine = this.readFirstLineOfFastq(path)

        def parts  = []
        def fields = []

        if (firstLine) {
            parts  = firstLine.split(/\s/)
            fields = parts[0].split(':')
        }

        def attrs = [flowcell: null, lane: null, barcode: null]

        if (fields.size() >= 7) {  // Illumina 1.8
            attrs.flowcell = fields[2]
            attrs.lane = fields[3]
            if (parts.size() > 1) {
                def fields2 = parts[1].split(':')
                if (fields2.size() >= 4) {
                    attrs.barcode = fields2[3]
                }
            }
        } else if (fields.size() == 5) {  // Broad 1.0, Illumina 1.0, Illumina 1.4
            def prefix = fields[0].substring(1)
            attrs.lane = fields[1]
            if (fields[4].contains('#')) { // Illumina 1.0, Illumina 1.4
                attrs.flowcell = prefix
                attrs.barcode = fields[4].split('#')[1].split('/')[0]
            } else { // Broad 1.0
                attrs.flowcell = prefix.substring(0, Math.min(5, prefix.length()))
                attrs.barcode = prefix.substring(Math.min(5, prefix.length()))
            }
        } // fields.size() != 0 NCBI SRA and others

        return attrs
    }

    private static String readFirstLineOfFastq(Path path) {
        def firstLine = null
        try {
            path.withInputStream { stream ->
                if (path.toString().endsWith('.gz')) {
                    stream = new java.util.zip.GZIPInputStream(stream)
                }
                Reader decoder = new InputStreamReader(stream, 'ASCII')
                BufferedReader buffered = new BufferedReader(decoder)
                firstLine = buffered.readLine()
                assert firstLine.startsWith('@')
            }
        } catch (Exception e) {
            throw new RuntimeException(
                "Error reading FASTQ file ${path}: ${e.message}", e
            )
        }
        return firstLine
    }

}

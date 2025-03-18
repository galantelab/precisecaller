class Utils {

    public static String extractReadGroup(meta, files, seq_center, seq_platform, log) {
        def CN = seq_center ? "CN:${seq_center}\\t" : ''

        if (!seq_platform) {
           log.error("--seq_platform is a mandatory option")
           return
        }

        def attrs1 = this.readNameAttrs(files[0], log)
        def attrs2 = this.readNameAttrs(files[1], log)

        if (attrs1.flowcell != attrs2.flowcell) {
           log.error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
           return
        }

        def ID = attrs1.flowcell ?
            "${attrs1.flowcell}.${meta.lane}" :
            "${meta.sample}.${meta.lane}"

        def PU = attrs1.barcode ?
            "${ID}.${attrs1.barcode}" :
            "${ID}"

        def RG = "\"@RG\\tID:${ID}\\t${CN}PU:${PU}\\tSM:${meta.sample}\\tLB:${meta.sample}_lib\\tPL:${seq_platform}\""
        return RG
    }

    private static Map readNameAttrs(path, log) {
        def firstLine = this.readFirstLineOfFastq(path, log)

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
        } else if (fields.size() != 0) {  // NCBI SRA and others
            log.warn "FASTQ file(${path}): Cannot extract flowcell ID from ${firstLine}"
        }

        return attrs
    }

    private static String readFirstLineOfFastq(path, log) {
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
            log.warn "FASTQ file(${path}): Error streaming"
            log.warn "${e.message}"
        }
        return firstLine
    }

}

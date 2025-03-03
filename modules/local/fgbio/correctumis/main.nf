process FGBIO_CORRECTUMIS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/87/87626ef674e2f19366ae6214575a114fe80ce598e796894820550731706a84be/data' :
        'community.wave.seqera.io/library/fgbio:2.4.0--913bad9d47ff8ddc' }"

    input:
    tuple val(meta), path(bam)
    path expected_umis

    output:
    tuple val(meta), path("*.correct_umis.bam") , emit: bam,  optional: true
    tuple val(meta), path("*.correct_umis.cram"), emit: cram, optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    def max_mismatches = args.contains("--max-mismatches") ? "" : "--max-mismatches 1"
    def min_distance = args.contains("--min-distance") ? "" : "--min-distance 2"

    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FastqToBam] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        CorrectUmis \\
        ${args} \\
        --input ${bam} \\
        --output ${prefix}.correct_umis.${suffix} \\
        --umi-files ${expected_umis} \\
        --dont-store-original-umis \\
        ${max_mismatches} \\
        ${min_distance}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    """
    touch ${prefix}.correct_umis.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

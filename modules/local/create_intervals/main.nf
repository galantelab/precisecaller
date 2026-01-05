process CREATE_INTERVALS {
    tag "${intervals.getName()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    path intervals

    output:
    path "*.bed",        emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${intervals.getBaseName()}.intervals"

    if (intervals.getName().endsWith('.fai')) {
        """
        awk -F '\\t' 'BEGIN {OFS=FS} {print \$1,0,\$2}' ${intervals} \\
            | sort -k 1,1 -k 2,2n -k 3,3n > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else if (intervals.getName().endsWith('.list') || intervals.getName().endsWith('.intervals')) {
        """
        awk -F '[:-]' 'BEGIN {OFS="\\t"} NF>=3 {print \$1,\$2-1,\$3}' ${intervals} \\
            | sort -k 1,1 -k 2,2n -k 3,3n > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else if (intervals.getName().endsWith('.interval_list')) {
        """
        awk 'BEGIN {OFS="\\t"} !/^@/ {print \$1,\$2-1,\$3}' ${intervals} \\
            | sort -k 1,1 -k 2,2n -k 3,3n > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else if (intervals.getName().endsWith('.bed')) {
        """
        sort -k 1,1 -k 2,2n -k 3,3n ${intervals} > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else {
        """
        echo "ERROR: Unsupported interval format: ${intervals}" >&2
        exit 1
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${intervals.getBaseName()}.intervals"

    """
    touch ${prefix}.stub.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}

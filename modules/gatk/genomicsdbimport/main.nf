process GATK_GENOMICSDBIMPORT {
    tag "${meta.label}"
    label 'gatk_genomicsdbimport'

    input:
    tuple val(meta), path(vcfs)
    path(sample_map)

    output:
    tuple val(meta), file("${meta.label}.db"), emit: db
    path  "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def avail_mem = (task.memory.giga).intValue() - 2
    def args = task.ext.args ?: ''
    """
    gatk  --java-options "-Xmx${avail_mem}g -XX:ConcGCThreads=${task.cpus}" \\
        GenomicsDBImport \\
        --genomicsdb-workspace-path ${meta.label}.db \\
        --batch-size 16 \\
        -L ${meta.interval} \\
        --sample-name-map ${sample_map}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.label}.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
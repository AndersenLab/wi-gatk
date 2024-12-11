process GATK_GENOTYPEGVCFS {
    label 'gatk_genotypegvcfs'

    input:
    tuple val(meta), path(contig_db)
    tuple val(meta2), path("ref.fa.gz"), path("ref.fa.gz.fai"), path("ref.dict"), path("ref.fa.gz.gzi")

    output:
    tuple val(meta), path("${meta.label}.vcf"), emit: vcf
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.giga*0.9).intValue()
    """
    gatk  --java-options "-Xmx${avail_mem}g -XX:ConcGCThreads=${task.cpus} -XX:-UsePerfData" \\
        GenotypeGVCFs \\
        -R ref.fa.gz \\
        -V gendb://${contig_db} \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        -G StandardHCAnnotation \\
        -L ${meta.interval} \\
        -O ${meta.label}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.label}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
 
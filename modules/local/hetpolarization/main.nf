process LOCAL_HETPOLARIZATION {
    tag "${meta.label}"
    label 'local_hetpolarization'
    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 32.GB * task.attempt }

    input:
    tuple val(meta), path(vcf)
    val mito_name

    output:
    tuple val(meta), path("${meta.label}.vcf.gz"), path("${meta.label}.vcf.gz.tbi"), emit: vcf
    path  "versions.yml",                                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    if [ "${meta.contig}" == "${mito_name}" ]
    then
        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} ${vcf} | \\
        vcffixup - | \\
        bcftools view -O z - > ${meta.label}.vcf.gz
    else
        bcftools view -O z ${vcf} | \\
        het_polarization | \\
        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} | \\
        vcffixup - | \\
        bcftools view -O z - > ${meta.label}.vcf.gz
    fi
    bcftools index -t ${meta.label}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.label}.vcf.gz
    touch ${meta.label}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
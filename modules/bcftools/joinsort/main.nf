process BCFTOOLS_JOINSORT {
    label 'bcftools_joinsort'
    errorStrategy 'retry'
    time { 40.minute * task.attempt }
    cpus { 1 * task.attempt }
    memory { 4.GB * task.attempt }

    input:
    tuple val(meta), path(vcf1), path(index1), path(vcf2), path(index2)

    output:
    tuple val(meta), path("${meta.label}.all.vcf.gz"), path("${meta.label}.all.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools concat -a -O b --threads ${task.cpus} ${vcf1} ${vcf2} |
    bcftools sort -O z -o ${meta.label}.all.vcf.gz
    bcftools index --tbi ${meta.label}.all.vcf.gz
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.label}.all.vcf.gz
    touch ${meta.label}.all.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
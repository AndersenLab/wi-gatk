process BCFTOOLS_COMPRESS {
    label 'bcftools_compress'
    errorStrategy 'retry'
    time { 20.minute * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${vcf}.gz"), path("${vcf}.gz.tbi"), emit: vcf
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view -O z ${vcf} -o ${vcf}.gz --threads ${task.cpus}
    bcftools index --tbi ${vcf}.gz
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${vcf}.gz
    touch ${vcf}.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
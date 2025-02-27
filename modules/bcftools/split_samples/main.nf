process BCFTOOLS_SPLIT_SAMPLES {

    tag "${sample}"
    label 'bcftools_split_samples'
    errorStrategy 'retry'
    time { 1.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    val sample
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools view -c1 -Oz -s ${sample} -o ${sample}.vcf.gz ${vcf}
    bcftools index --tbi ${sample}.vcf.gz
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample}.vcf.gz
    touch ${sample}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
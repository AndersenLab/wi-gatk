process LOCAL_SOFT_STATS {
    label 'local_soft_stats'

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    path "soft_filter.stats", emit: stats
    path "versions.yml",      emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools stats --threads ${task.cpus} \\
                    -s- --verbose ${vcf} > soft_filter.stats

    {
        echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' ${vcf};
    }     > soft_filter.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch soft_filter.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}

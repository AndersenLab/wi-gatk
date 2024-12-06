process LOCAL_VCF_STATS {
    label 'local_vcf_stats'

    input:
    tuple val(meta), path(soft_vcf), path(soft_index)
    tuple val(meta2), path(hard_vcf), path(hard_index)

    output:
    path "soft_filter.stats.txt",        emit: soft_stats
    path "hard_filter.stats.txt",        emit: hard_stats
    path "soft_filter.filter_stats.txt", emit: soft_filter_stats
    path "versions.yml",                 emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bcftools stats --threads ${task.cpus} \\
                    -s- --verbose ${soft_vcf} > soft_filter.stats.txt

    {
        echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' ${soft_vcf};
    }     > soft_filter.filter_stats.txt

    bcftools stats -s- --verbose ${hard_vcf} > hard_filter.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch soft_filter.stats.txt
    touch hard_filter.stats.txt
    touch soft_filter.filter_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}

process BCFTOOLS_CONCAT_OVERLAPPING_VCFS {
    tag "${meta.contig} ${meta.id}"
    label 'bcftools_concat_vcfs'
    errorStrategy 'retry'
    time { 4.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
    tuple val(meta), path(vcfs), path(indices)
    path partitions

    output:
    tuple val(meta), path("${meta.contig}.${meta.id}.vcf.gz"), path("${meta.contig}.${meta.id}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    grep -w ${meta.contig} ${partitions} | \\
        awk '{ printf "%s_%i_%i.${meta.id}.vcf.gz\\n", \$1, \$3, \$4 }' > partition_set.tsv
    bcftools concat -a -d all -O z --file-list partition_set.tsv > ${meta.contig}.${meta.id}.vcf.gz
    bcftools index --tbi ${meta.contig}.${meta.id}.vcf.gz
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.contig}.${meta.id}.vcf.gz
    touch ${meta.contig}.${meta.id}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
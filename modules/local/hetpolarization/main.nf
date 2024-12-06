process LOCAL_HETPOLARIZATION {
    label 'local_hetpolarization'
    errorStrategy 'retry'
    time { 2.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
    tuple val(meta), path(vcfs), path(indices)
    path partitions
    val mito_name

    output:
    tuple val(meta), path("${meta.contig}.${meta.id}.vcf.gz"), path("${meta.contig}.${meta.id}.vcf.gz.tbi"), emit: vcf
    path  "versions.yml",                                                                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    grep -w ${meta.contig} ${partitions} | awk '{printf "%s_%i_%i.vcf.gz\\n", \$1, \$3, \$4}' > sample_list.txt


    if [ "${meta.contig}" == "${mito_name}" ]
    then
        #bcftools index -c tmp.bcf
        bcftools concat -a -d all -f sample_list.txt -O b | \\
        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} - | \\
        vcffixup - | \\
        bcftools view -O z - > ${meta.contig}.${meta.id}.vcf.gz
    else
        #bcftools index -c tmp.bcf
        bcftools concat -a -d all -f sample_list.txt -O z | \\
        het_polarization | \\
        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} | \\
        vcffixup - | \\
        bcftools view -O z - > ${meta.contig}.${meta.id}.vcf.gz
    fi
    bcftools index -t ${meta.contig}.${meta.id}.vcf.gz

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
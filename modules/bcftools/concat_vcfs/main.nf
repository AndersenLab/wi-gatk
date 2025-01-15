process BCFTOOLS_CONCAT_VCFS {
    tag "${meta.id}"
    label 'bcftools_concat_vcfs'
    errorStrategy 'retry'
    time { 3.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    tuple val(meta), path(vcfs)
    path(contigs)

    output:
    tuple val(meta), path("${meta.id}*.vcf.gz"), path("${meta.id}*.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    CONTIGS=(`cat ${contigs}`)
    for C in \${CONTIGS[*]}; do
        ls \${C}.${meta.id}*.vcf* >> contig_set.tsv
    done
    SUFFIX=`ls \${CONTIGS[0]}.${meta.id}*.vcf* | sed "s/\${CONTIGS[0]}\\.${meta.id}\\.//" | sed "s/\\.gz//"`
    SUFFIX="\${SUFFIX}.gz"
    bcftools concat  -O z --file-list contig_set.tsv > ${meta.id}.\${SUFFIX}
    bcftools index --tbi ${meta.id}.\${SUFFIX}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    CONTIGS=(`cat ${contigs}`)
    SUFFIX=`ls \${CONTIGS[0]}.${meta.id}*.vcf* | sed "s/\${CONTIGS[0]}\\.${meta.id}\\.//"`
    touch ${meta.id}.\${SUFFIX}
    touch ${meta.id}.\${SUFFIX}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
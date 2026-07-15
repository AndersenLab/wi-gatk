process LOCAL_HETPOLARIZATION {
    tag "${meta.label}"
    label 'local_hetpolarization'
    errorStrategy 'retry'
    time { 6.hour * task.attempt }
    cpus { 1 * task.attempt }
    memory { 31.GB * task.attempt }

    input:
    tuple val(meta), path(vcf)
    val mito_name

    output:
    tuple val(meta), path("${meta.label}.*variant.vcf.gz"), path("${meta.label}.*variant.vcf.gz.tbi"), emit: vcf
    path  "versions.yml",                                                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threads = task.cpus > 1 ? "--threads=${task.cpus - 1}" : ""
    """
    if [ "${meta.sites}" == "invariant" ]; then
        bcftools view -O v --max-ac 1 ${threads} ${vcf} | \\
        vcffixup - | \\
        bcftools view -O z - > ${meta.label}.invariant.vcf.gz
        bcftools index -t ${meta.label}.invariant.vcf.gz
    else
        if [ "${meta.contig}" == "${mito_name}" ]
        then
            bcftools view -O v --min-af 0.000001 ${threads} ${vcf} | \\
            vcffixup - | \\
            bcftools view -O z - > ${meta.label}.variant.vcf.gz
        else
            bcftools view -O z --min-ac 2 ${vcf} | \\
            het_polarization | \\
            bcftools view -O v --min-af 0.000001 ${threads} | \\
            vcffixup - | \\
            bcftools view -O z - > ${meta.label}.variant.vcf.gz
        fi
        bcftools index -t ${meta.label}.variant.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    if [ "${meta.sites}" == "invariant" ]; then
        touch ${meta.label}.invariant.vcf.gz
        touch ${meta.label}.invariant.vcf.gz.tbi
    else
        touch ${meta.label}.variant.vcf.gz
        touch ${meta.label}.variant.vcf.gz.tbi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
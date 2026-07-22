process GTCHECK {

    label 'bcftools_gtcheck'
    errorStrategy 'retry'
    time { 6.hour * task.attempt }
    cpus { 1 * task.attempt }
    memory { 31.GB * task.attempt }

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(strategy), path(full_vcf), path(full_vcf_index)

    output:
    path "gtcheck.tsv"  , emit: gt
    path "versions.yml" , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def strain_option = strategy == "full" ? "-g ${full_vcf}" : ""
    """
    bcftools gtcheck --no-HWE-prob -e 0 ${vcf} ${strain_option} -o gtcheck.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_gtcheck.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}

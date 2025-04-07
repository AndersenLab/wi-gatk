process GATK_SOFT_FILTER {
    tag "${meta.label}"
    label 'gatk_soft_filter'
    errorStrategy 'retry'
    time { 2.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
    val filter_params
    tuple val(meta), path(vcf), path(vcf_index)
    tuple val(meta2), path("ref.fa.gz"), path("ref.fa.gz.fai"), path("ref.dict"), path("ref.fa.gz.gzi")

    output:
    tuple val(meta), path("${meta.label}.gatksoft.vcf"), emit: vcf
    path  "versions.yml",                                emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.giga*0.9).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}g -XX:+UseSerialGC" \\
        VariantFiltration \\
        -R ref.fa.gz \\
        --variant ${vcf} \\
        --genotype-filter-expression "DP < ${filter_params.min_depth}"    --genotype-filter-name "DP_min_depth" \\
        --filter-expression "QUAL < ${filter_params.qual}"                --filter-name "QUAL_quality" \\
        --filter-expression "FS > ${filter_params.fisherstrand}"          --filter-name "FS_fisherstrand" \\
        --filter-expression "QD < ${filter_params.quality_by_depth}"      --filter-name "QD_quality_by_depth" \\
        --filter-expression "SOR > ${filter_params.strand_odds_ratio}"    --filter-name "SOR_strand_odds_ratio" \\
        --genotype-filter-expression "isHet == 1"                         --genotype-filter-name "is_het" \\
        -O ${meta.label}.gatksoft.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.label}.gatksoft.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

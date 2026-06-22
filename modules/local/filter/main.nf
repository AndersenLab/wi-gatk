process LOCAL_FILTER {
    tag "${meta.label}"
    label 'local_filter'
    errorStrategy 'retry'
    time { 1.hour * task.attempt }
    cpus { 1 * task.attempt }
    memory { 4.GB * task.attempt }

    input:
    val filter_params
    tuple val(meta), path(vcf)
    tuple val(meta2), path("ref.fa.gz"), path("ref.fa.gz.fai"), path("ref.dict"), path("ref.fa.gz.gzi")
    val mito_name
    val all_sites

    output:
    tuple val(meta), path("${meta.label}.soft.vcf.gz"), path("${meta.label}.soft.vcf.gz.tbi"), emit: soft, optional: true
    tuple val(meta), path("${meta.label}.hard.vcf.gz"), path("${meta.label}.hard.vcf.gz.tbi"), emit: hard, optional: true
    tuple val(meta), path("${meta.label}.allsites.vcf.gz"), path("${meta.label}.allsites.vcf.gz.tbi"), emit: allsites, optional: true
    path  "versions.yml",                                                                      emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.giga*0.9).intValue()
    def all_sites_arg = all_sites ? '' : '--min-af 0.000001 --max-af 0.999999'
    def min_alleles_arg = all_sites ? '-m1' : '-m2'
    """
    bcftools view -O z ${vcf} > out.vcf.gz
    bcftools index --tbi out.vcf.gz
    
    # Apply high missing and high heterozygosity filters
    bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING  <= ${filter_params.high_missing}' out.vcf.gz |\\
    bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${filter_params.high_heterozygosity}' -O z > ${meta.label}.soft.vcf.gz
   
    # Create hard-filtered version by removing marked genotypes/variants
    if [ ${meta.contig} == "${mito_name}" ]
    then
        bcftools view ${min_alleles_arg} -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
        bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' | \\
        bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
        bcftools view -O v ${all_sites_arg} | \\
        vcffixup - | \\
        bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
    else
        bcftools view ${min_alleles_arg} -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
        bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' | \\
        bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
        bcftools view -O v ${all_sites_arg} | \\
        vcffixup - | \\
        bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
    fi

    if [[ ${all_sites} == "true" ]]; then
        bcftools annotate -x INFO,^FORMAT/GT -O z ${meta.label}.hard.vcf.gz > ${meta.label}.allsites.vcf.gz
        bcftools index --tbi ${meta.label}.allsites.vcf.gz
        rm ${meta.label}.soft.vcf.gz
        rm ${meta.label}.hard.vcf.gz
    else
        bcftools index --tbi ${meta.label}.soft.vcf.gz
        bcftools index --tbi ${meta.label}.hard.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    if [[ ${all_sites} == "true" ]]; then
        touch ${meta.label}.allsites.vcf.gz
        touch ${meta.label}.allsites.vcf.gz.tbi
    else
        touch ${meta.label}.soft.vcf.gz
        touch ${meta.label}.soft.vcf.gz.tbi
        touch ${meta.label}.hard.vcf.gz
        touch ${meta.label}.hard.vcf.gz.tbi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}

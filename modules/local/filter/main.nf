process LOCAL_FILTER {
    label 'local_filter'
    errorStrategy 'retry'
    time { 2.hour * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 8.GB * task.attempt }

    input:
    val filter_params
    tuple val(meta), path(vcf)
    tuple val(meta2), path("ref.fa.gz"), path("ref.fa.gz.fai"), path("ref.dict"), path("ref.fa.gz.gzi")
    val mito_name

    output:
    tuple val(meta), path("${meta.contig}.soft.vcf.gz"), emit: soft
    tuple val(meta), path("${meta.contig}.hard.vcf.gz"), emit: hard
    path  "versions.yml",                                emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.giga*0.9).intValue()
    """
    bcftools view -O z ${vcf} > out.vcf.gz
    bcftools index --tbi out.vcf.gz
    
    # Apply high missing and high heterozygosity filters
    bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING  <= ${filter_params.high_missing}' out.vcf.gz |\\
    bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${filter_params.high_heterozygosity}' -O z > ${meta.contig}.soft.vcf.gz
   
    # Create hard-filtered version by removing marked genotypes/variants
    if [ ${meta.contig} == "${mito_name}" ]
    then
        bcftools view -m2 -M2 --trim-alt-alleles -O u ${meta.contig}.soft.vcf.gz | \\
        bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' | \\
        bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
        bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
        vcffixup - | \\
        bcftools view -O z --trim-alt-alleles > ${meta.contig}.hard.vcf.gz
    else
        bcftools view -m2 -M2 --trim-alt-alleles -O u ${meta.contig}.soft.vcf.gz | \\
        bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' | \\
        bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
        bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
        vcffixup - | \\
        bcftools view -O z --trim-alt-alleles > ${meta.contig}.hard.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.contig}.soft.vcf
    touch ${meta.contig}.hard.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}

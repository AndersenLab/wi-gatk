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

    output:
    tuple val(meta), path("${meta.label}.soft.vcf.gz"), path("${meta.label}.soft.vcf.gz.tbi"),                   emit: soft, optional: true
    tuple val(meta), path("${meta.label}.hard.vcf.gz"), path("${meta.label}.hard.vcf.gz.tbi"),                   emit: hard, optional: true
    tuple val(meta), path("${meta.label}.${meta.sites}.vcf.gz"), path("${meta.label}.${meta.sites}.vcf.gz.tbi"), emit: allsites, optional: true
    path  "versions.yml",                                                                                        emit: versions
    
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
    bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${filter_params.high_heterozygosity}' -O z > ${meta.label}.soft.vcf.gz
   
    # Create hard-filtered version by removing marked genotypes/variants
    if [ "${meta.sites}" == "invariant" ]; then
        if [ ${meta.contig} == "${mito_name}" ]
        then
            bcftools view -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' | \\
            bcftools filter -O v --exclude 'FILTER != "PASS"' | \\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
        else
            bcftools view -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' | \\
            bcftools filter -O v --exclude 'FILTER != "PASS"' | \\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
        fi
    else
        if [ ${meta.contig} == "${mito_name}" ]
        then
            bcftools view -m2 -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' | \\
            bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
            bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
        else
            bcftools view -m2 -M2 --trim-alt-alleles -O u ${meta.label}.soft.vcf.gz | \\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' | \\
            bcftools filter -O u --exclude 'FILTER != "PASS"' | \\
            bcftools view -O v --min-af 0.000001 --max-af 0.999999 | \\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > ${meta.label}.hard.vcf.gz
        fi
    fi

    if [ "${meta.all_sites}" == "true" ]; then
        bcftools annotate -x INFO,^FORMAT/GT -O z ${meta.label}.hard.vcf.gz > ${meta.label}.${meta.sites}.vcf.gz
        bcftools index --tbi ${meta.label}.${meta.sites}.vcf.gz
    fi
    if [ "${meta.sites}" == "invariant" ]; then
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
    if [ ${meta.all_sites} == "true" ]; then
        touch ${meta.label}.${meta.sites}.vcf.gz
        touch ${meta.label}.${meta.sites}.vcf.gz.tbi
    fi
    if [ ${meta.sites} != "invariant" ]; then
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

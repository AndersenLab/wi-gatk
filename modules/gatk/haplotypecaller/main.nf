process GATK_HAPLOTYPECALLER {
    tag "${meta.id} ${meta.contig}"
    label 'gatk_haplotypecaller'
    //errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    errorStrategy 'retry'
    time { 120.minute * task.attempt }
    cpus = { 2 * task.attempt }
    memory = { 4.GB * task.attempt - 1.GB }
    input:
    tuple val(meta), path(bam), path(bam_index) 
    tuple val(meta2), path("ref.fa.gz"), path("ref.fa.gz.fai"), path("ref.dict"), path("ref.fa.gz.gzi")

    output:
    tuple val(meta), path("${meta.contig}.${meta.id}.g.vcf"), emit: vcf
    path  "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.giga - 1).intValue()
    """
    gatk HaplotypeCaller \\
        --java-options "-Xmx${avail_mem}g \\
                        -XX:ConcGCThreads=${task.cpus} \\
                        -XX:ParallelGCThreads=${task.cpus}" \\
        --emit-ref-confidence GVCF \\
        --annotation DepthPerAlleleBySample \\
        --annotation Coverage \\
        --annotation GenotypeSummaries \\
        --annotation TandemRepeat \\
        --annotation StrandBiasBySample \\
        --annotation ChromosomeCounts \\
        --annotation ReadPosRankSumTest \\
        --annotation AS_ReadPosRankSumTest \\
        --annotation AS_QualByDepth \\
        --annotation AS_StrandOddsRatio \\
        --annotation AS_MappingQualityRankSumTest \\
        --annotation DepthPerSampleHC \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --annotation-group StandardHCAnnotation \\
        --do-not-run-physical-phasing \\
        -R ref.fa.gz \\
        -I ${meta.id}.bam \\
        -L ${meta.contig} \\
        -O ${meta.contig}.${meta.id}.g.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.contig}.${meta.id}.g.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

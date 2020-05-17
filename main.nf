#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Variant Calling Pipeline
    Authors:
    - Stefan Zdraljevic
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2

// For now, this pipeline requires NXF_VER 20.01.0-rc1
// Prefix this version when running
// e.g.
// NXF_VER=20.01.0-rc1 nextflow run ...
//assert System.getenv("NXF_VER") == "20.01.0-rc1"

/*
    Params
*/

date = new Date().format( 'yyyyMMdd' )
params.debug = false
params.help = false
params.bam_location = "" // Use this to specify the directory for bams

// Check that reference exists
params.reference = ""
reference = file(params.reference, checkIfExists: true)

// Debug
if (params.debug.toString() == "true") {
    params.output = "release-debug"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bam_location = "${workflow.projectDir}" // Use this to specify the directory for bams
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "WI-${date}"
    params.sample_sheet = "${workflow.projectDir}/sample_sheet.tsv"
}


// Variant Filtering
params.min_depth = 5
params.qual = 30.0
params.strand_odds_ratio = 5.0
params.dv_dp = 0.5
params.quality_by_depth = 5.0
params.fisherstrand = 50.0
params.high_missing = 0.95
params.high_heterozygosity = 0.10


def log_summary() {

    out =  '''
 _______ _______ _______ __  __      _______ _______ 
|     __|   _   |_     _|  |/  |    |    |  |    ___|
|    |  |       | |   | |     <     |       |    ___|
|_______|___|___| |___| |__|\\__|    |__|____|___|    
                                              
'''

out += """

    parameters               description                Set/Default
    ==========               ===========                ========================
    output                   Release Directory          ${params.output}
    sample_sheet             sample sheet               ${params.sample_sheet}
    reference                Reference Genome           ${reference}
    username                                            ${"whoami".execute().in.text}

    Nextflow Run
    ---------------
    ${workflow.commandLine}
    run name                                            ${workflow.runName}
    scriptID                                            ${workflow.scriptId}
    git commit                                          ${workflow.commitId}
    container                                           ${workflow.container}

    Reference Genome
    ---------------
    reference_base          location of ref genomes     ${params.reference_base}
    species/project/build                               ${params.species} / ${params.project} / ${params.ws_build}

    Variant Filters         
    ---------------           
    min_depth                Minimum variant depth      ${params.min_depth}
    qual                     Variant QUAL score         ${params.qual}
    ad_dp                    Good ALT reads / depth     ${params.dv_dp}
    strand_odds_ratio        SOR_strand_odds_ratio      ${params.strand_odds_ratio} 
    quality_by_depth         QD_quality_by_depth        ${params.quality_by_depth} 
    fisherstrand             FS_fisher_strand           ${params.fisherstrand}
    missing_max              % missing genotypes        ${params.high_missing}
    heterozygosity_max       % max heterozygosity       ${params.high_heterozygosity}

---
"""
out
}

log.info(log_summary())

if (params.help) {
    exit 1
}

if (workflow.profile == "") {
    println "Must set -profile: local, quest, gcp"
    exit 1
}


// Read sample sheet
sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> 
                                    // Optionally allow user to specify a bam location
                                    if (params.bam_location != "") {
                                        row.bam = "${params.bam_location}/${row.bam}"
                                    }           
                                    [row.strain,
                                     file("${row.bam}", checkExists: true),
                                     file("${row.bam}.bai", checkExists: true)] }


workflow {

    // Generate a summary of the current run
    summary(Channel.from("run"))

    // Get contigs from first bam
    sample_sheet.first() | get_contigs
    contigs = get_contigs.out.splitText { it.strip() }

    // Call individual variants
    sample_sheet.combine(contigs) | call_variants_individual

    call_variants_individual.out.groupTuple()
                                .map { strain, vcf -> [strain, vcf]}
                                .combine(get_contigs.out) | \
                                concat_strain_gvcfs
    
    // gatk genomics db
    sample_map = sample_sheet.map { "${it[0]}\t${it[0]}.g.vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)
    concat_strain_gvcfs.out.flatten()
                           .toList()
                           .map { [it] }
                           .combine(contigs)
                           .combine(sample_map) | \
                        import_genomics_db | \
                        genotype_cohort_gvcf_db | \
                        annotate_vcf

    vcfs = annotate_vcf.out.anno_vcf.collect().map { [it] }.combine(get_contigs.out)
    vcfs | concatenate_vcf
    concatenate_vcf.out.vcf | soft_filter
    soft_filter.out.soft_filter_vcf.combine(get_contigs.out) | hard_filter

    // Generate Strain-level TSV and VCFs
    soft_filter.out.soft_filter_vcf | strain_list
    strain_set = strain_list.out.splitText()
                   .map { it.trim() }

    strain_set.combine( soft_filter.out.soft_filter_vcf ) | generate_strain_tsv_vcf
    
    // Extract SNPeff severity tracks
    mod_tracks = Channel.from(["LOW", "MODERATE", "HIGH", "MODIFIER"])
    soft_filter.out.soft_filter_vcf.spread(mod_tracks) | generate_severity_tracks

    // MultiQC Report
    soft_filter.out.soft_vcf_stats.concat(
        hard_filter.out.hard_vcf_stats,
        annotate_vcf.out.snpeff_csv
    ).collect() | multiqc_report

}


process summary {
    // Generates a summary of the run for the release directory.
    
    executor 'local'

    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        path("sample_sheet.tsv")
        path("summary.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        cat ${params.sample_sheet} > sample_sheet.tsv
    """

}

/*=========================================
~ ~ ~ > * Generate Interval List  * < ~ ~ ~ 
=========================================*/

process get_contigs {

    label 'sm'

    input:
        tuple strain, path(bam), path(bai)

    output:
        path("contigs.txt")

    """
        samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    """

}

/*===========================================
~ ~ ~ > * Run GATK4 HaplotypeCaller * < ~ ~ ~ 
===========================================*/

process call_variants_individual {

    label 'md'

    tag { "${strain}:${region}" }

    input:
        tuple strain, path(bam), path(bai), val(region)

    output:
        tuple strain, path("${region}.g.vcf.gz")

    """
        gatk HaplotypeCaller --java-options "-Xmx${task.memory.toGiga()}g -Xms1g -XX:ConcGCThreads=${task.cpus}" \\
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
            -R ${reference} \\
            -I ${bam} \\
            -L ${region} \\
            -O ${region}.g.vcf   
        bcftools view -O z ${region}.g.vcf > ${region}.g.vcf.gz
        rm ${region}.g.vcf
    """
}

/*
=============================================
~ ~ ~ > *  Merge Sample gVCF files  * < ~ ~ ~ 
=============================================
*/

process concat_strain_gvcfs {

    label 'sm'
    tag { "${strain}" }

    input:
        tuple strain, path("*"), path(contigs)

    output:
        tuple path("${strain}.g.vcf.gz"), path("${strain}.g.vcf.gz.tbi")

    """
        awk '{ print \$0 ".g.vcf.gz" }' ${contigs} > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > ${strain}.g.vcf.gz
        bcftools index --tbi ${strain}.g.vcf.gz
    """
}

 process import_genomics_db {

    tag { "${contig}" }
    label 'lg'

    input:
        tuple path(vcfs), contig, path(sample_map)

    output:
        tuple val(contig), file("${contig}.db")

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()-3}g -XX:ConcGCThreads=${task.cpus}" \\
            GenomicsDBImport \\
            --genomicsdb-workspace-path ${contig}.db \\
            --batch-size 16 \\
            -L ${contig} \\
            --sample-name-map ${sample_map} \\
            --reader-threads ${task.cpus}
    """
}

/*==================================
~ ~ ~ > *  Genotype gVCFs  * < ~ ~ ~ 
==================================*/

process genotype_cohort_gvcf_db {
    // Heterozygous polarization is also performed here.

    tag { "${contig}" }
    label 'lg'

    input:
        tuple val(contig), file("${contig}.db")

    output:
        tuple val(contig), file("${contig}_cohort.bcf"), file("${contig}_cohort.bcf.csi")

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            GenotypeGVCFs \\
            -R ${reference} \\
            -V gendb://${contig}.db \\
            -G StandardAnnotation \\
            -G AS_StandardAnnotation \\
            -G StandardHCAnnotation \\
            -L ${contig} \\
           -O ${contig}_cohort.vcf

        bcftools view -O z ${contig}_cohort.vcf | \
        het_polarization > ${contig}_cohort.bcf
        bcftools index ${contig}_cohort.bcf
    """
}

/*=====================================================
~ ~ ~ > *  SnpEff Annotate Cohort Chrom VCFs  * < ~ ~ ~  
=====================================================*/

process annotate_vcf {

    label 'lg' 
    
    tag { contig }

    cache 'lenient'

    input:
        tuple val(contig), file("${contig}.bcf"), file("${contig}.bcf.csi")

    output:
        tuple path("${contig}.annotated.vcf.gz"), path("${contig}.annotated.vcf.gz.tbi"), emit: 'anno_vcf'
        path "${contig}.${date}.snpeff.csv", emit: 'snpeff_csv'


    """
      bcftools view -O v --threads=${task.cpus-1} ${contig}.bcf | \\
      snpEff eff -csvStats ${contig}.${date}.snpeff.csv \\
                 -no-downstream \\
                 -no-intergenic \\
                 -no-upstream \\
                 -nodownload \\
      -dataDir ${params.snpeff_dir} \\
      -config ${params.snpeff_dir}/snpEff.config \\
      ${params.snpeff_reference} | \\
      bcftools view --threads=${task.cpus-1} -O z > ${contig}.annotated.vcf.gz
      bcftools index --tbi ${contig}.annotated.vcf.gz
    """

}

/*===============================================
~ ~ ~ > *   Concatenate Annotated VCFs  * < ~ ~ ~
===============================================*/

process concatenate_vcf {

    label 'lg'

    input:
        tuple path("*"), path("contigs.txt")

    output:
        tuple path("WI.annotated.vcf.gz"), path("WI.annotated.vcf.gz.tbi"), emit: 'vcf'

    """
        awk '{ print \$0 ".annotated.vcf.gz" }' contigs.txt > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > WI.annotated.vcf.gz
        bcftools index --tbi WI.annotated.vcf.gz
    """
}

/*===========================================
~ ~ ~ > *   Apply SNV Soft Filters  * < ~ ~ ~
===========================================*/

process soft_filter {

    publishDir "${params.output}/variation", mode: 'copy'

    label 'lg'

    input:
        tuple path(vcf), path(vcf_index)

    output:
        tuple path("WI.${date}.soft-filter.vcf.gz"), path("WI.${date}.soft-filter.vcf.gz.csi"), emit: soft_filter_vcf
        path "WI.${date}.soft-filter.vcf.gz.tbi"
        path "WI.${date}.soft-filter.stats.txt", emit: 'soft_vcf_stats'
        path "WI.${date}.soft-filter.filter_stats.txt"
    
    /*
        ad_dp is a binary that adds the ad_dp filter. Do not remove.
        het_polarization polarizes het-variants to REF or ALT
    */
    """
        function cleanup {
            rm ad_dp.filtered.vcf.gz
        }
        trap cleanup EXIT

        gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            VariantFiltration \\
            -R ${reference} \\
            --variant ${vcf} \\
            --genotype-filter-expression "DP < ${params.min_depth}"    --genotype-filter-name "DP_min_depth" \\
            --filter-expression "QUAL < ${params.qual}"                --filter-name "QUAL_quality" \\
            --filter-expression "FS > ${params.fisherstrand}"          --filter-name "FS_fisherstrand" \\
            --filter-expression "QD < ${params.quality_by_depth}"      --filter-name "QD_quality_by_depth" \\
            --filter-expression "SOR > ${params.strand_odds_ratio}"    --filter-name "SOR_strand_odds_ratio" \\
            --genotype-filter-expression "isHet == 1"                  --genotype-filter-name "is_het" \\
            -O /dev/stdout | \\
            ad_dp
        
        # ad_dp filter
        bcftools index --tbi ad_dp.filtered.vcf.gz
        
        # Apply high missing and high heterozygosity filters
        bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING  <= ${params.high_missing}' ad_dp.filtered.vcf.gz |\\
        bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > WI.${date}.soft-filter.vcf.gz

        bcftools index WI.${date}.soft-filter.vcf.gz
        bcftools index --tbi WI.${date}.soft-filter.vcf.gz
        bcftools stats --threads ${task.cpus} \\
                       -s- --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt

        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.filter_stats.txt

    """
}

/*==============================================
~ ~ ~ > *   Apply Hard Filters on VCF  * < ~ ~ ~
==============================================*/

process hard_filter {
    /*
        !! Important
        CSQ Annotations take place after hard filtering because they create a haplotype-based prediction.
        THerefore, it is essential that poor-quality variants are removed.
    */


    label 'lg'

    publishDir "${params.output}/variation", mode: 'copy'

    input:
        tuple path(vcf), path(vcf_index), path(contigs)

    output:
        tuple path("WI.${date}.hard-filter.vcf.gz"), path("WI.${date}.hard-filter.vcf.gz.csi"), emit: 'vcf'
        path "WI.${date}.hard-filter.vcf.gz.tbi"
        path "WI.${date}.hard-filter.stats.txt", emit: 'hard_vcf_stats'


    """
        function cleanup {
            # cleanup files on completion
            rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
        }
        trap cleanup EXIT
        # Generate hard-filtered VCF
        function generate_hard_filter {
            bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT != "PASS"' |\\
            bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
            bcftools view -O v --min-af 0.0000000000001 --max-af 0.999999999999 |\\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
        }

        export -f generate_hard_filter

        parallel --verbose generate_hard_filter :::: ${contigs}

        # Remove transposons from GFF3 as they throw errors with CSQ
        gzip -dc ${file(params.csq_gff, checkIfExists: true)} | grep -v 'transposon' | bgzip > csq.gff.gz
        tabix -p gff csq.gff.gz


        awk '{ print \$0 ".vcf.gz" }' ${contigs} > contig_set.tsv
        bcftools concat  -O u --file-list contig_set.tsv | \\
        bcftools csq -O z --fasta-ref ${reference} \\
                     --gff-annot csq.gff.gz \\
                     --phase a > WI.${date}.hard-filter.vcf.gz

        bcftools index WI.${date}.hard-filter.vcf.gz
        bcftools index --tbi WI.${date}.hard-filter.vcf.gz
        bcftools stats -s- --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt
    """
}

process strain_list {

    input:
        tuple path(vcf), path(vcf_index)
    
    output:
        file("samples.txt")

    """
        bcftools query --list-samples ${vcf} > samples.txt
    """
}


process generate_strain_tsv_vcf {
    // Generate a single TSV and VCF for every strain.

    tag { strain }

    publishDir "${params.output}/strain/vcf", mode: 'copy', pattern: "*.vcf.gz*"
    publishDir "${params.output}/strain/tsv", mode: 'copy', pattern: "*.tsv.gz*"

    input:
        tuple val(strain), path(vcf), file(vcf_index)

    output:
        tuple path("${strain}.${date}.vcf.gz"),  path("${strain}.${date}.vcf.gz.tbi")
        tuple path("${strain}.${date}.vcf.gz"),  path("${strain}.${date}.vcf.gz.csi")
        tuple path("${strain}.${date}.tsv.gz"),  path("${strain}.${date}.tsv.gz.tbi")

    """
        # Generate VCF
        bcftools view -O z --samples ${strain} \\
                           --exclude-uncalled \\
                           ${vcf}  > ${strain}.${date}.vcf.gz
        bcftools index ${strain}.${date}.vcf.gz
        bcftools index --tbi ${strain}.${date}.vcf.gz

        # Generate TSV
        {
            echo 'CHROM\\tPOS\\tREF\\tALT\\tFILTER\\tFT\\tGT';
            bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\t%FILTER\\t%FT\\t%TGT]\\n' --samples ${strain} ${vcf};
        } > ${strain}.${date}.tsv
        bgzip ${strain}.${date}.tsv
        tabix -S 1 -s 1 -b 2 -e 2 ${strain}.${date}.tsv.gz
    """

}

process generate_severity_tracks {
    /*
        The severity tracks are bedfiles with annotations for
        LOW
        MODERATE
        HIGH
        MODIFIER
        variants as annotated with SNPEff

        They are used on the CeNDR browser.
    */

    publishDir "${params.output}/tracks", mode: 'copy'

    tag { severity }

    input:
        tuple path("in.vcf.gz"), path("in.vcf.gz.csi"), val(severity)
    output:
        set file("${date}.${severity}.bed.gz"), file("${date}.${severity}.bed.gz.tbi")

    """
        bcftools view --apply-filters PASS in.vcf.gz | \
        grep ${severity} | \
        awk '\$0 !~ "^#" { print \$1 "\\t" (\$2 - 1) "\\t" (\$2)  "\\t" \$1 ":" \$2 "\\t0\\t+"  "\\t" \$2 - 1 "\\t" \$2 "\\t0\\t1\\t1\\t0" }' | \\
        bgzip  > ${date}.${severity}.bed.gz
        tabix -p bed ${date}.${severity}.bed.gz
        fsize=\$(zcat ${date}.${severity}.bed.gz | wc -c)
        if [ \${fsize} -lt 2000 ]; then
            exit 1
        fi;
    """
}

process multiqc_report {

    publishDir "${params.output}/report", mode: 'copy'

    input:
        path("*")

    output:
        file("multiqc_data/*.json")
        file("multiqc.html")

    """
        multiqc -k json --filename multiqc.html .
    """

}

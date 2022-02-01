#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Variant Calling Pipeline
    Authors:
    - Stefan Zdraljevic
    - Daniel Cook <danielecook@gmail.com>
    - Dan Lu
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
params.bam_location = "/projects/b1059/data/${params.species}/WI/alignments/" // Use this to specify the directory for bams
params.mito_name = "MtDNA" // Name of contig to skip het polarization


// Check that reference exists
// params.reference = ""
//reference = file(params.reference, checkIfExists: true)

// Debug
if (params.debug) {
    params.output = "release-debug"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bam_location = "${workflow.projectDir}" // Use this to specify the directory for bams
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "WI-${date}"
    params.sample_sheet = "${workflow.projectDir}/sample_sheet.tsv"
}

// set default project and ws build for species
if(params.species == "c_elegans") {
    params.project="PRJNA13758"
    params.ws_build="WS283"
} else if(params.species == "c_briggsae") {
    params.project="QX1410_nanopore"
    params.ws_build="Feb2020"
} else if(params.species == "c_tropicalis") {
    params.project="NIC58_nanopore"
    params.ws_build="June2021"
}

// check reference
if(params.species == "c_elegans" | params.species == "c_briggsae" | params.species == "c_tropicalis") {
    params.reference = "/projects/b1059/data/${params.species}/genomes/${params.project}/${params.ws_build}/${params.species}.${params.project}.${params.ws_build}.genome.fa.gz"
} else if (params.species == null) {
    if (params.reference == null) {
        if (params.help) {
        } else { 
        println """

        Please specify a species: c_elegans c_brigssae c_tropicalis with option --species, or a ref genome with --reference"

        """
        exit 1
        }
    }
}

reference = file(params.reference, checkIfExists: true)


// Variant Filtering
params.min_depth = 5
params.qual = 30.0
params.strand_odds_ratio = 5.0
params.quality_by_depth = 20.0
params.fisherstrand = 100.0
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

To run the pipeline:

nextflow main.nf --debug
nextflow main.nf -profile debug
nextflow main.nf --sample_sheet=/path/sample_sheet_GATK.tsv --bam_location=/projects/b1059/workflows/alignment-nf/

    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    ${params.debug}
    --output                   Release Directory                     ${params.output}
    --sample_sheet             Sample sheet                          ${params.sample_sheet}
    --bam_location             Directory of bam files                ${params.bam_location}
    --reference                Reference Genome                      ${params.reference}
    --mito_name                Contig not to polarize hetero sites   ${params.mito_name}
    --username                                                       ${"whoami".execute().in.text}

    Reference Genome
    --------------- 
    --reference_base           Location of ref genomes               ${params.reference_base}
    --species/project/build    These 4 params form --reference       ${params.species} / ${params.project} / ${params.ws_build}

    Variant Filters         
    ---------------           
    --min_depth                Minimum variant depth                 ${params.min_depth}
    --qual                     Variant QUAL score                    ${params.qual}
    --strand_odds_ratio        SOR_strand_odds_ratio                 ${params.strand_odds_ratio} 
    --quality_by_depth         QD_quality_by_depth                   ${params.quality_by_depth} 
    --fisherstrand             FS_fisher_strand                      ${params.fisherstrand}
    --high_missing             Max % missing genotypes               ${params.high_missing}
    --high_heterozygosity      Max % max heterozygosity              ${params.high_heterozygosity}

---
"""
out
}

log.info(log_summary())

if (params.help) {
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
                                     file("${row.bam}", checkIfExists: true),
                                     file("${row.bam}.bai", checkIfExists: true)] }


workflow {

    // Generate a summary of the current run
    summary(Channel.from("run").combine(Channel.from("${params.sample_sheet}")))

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
                           .collect() // .toList() might be causing this process to always repeat, try switching to collect (https://gitter.im/nextflow-io/nextflow/archives/2018/10/08)
                           .map { [it] }
                           .combine(contigs)
                           .combine(sample_map) | \
                        import_genomics_db | \
                        genotype_cohort_gvcf_db


    // Combine VCF and filter
    genotype_cohort_gvcf_db.out.collect() | concatenate_vcf
    concatenate_vcf.out.vcf | soft_filter
    soft_filter.out.soft_filter_vcf.combine(get_contigs.out) | hard_filter


    // MultiQC Report
    soft_filter.out.soft_vcf_stats.concat(
        hard_filter.out.hard_vcf_stats
    ).collect() | multiqc_report

}


process summary {
    // Generates a summary of the run for the release directory.
    
    executor 'local'

    publishDir "${params.output}", mode: 'copy'
    
    input:
        tuple val(run), path("sample_sheet")

    output:
        path("sample_sheet.tsv")
        path("summary.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        cat ${sample_sheet} > sample_sheet.tsv
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
            --do-not-run-physical-phasing \\
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
    label 'xl'

    input:
        tuple path(vcfs), contig, path(sample_map)

    output:
        tuple val(contig), file("${contig}.db")

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()-3}g -Xms${task.memory.toGiga()-4}g -XX:ConcGCThreads=${task.cpus}" \\
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
    label 'xl'

    input:
        tuple val(contig), file("${contig}.db")

    output:
        tuple file("${contig}_cohort_pol.vcf.gz"), file("${contig}_cohort_pol.vcf.gz.csi")


    /*
        het_polarization polarizes het-variants to REF or ALT (except for mitochondria)
    */

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

        if [ "${contig}" == "${params.mito_name}" ]
        then
            bcftools view -O b ${contig}_cohort.vcf > ${contig}_cohort.bcf
            bcftools index ${contig}_cohort.bcf

        else
        
            bcftools view -O z ${contig}_cohort.vcf | \\
            het_polarization > ${contig}_cohort.bcf
            bcftools index ${contig}_cohort.bcf
        
        fi

        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} ${contig}_cohort.bcf | \\
        vcffixup - | \\
        bcftools view -O z > ${contig}_cohort_pol.vcf.gz
        bcftools index ${contig}_cohort_pol.vcf.gz
    """
}


/*===============================================
~ ~ ~ > *   Concatenate VCFs  * < ~ ~ ~
===============================================*/

process concatenate_vcf {

    label 'lg'

    input: 
      path("*")

    output:
        tuple path("WI.raw.vcf.gz"), path("WI.raw.vcf.gz.tbi"), emit: 'vcf'

    """
        ls *_cohort_pol.vcf.gz > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > WI.raw.vcf.gz
        bcftools index --tbi WI.raw.vcf.gz
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
    

    """
        function cleanup {
            rm out.vcf.gz
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
            -O out.vcf
        
        bgzip out.vcf
        bcftools index --tbi out.vcf.gz
        
        # Apply high missing and high heterozygosity filters
        bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING  <= ${params.high_missing}' out.vcf.gz |\\
        bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > WI.${date}.soft-filter.vcf.gz

        bcftools index WI.${date}.soft-filter.vcf.gz
        bcftools index --tbi WI.${date}.soft-filter.vcf.gz
        bcftools stats --threads ${task.cpus} \\
                       -s- --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt

        {
            echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
            bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' WI.${date}.soft-filter.vcf.gz;
        }     > WI.${date}.soft-filter.filter_stats.txt

    """
}

/*==============================================
~ ~ ~ > *   Apply Hard Filters on VCF  * < ~ ~ ~
==============================================*/

process hard_filter {


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

        if [ \${1} == "${params.mito_name}" ]

        then
            bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' |\\
            bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
            bcftools view -O v --min-af 0.000001 --max-af 0.999999 |\\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz

        else
            bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' |\\
            bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
            bcftools view -O v --min-af 0.000001 --max-af 0.999999 |\\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
        fi

        }

        export -f generate_hard_filter

        parallel --verbose generate_hard_filter :::: ${contigs}


        awk '{ print \$0 ".vcf.gz" }' ${contigs} > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > WI.${date}.hard-filter.vcf.gz

        bcftools index WI.${date}.hard-filter.vcf.gz
        bcftools index --tbi WI.${date}.hard-filter.vcf.gz
        bcftools stats -s- --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt
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

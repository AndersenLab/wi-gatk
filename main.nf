#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Variant Calling Pipeline
    Authors:
    - Stefan Zdraljevic
    - Daniel Cook <danielecook@gmail.com>
    - Dan Lu
    - Mike Sauria <mike.sauria@jhu.edu>
*/

include { CREATE_BAM_SAMPLE_SHEET                           } from "./modules/sample_sheet/create_bam_sample_sheet/main"
include { SAMTOOLS_GET_CONTIGS                              } from "./modules/samtools/get_contigs/main"
include { GATK_HAPLOTYPECALLER                              } from "./modules/gatk/haplotypecaller/main"
include { BCFTOOLS_CONCAT_VCFS as BCFTOOLS_CONCAT_GVCFS     } from "./modules/bcftools/concat_vcfs/main"
include { BCFTOOLS_CONCAT_VCFS as BCFTOOLS_CONCAT_SOFT_VCFS } from "./modules/bcftools/concat_vcfs/main"
include { BCFTOOLS_CONCAT_VCFS as BCFTOOLS_CONCAT_HARD_VCFS } from "./modules/bcftools/concat_vcfs/main"
include { BCFTOOLS_CONCAT_OVERLAPPING_VCFS as BCFTOOLS_CONCAT_OVERLAPPING_SOFT } from "./modules/bcftools/concat_overlapping_vcfs/main"
include { BCFTOOLS_CONCAT_OVERLAPPING_VCFS as BCFTOOLS_CONCAT_OVERLAPPING_HARD } from "./modules/bcftools/concat_overlapping_vcfs/main"
include { GATK_GENOMICSDBIMPORT                             } from "./modules/gatk/genomicsdbimport/main"
include { GATK_GENOTYPEGVCFS                                } from "./modules/gatk/genotypegvcfs/main"
include { LOCAL_HETPOLARIZATION                             } from "./modules/local/hetpolarization/main"
include { GATK_SOFT_FILTER                                  } from "./modules/gatk/soft_filter/main"
include { LOCAL_FILTER                                      } from "./modules/local/filter/main"
include { LOCAL_VCF_STATS                                   } from "./modules/local/vcf_stats/main"
include { MULTIQC_REPORT                                    } from "./modules/multiqc/report/main"
include { LOCAL_REPORT                                      } from "./modules/local/report/main"

// Needed to publish results
nextflow.preview.output = true

date = new Date().format( 'yyyyMMdd' )

// Debug
if (params.debug) {
    species = "c_elegans"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.txt"
    bam_folder = "${workflow.projectDir}/test_data/bams"
    gvcf_folder = "${workflow.projectDir}/test_data/gVCFs"
}

if (params.help == false & params.debug == false) {
    if (params.species == null) {
        println """
        Please specify a species with option --species
        """
        exit 1
    } else {
        species = params.species
    }
    if (params.sample_sheet == null) {
        println """
        Please specify a sample sheet with option --sample_sheet
        """
        exit 1
    }
    if (species != "c_elegans" & species != "c_briggsae" & species != "c_tropicalis"){
        if (params.reference == null) {
            println """
            When using a species other than C. elegans, C. briggsae, or C. tropicalis,
            a reference genome must be specified with option --reference
            """
            exit 1
        }
        if (params.bam_location == null) {
            println """
            When using a species other than C. elegans, C. briggsae, or C. tropicalis,
            a BAM directory must be specified with option --bam_location
            """
            exit 1
        }
        if (params.gvcf_location == null) {
            println """
            When using a species other than C. elegans, C. briggsae, or C. tropicalis,
            a gVCF directory must be specified with option --gvcf_location
            """
            exit 1
        }
    } else {
        if(params.bam_location != null) {
            bam_folder = "${params.bam_location}"
        } else {
            bam_folder = "${params.data_path}/${species}/WI/alignments/"
        }
        if(params.gvcf_location != null) {
            gvcf_folder = "${params.gvcf_location}"
        } else {
            gvcf_folder = "${params.data_path}/${species}/WI/gVCFs/"
        }
    }
} else if (params.debug == false) {
    if(params.bam_location != null) {
        bam_folder = "${params.bam_location}"
    } else if (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis"){
        bam_folder = "${params.data_path}/${species}/WI/alignments/"
    }
    if(params.gvcf_location != null) {
        gvcf_folder = "${params.gvcf_location}"
    } else if (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis"){
        gvcf_folder = "${params.data_path}/${species}/WI/gVCFs/"
    }
}

// set default project and ws build for species
if(species == "c_elegans") {
    params.project="PRJNA13758"
    params.ws_build="WS283"
} else if(species == "c_briggsae") {
    params.project="QX1410_nanopore"
    params.ws_build="Feb2020"
} else if(species == "c_tropicalis") {
    params.project="NIC58_nanopore"
    params.ws_build="June2021"
}

// check reference
if (params.reference == null){
    if(params.data_path != null && (species == "c_elegans" | species == "c_briggsae" | species == "c_tropicalis")) {
        reference = "${params.data_path}/${species}/genomes/${params.project}/${params.ws_build}/${species}.${params.project}.${params.ws_build}.genome.fa.gz"
    } else if (params.help) {
        reference = null
    } else { 
        println """
        Please specify c_elegans, c_brigssae, or c_tropicalis as the species with option --species
        or a reference genome with --reference
        """
        exit 1
    }
} else {
    reference = params.reference
}

if (reference != null){
    ref_base = reference.take(reference.take(reference.lastIndexOf('.') - 1).lastIndexOf("."))
    reference_fa = "${reference}"
    reference_index = "${reference}.fai"
    reference_dict = "${ref_base}.dict"
    reference_gzi = "${reference}.gzi"
}


def log_summary() {

    out =  '''
 _______ _______ _______ __  __      _______ _______ 
|     __|   _   |_     _|  |/  |    |    |  |    ___|
|    |  |       | |   | |     <     |       |    ___|
|_______|___|___| |___| |__|\\__|    |__|____|___|    
                                              
'''

out += """

To run the pipeline:

nextflow main.nf --help
nextflow main.nf --debug
nextflow main.nf --sample_sheet=/path/sample_sheet.txt --species c_elegans --bam_location=/path/to/bams --gvcf_location=/path/to/gvcfs

    parameters                 description                           Set/Default
    ==========                 ===========                           ========================
    --debug                    Use --debug to indicate debug mode    ${params.debug}
    --species                  Species to call variants from         ${species}
    --sample_sheet             Sample sheet                          ${params.sample_sheet}
    --bam_location             Directory of BAM files                ${bam_folder}
    --gvcf_location            Directory of gVCF files               ${gvcf_folder}
    --mito_name                Contig not to polarize hetero sites   ${params.mito_name}
    --partition                Partition size in bp for subsetting   ${params.partition}
    --gvcf_only                Create sample gVCFs and stop          ${params.gvcf_only}
    --username                                                       ${"whoami".execute().in.text}

    Reference Genome
    --------------- 
    --reference                The fa.gz reference file to use       ${reference}

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

if (params.help == true) {
    exit 1
}

now = new Date()
timestamp = now.format("yyyyMMdd-HH-mm-ss")
log.info("Started running ${now}")

workflow {
    main:
    ch_versions = Channel.empty()

    // Read sample sheet
    sample_sheet_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)
        .ifEmpty { exit 1, "sample sheet not found" }

    // Make channel for reference
    reference_ch = Channel.of( ["id": species] )
        .combine(Channel.fromPath(reference_fa, checkIfExists: true))
        .combine(Channel.fromPath(reference_index, checkIfExists: true))
        .combine(Channel.fromPath(reference_dict, checkIfExists: true))
        .combine(Channel.fromPath(reference_gzi, checkIfExists: true))
        .first()

    // Make channel of filter parameters
    filter_ch = Channel.of( [min_depth: params.min_depth, qual: params.qual, fisherstrand: params.fisherstrand,
                             quality_by_depth: params.quality_by_depth, strand_odds_ratio: params.strand_odds_ratio,
                             high_missing: params.high_missing, high_heterozygosity: params.high_heterozygosity] )
        .first()

    // Determine which samples are missing gVCF files and create them
    CREATE_BAM_SAMPLE_SHEET( sample_sheet_ch,
                             Channel.fromPath(bam_folder, checkIfExists: true),
                             Channel.fromPath(gvcf_folder, checkIfExists: true) )

    // Create channel of bams and indices
    bam_ch = CREATE_BAM_SAMPLE_SHEET.out.bam
        .splitCsv( strip: true )
        .map{ it: [ [id: it[0]], file("${bam_folder}/${it[0]}.bam"), file("${bam_folder}/${it[0]}.bam.bai") ] }

    // Get contigs from bam file
    SAMTOOLS_GET_CONTIGS( sample_sheet_ch.splitCsv( strip: true )
        .first()
        .map{ it: [ [id: it[0]], file("${bam_folder}/${it[0]}.bam"), file("${bam_folder}/${it[0]}.bam.bai") ] },
        params.partition )
    ch_versions = ch_versions.mix(SAMTOOLS_GET_CONTIGS.out.versions)

    partitions_ch = SAMTOOLS_GET_CONTIGS.out.partitions.splitCsv( strip: true, sep: "\t" )
        .map{ it: [interval: "${it[0]}:${it[2]}-${it[3]}", label:  "${it[0]}_${it[2]}_${it[3]}", contig: it[0], start: it[2], end: it[3], size: it[1]] }

    ///////////////////////////////////////////////////////
    // This section is for samples missing premade gVCFs //
    ///////////////////////////////////////////////////////

    // Create sample by contig channel
    bam_contig_ch = bam_ch
        .combine(SAMTOOLS_GET_CONTIGS.out.contigs.splitCsv( strip: true ))
        .map{ it: [[id: it[0].id, contig: it[3]], it[1], it[2]] }

    // Call variants in each sample/contig
    GATK_HAPLOTYPECALLER( bam_contig_ch,
                          reference_ch )
    ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions)

    // Group contig variant calls by sample
    gvcf_contig_ch = GATK_HAPLOTYPECALLER.out.vcf
        .map{ it: [ it[0].id, it[1] ] }
        .groupTuple()
        .map{ it: [ [id: it[0]], it[1] ] }

    // Combine contigs for each sample into single gVCF samples
    BCFTOOLS_CONCAT_GVCFS( gvcf_contig_ch,
                           SAMTOOLS_GET_CONTIGS.out.contigs )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_GVCFS.out.versions)

    // Collect concatenated gVCFs
    new_gvcf_ch = BCFTOOLS_CONCAT_GVCFS.out.vcf

    if (params.gvcf_only == false){
        ///////////////////////////////////////////////////////
        // This section uses all gVCFs for joint genotyping  //
        ///////////////////////////////////////////////////////

        // Create gVCF sample map
        premade_gvcf_ch = CREATE_BAM_SAMPLE_SHEET.out.gvcf
            .splitCsv( strip: true )
            .map{ it: [ [id: it[0]], file("${gvcf_folder}/${it[0]}.g.vcf.gz"), file("${gvcf_folder}/${it[0]}.g.vcf.gz.tbi") ] }

        sample_map_ch = bam_ch
            .map{ it: "${it[0].id}\tnew_gvcfs/${it[0].id}.g.vcf.gz" }
            .concat(premade_gvcf_ch
                .map{ it: "${it[0].id}\tpremade_gvcfs/${it[0].id}.g.vcf.gz" })
            .collectFile(name: "sample_map.tsv", newLine: true)
            .first()

        flat_new_gvcf_ch = new_gvcf_ch
            .map{ it: it[1] }
            .collect()
            .concat(new_gvcf_ch
                .map{ it: it[2] }
                .collect())
            .collect()
            .map{ it: [it] }

        // Combine allgVCFs with partitions
        // all_gvcf_contig_ch = partitions_ch
        //     .combine(all_gvcfs_ch)

        // Create a genomics db for each contig
        GATK_GENOMICSDBIMPORT( Channel.fromPath(gvcf_folder),
                               flat_new_gvcf_ch,
                               partitions_ch,
                               sample_map_ch )
        ch_versions = ch_versions.mix(GATK_GENOMICSDBIMPORT.out.versions)

        // Genotype cohort contig db
        GATK_GENOTYPEGVCFS( GATK_GENOMICSDBIMPORT.out.db,
                            reference_ch )
        ch_versions = ch_versions.mix(GATK_GENOTYPEGVCFS.out.versions)

        // Concatenate VCFs by contig and perform het polarization
        LOCAL_HETPOLARIZATION( GATK_GENOTYPEGVCFS.out.vcf,
                               params.mito_name )
        ch_versions = ch_versions.mix(LOCAL_HETPOLARIZATION.out.versions)

        // Mark filtered genotypes/variants with GATK
        GATK_SOFT_FILTER( filter_ch,
                          LOCAL_HETPOLARIZATION.out.vcf,
                          reference_ch )
        ch_versions = ch_versions.mix(GATK_SOFT_FILTER.out.versions)

        // Mark filtered genotypes/variants with bcftools and remove them
        LOCAL_FILTER( filter_ch,
                      GATK_SOFT_FILTER.out.vcf,
                      reference_ch,
                      params.mito_name )
        ch_versions = ch_versions.mix(LOCAL_FILTER.out.versions)

        // Collect cohort VCFs by contig for concatenation
        soft_cohort_contig_ch = LOCAL_FILTER.out.soft
            .map{ it: [it[0].contig, it[1], it[2]] } 
            .groupTuple()
            .map{ it: [[contig: it[0], id: "soft"], it[1], it[2]] }

        hard_cohort_contig_ch = LOCAL_FILTER.out.hard
            .map{ it: [it[0].contig, it[1], it[2]] } 
            .groupTuple()
            .map{ it: [[contig: it[0], id: "hard"], it[1], it[2]] }

        // Concat partitioned filtered vcfs
        BCFTOOLS_CONCAT_OVERLAPPING_SOFT( soft_cohort_contig_ch,
                                          SAMTOOLS_GET_CONTIGS.out.partitions )

        BCFTOOLS_CONCAT_OVERLAPPING_HARD( hard_cohort_contig_ch,
                                          SAMTOOLS_GET_CONTIGS.out.partitions )

        // Collect genotyped concatenated contigs together
        soft_vcf_ch = BCFTOOLS_CONCAT_OVERLAPPING_SOFT.out.vcf
            .map{ it: it[1] }
            .toSortedList()
            .map{ it: [[id: "soft"], it] }
            .view()

        hard_vcf_ch = BCFTOOLS_CONCAT_OVERLAPPING_HARD.out.vcf
            .map{ it: it[1] }
            .toSortedList()
            .map{ it: [[id: "hard"], it] }
            .view()

        // Combine contigs for soft-filtered variant calls
        BCFTOOLS_CONCAT_SOFT_VCFS( soft_vcf_ch,
                                  SAMTOOLS_GET_CONTIGS.out.contigs )
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_SOFT_VCFS.out.versions)

        // Combine contigs for hard-filtered variant calls
        BCFTOOLS_CONCAT_HARD_VCFS( hard_vcf_ch,
                                  SAMTOOLS_GET_CONTIGS.out.contigs )
        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT_HARD_VCFS.out.versions)

        // Get filtering stats for filtered VCFs
        LOCAL_VCF_STATS( BCFTOOLS_CONCAT_SOFT_VCFS.out.vcf,
                         BCFTOOLS_CONCAT_HARD_VCFS.out.vcf )
        ch_versions = ch_versions.mix(LOCAL_VCF_STATS.out.versions)

        // Create multiqc report
        MULTIQC_REPORT( LOCAL_VCF_STATS.out.soft_stats
                            .concat(LOCAL_VCF_STATS.out.hard_stats)
                            .collect() )
        ch_versions = ch_versions.mix(MULTIQC_REPORT.out.versions)

        // Create GATK report
        LOCAL_REPORT( MULTIQC_REPORT.out.report,
                      LOCAL_VCF_STATS.out.soft_filter_stats,
                      LOCAL_VCF_STATS.out.soft_stats,
                      LOCAL_VCF_STATS.out.hard_stats,
                    "${workflow.projectDir}/bin/gatk_report.Rmd",
                    params.timestamp )
        ch_versions = ch_versions.mix(LOCAL_REPORT.out.versions)

        ch_soft_vcf       = BCFTOOLS_CONCAT_SOFT_VCFS.out.vcf
        ch_hard_vcf       = BCFTOOLS_CONCAT_HARD_VCFS.out.vcf
        ch_filter_stats   = LOCAL_VCF_STATS.out.soft_filter_stats
        ch_soft_stats     = LOCAL_VCF_STATS.out.soft_stats
        ch_hard_stats     = LOCAL_VCF_STATS.out.hard_stats
        ch_multiqc_report = MULTIQC_REPORT.out.report
        ch_multiqc_json   = MULTIQC_REPORT.out.json
        ch_local_report   = LOCAL_REPORT.out.html

    } else {
        ch_soft_vcf       = Channel.empty()
        ch_hard_vcf       = Channel.empty()
        ch_filter_stats   = Channel.empty()
        ch_soft_stats     = Channel.empty()
        ch_hard_stats     = Channel.empty()
        ch_multiqc_report = Channel.empty()
        ch_multiqc_json   = Channel.empty()
        ch_local_report   = Channel.empty()
    }

    // Collate and save software versions
        ch_versions
            .collectFile(name: 'workflow_software_versions.txt', sort: true, newLine: true)
            .set { ch_collated_versions }
        
    publish:
    BCFTOOLS_CONCAT_GVCFS.out.vcf >> "gVCFs"
    ch_soft_vcf                   >> "variation"
    ch_hard_vcf                   >> "variation"
    ch_filter_stats               >> "variation"
    ch_soft_stats                 >> "variation"
    ch_hard_stats                 >> "variation"
    ch_multiqc_report             >> "report"
    ch_multiqc_json               >> "report"
    ch_local_report               >> "report"
    ch_collated_versions          >> "."
}

// Current bug that publish doesn't work without an output closure
output {
    "gVCFs" {
        mode "copy"
    }
    "variation" {
        mode "copy"
    }
    "report" {
        mode "copy"
    }
    "." {
        mode "copy"
    }
}

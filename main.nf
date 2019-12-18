#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Variant Calling Pipeline
    Authors:
    - Stefan Zdraljevic
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2

// For now, this pipeline requires NXF_VER 19.09.0
// Prefix this version when running
// e.g.
// NXF_VER=19.09.0-edge nextflow run ...
assert System.getenv("NXF_VER") == "19.09.0-edge"

/*
    Params
*/

date = new Date().format( 'yyyyMMdd' )
params.debug = false
params.email = ""
params.reference = "${workflow.projectDir}/WS245/WS245.fa.gz"
reference_uncompressed = file(params.reference.replace(".gz", ""), checkExists: true)
parse_conda_software = file("${workflow.projectDir}/scripts/parse_conda_software.awk")

// Debug
if (params.debug.toString() == "true") {
    params.output = "release-${date}-debug"
    params.sample_sheet = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bamdir = "${workflow.projectDir}/test_data"
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "alignment-${date}"
    params.sample_sheet = "sample_sheet.tsv"
}

// defaults
params.bamdir                 = null
params.out_base               = null
params.gff                    = null
params.strains                = false
params.annotation_reference   = null

params.out                    = "${date}-${params.out_base}"

def log_summary() {

    out =  '''
 _______ _______ _______ __  __      _______ _______ 
|     __|   _   |_     _|  |/  |    |    |  |    ___|
|    |  |       | |   | |     <     |       |    ___|
|_______|___|___| |___| |__|\\__|    |__|____|___|    
                                              
'''
out += """
    ----------------------------------------------------------------
                          USAGE                                     
    ----------------------------------------------------------------
    
    nextflow ID_diverged-regions.nf --out_base Analysis
    
    Mandatory arguments:
    --out_base             String                Name of folder to output results
    --strains              String                Name of folder to output results 
    --bamdir               String                Name of folder where bam files are
    
    --------------------------------------------------------
    Optional arguments:
    Information describing the stucture of the input files can be located in input_files/README.txt
    
    --cores                INTEGER              Number of cpu to use (default=2)
    --annotation_reference STRING               Wormbase build for SnpEff annotations
    --email                STRING               email address for job notifications
    
    Flags:
    --help                                      Display this message
    
    --------------------------------------------------------
    
     Required software packages to be in users path
    BCFtools               v1.9
    GATK4                  v4.1
    --------------------------------------------------------

    parameters              description                    Set/Default
    ==========              ===========                    ========================
    --bamdir                bam directory                  ${params.bamdir}
    --sample_sheet          sample sheet                   ${params.sample_sheet}
    CPUs                                = ${params.cores}
    debug                               = ${params.debug}
    Region file.                        = ${params.interval_bed}
    Output folder                       = ${params.out}

"""

out
}

if (params.help) {
    log_summary()
    exit 1
}

if (workflow.profile == "") {
    println "Must set -profile: local, quest, gcp"
    exit 1
}

strains = params.strains ? params.strains.split(",") : false

// // Define contigs here!
// CONTIG_LIST = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
// contigs = Channel.from(CONTIG_LIST)

// /*
// ~ ~ ~ > * gff file for snpeff
// */

// cegff = Channel.fromPath(params.gff)

/*
===========================================
~ ~ ~ > * Generate Interval List  * < ~ ~ ~ 
===========================================
*/

process split_genome {

    label 'md'

    conda "gatk4=4.1.4.0"

    output:
        path "scatter/*-scattered.interval_list"

    script:
        intervals = params.interval_bed != "" ?  "-L ${params.interval_bed}" : ""

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
        SplitIntervals \\
        -R ${reference_uncompressed} \\
        ${intervals} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
        --scatter-count 20 \\
        -ip 250 \\
        -O scatter
    """
}

process get_contigs {

    label 'sm'

    input:
        tuple strain, ref_strain, path(bam), path(bai)

    output:
        path("contigs.txt")

    """
        samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    """

}

// Read sample sheet
sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                      .map { row -> [row.strain,
                                     row.reference_strain == "TRUE",
                                     file("${params.bamdir}/${row.strain}.bam", checkExists: true),
                                     file("${params.bamdir}/${row.strain}.bam.bai", checkExists: true)] }
                      .distinct()
                      .filter { params.strains ? it[0] in strains : true }

workflow {

    // Get contigs from first bam
    sample_sheet.first() | get_contigs
    contigs = get_contigs.out.splitText { it.strip() }

    // Call individual variants
    sample_sheet.combine(contigs) | call_variants_individual

    call_variants_individual.out.groupTuple()
                                .map { strain, ref_strain, vcf -> [strain, ref_strain[0], vcf]}
                                .combine(get_contigs.out) | \
                                concat_strain_gvcfs
    
    // gatk genomics db
    sample_map = sample_sheet.map { "${it[0]}\t${it[0]}.g.vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)
    concat_strain_gvcfs.out.flatten()
                           .toList()
                           .map { [it] }
                           .combine(contigs)
                           .combine(sample_map)
                           .view() | import_genomics_db


}

/*
=============================================
~ ~ ~ > * Run GATK4 HaplotypeCaller * < ~ ~ ~ 
=============================================
*/

process call_variants_individual {

    label 'lg'

    tag { "${strain}:${region}" }

    input:
        tuple strain, ref_strain, path(bam), path(bai), val(region)

    output:
        tuple strain, ref_strain, path("${region}.g.vcf.gz")

    """
        gatk HaplotypeCaller --java-options "-Xmx${task.cpus}g -Xms1g" \\
            --emit-ref-confidence GVCF \\
            --max-genotype-count 3000 \\
            --max-alternate-alleles 100 \\
            --annotation DepthPerAlleleBySample \\
            --annotation Coverage \\
            --annotation GenotypeSummaries \\
            --annotation TandemRepeat \\
            --annotation StrandBiasBySample \\
            --annotation ChromosomeCounts \\
            --annotation AS_QualByDepth \\
            --annotation AS_StrandOddsRatio \\
            --annotation AS_MappingQualityRankSumTest \\
            --annotation DepthPerSampleHC \\
            --annotation-group StandardAnnotation \\
            --annotation-group AS_StandardAnnotation \\
            --annotation-group StandardHCAnnotation \\
            -R ${reference_uncompressed} \\
            -I ${bam} \\
            -L ${region} \\
            -O ${region}.g.vcf   
        bcftools view -O z ${region}.g.vcf > ${region}.g.vcf.gz
    """
}

// individual_sites
//   .groupTuple()
//   .join(merged_sample_maps)
//   .into{merge_sm_int_gvcfs;
//        print_sm_int_gvcfs}

/*
=============================================
~ ~ ~ > *  Merge Sample gVCF files  * < ~ ~ ~ 
=============================================
*/

process concat_strain_gvcfs {

    label 'sm'
    tag { "${strain}" }
    //conda 'bcftools=1.9'

    input:
        tuple strain, ref_strain, path("*"), path(contigs)

    output:
        tuple path("${strain}.g.vcf.gz"), path("${strain}.g.vcf.gz.csi")

    """
        awk '{ print \$0 ".g.vcf.gz" }' ${contigs} > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > ${strain}.g.vcf.gz
        bcftools index ${strain}.g.vcf.gz
    """
}

// merged_sample_gvcf
//   .toSortedList()
//   .into{merged_sample_gvcf_to_sample_map;
//         merged_sample_gvcf_to_db;
//         merged_sample_gvcf_other;}


// merged_sample_gvcf_index
//   .toSortedList()
//   .into{merged_sample_gvcf_index_to_sample_map;
//         merged_sample_gvcf_index_to_db;
//         merged_sample_gvcf_index_to_other;}


/*
===============================================
~ ~ ~ > *  Combine All Samples to DB  * < ~ ~ ~ 
===============================================
*/

//  process cohort_to_sample_map {

//     memory '64 GB'

//     publishDir "${params.out}/Isotype_gVCF", mode: 'copy'

//     executor 'local'

//     input:
//       file(gvcfs) from merged_sample_gvcf_to_sample_map
//       file(indices) from merged_sample_gvcf_index_to_sample_map

//     output:
//       file("cohort.sample_map") into cohort_map


//     """
//       ls *_merged-intervals.g.vcf > cohort.sample_map_tmp
//       awk -F"_" 'BEGIN{OFS="\t"}; \$1=\$1' OFS="\t" cohort.sample_map_tmp | \\
//       awk '{print \$1, \$1"_"\$2}' OFS="\t" > cohort.sample_map
//     """
// }

 process import_genomics_db {

    tag { "${contig}" }
    label 'lg'

    input:
        tuple path(vcfs), contig, path(sample_map)

    output:
      set val(chr), file("${contig}_database.tar")

    """
        gatk  --java-options "-Xmx${task.memory.toGiga()-1}g -Xms${task.memory.toGiga()-2}g" \\
            GenomicsDBImport \\
            --genomicsdb-workspace-path ${contig}_database \\
            --batch-size 16 \\
            -L ${contig} \\
            --sample-name-map ${sample_map} \\
            --reader-threads ${task.cpus}

      tar -cf ${contig}_database.tar ${contig}_database
    """
}

/*
====================================
~ ~ ~ > *  Genotype gVCFs  * < ~ ~ ~ 
====================================
*/

//  process genotype_cohort_gvcf_db {

//     memory '64 GB'

//     tag { "${chr}" }

//     cpus 12

//     input:
//       set val(chr), file(chr_db) from chromosomal_db

//     output:
//       set val(chr), file("${chr}_cohort.vcf"), file("${chr}_cohort.vcf.idx") into cohort_vcf

//     """
//       tar -xf ${chr_db}
//       WORKSPACE=\$( basename ${chr_db} .tar)

//       gatk --java-options "-Xmx48g -Xms48g" \\
//         GenotypeGVCFs \\
//         -R ${reference_uncompressed} \\
//         -V gendb://\$WORKSPACE \\
//         -G StandardAnnotation \\
//         -G AS_StandardAnnotation \\
//         -G StandardHCAnnotation \\
//         -L ${chr} \\
//         --use-new-qual-calculator \\
//         -O ${chr}_cohort.vcf
//     """
// }

// /*
// =======================================================
// ~ > *                                             * < ~
// ~ ~ > *                                         * < ~ ~
// ~ ~ ~ > *  SnpEff Annotate Cohort Chrom VCFs  * < ~ ~ ~ 
// ~ ~ > *                                         * < ~ ~ 
// ~ > *                                             * < ~
// =======================================================
// */


// process annotate_vcf {

//     memory '64 GB'

//     cpus 12

//     tag { chrom }

//     input:
//         set val(chrom), file(vcf), file(vcfindex) from cohort_vcf
//         file 'snpeff.config' from file("${baseDir}/snpeff_data/snpEff.config")
//         file 'snpeff_data' from file("${baseDir}/snpeff_data")

//     output:
//         set file("${chrom}.annotated.vcf.gz"), file("${chrom}.annotated.vcf.gz.tbi") into annotated_vcf
//         file("snpeff_out.csv") into snpeff_multiqc


//     """
//       bcftools view -Oz -o ${vcf}.gz ${vcf}
//       tabix -p vcf ${vcf}.gz

//       # bcftools csq
//       bcftools view --threads=${task.cpus-1} --regions ${chrom} ${vcf}.gz | \\
//       snpEff eff -csvStats snpeff_out.csv \\
//                  -no-downstream \\
//                  -no-intergenic \\
//                  -no-upstream \\
//                  -nodownload \\
//       -dataDir . \\
//       -config snpeff.config \\
//       ${params.annotation_reference} | \\
//       bcftools view --threads=${task.cpus-1} -O z > ${chrom}.annotated.vcf.gz
      
//       tabix -p vcf ${chrom}.annotated.vcf.gz
//     """

// }

// /*
// =================================================
// ~ > *                                       * < ~
// ~ ~ > *                                   * < ~ ~
// ~ ~ ~ > *   Concatenate Annotated VCFs  * < ~ ~ ~
// ~ ~ > *                                   * < ~ ~
// ~ > *                                       * < ~
// =================================================
// */

// contig_raw_vcf = CONTIG_LIST*.concat(".annotated.vcf.gz")

// process concatenate_annotated_vcf {

//     memory '64 GB'

//     cpus 12

//     input:
//         file(merge_vcf) from annotated_vcf.collect()

//     output:
//         set file("WI.annotated.vcf.gz"), file("WI.annotated.vcf.gz.tbi") into annotated_concatenated_vcf
//         set val("soft"), file("WI.annotated.vcf.gz"), file("WI.annotated.vcf.gz.tbi") into annotated_sample_summary
//         file("WI.annotated.stats.txt") into annotated_stats

//     """
//       bcftools concat --threads ${task.cpus-1} -O z ${contig_raw_vcf.join(" ")} > WI.annotated.vcf.gz
//       tabix -p vcf WI.annotated.vcf.gz
//       bcftools stats --verbose WI.annotated.vcf.gz > WI.annotated.stats.txt
//     """
// }

// annotated_concatenated_vcf
//   .into{ann_vcf_to_soft_filter;
//         ann_vcf_to_hard_filter;
//         ann_vcf_to_vqsr;
//         ann_vcf_to_recal;
//         ann_vcf_to_cnn_train;
//         ann_vcf_to_cnn_apply;}


// /*
// ================================================
// ~ > *                                      * < ~
// ~ ~ > *                                  * < ~ ~
// ~ ~ ~ > *   Apply Soft Filters on VCF  * < ~ ~ ~
// ~ ~ > *                                  * < ~ ~
// ~ > *                                      * < ~
// ================================================
// */

// /*
// ============================================
// ~ ~ ~ > *   Split Indels and SNVs  * < ~ ~ ~
// ============================================
// */

// process split_snv_indel {

//     memory '64 GB'

//     cpus 12

//     input:
//       set file(merge_ann_vcf), file(merge_ann_index) from ann_vcf_to_soft_filter

//     output:
//       set file("merged_annotated_SNV.vcf"), file("merged_annotated_SNV.vcf.idx") into split_snv_vcf
//       set file("merged_annotated_INDEL.vcf"), file("merged_annotated_INDEL.vcf.idx") into split_indel_vcf


//   """
//     gatk --java-options "-Xmx48g -Xms48g" \\
//         SelectVariants \\
//         -R ${reference_uncompressed} \\
//         -V ${merge_ann_vcf} \\
//         --select-type-to-include SNP \\
//         -O merged_annotated_SNV.vcf


//     gatk --java-options "-Xmx48g -Xms48g" \\
//         SelectVariants \\
//         -R ${reference_uncompressed} \\
//         -V ${merge_ann_vcf} \\
//         --select-type-to-include INDEL \\
//         -O merged_annotated_INDEL.vcf
//   """

// }

// /*
// ================================================
// ~ ~ ~ > *   Apply Indels Soft Filters  * < ~ ~ ~
// ================================================
// */

// process apply_soft_filters_indel {

//     cpus params.cpu

//     input:
//       set file(indel_vcf), file(indel_index) from split_indel_vcf

//     output:
//       set file("soft_filtered_indels.vcf"), file("soft_filtered_indels.vcf.idx") into soft_filter_indels


//     """
//       bcftools norm --threads ${task.cpus} -m -any -Ov -o norm_indels.vcf ${indel_vcf} 

//       gatk --java-options "-Xmx48g -Xms48g" \\
//           VariantFiltration \\
//           -R ${reference_uncompressed} \\
//           --variant norm_indels.vcf \\
//           --genotype-filter-expression "DP < ${params.min_depth}" \\
//           --genotype-filter-name "depth" \\
//           --filter-expression "QD < ${params.quality_by_depth}" \\
//           --filter-name "qd" \\
//           --filter-expression "SOR > ${params.strand_odds_ratio}" \\
//           --filter-name "sor" \\
//           --filter-expression "QUAL < ${params.qual}" \\
//           --filter-name "quality" \\
//           --filter-expression "ReadPosRankSum < ${params.readbias}" \\
//           --filter-name "readend" \\
//           --filter-expression "FS > ${params.fisherstrand}" \\
//           --filter-name "fisherstrand" \\
//           -O soft_filtered_indels.vcf
//     """
// }

// /*
// =============================================
// ~ ~ ~ > *   Apply SNV Soft Filters  * < ~ ~ ~
// =============================================
// */

// process apply_soft_filters_snv {

//     cpus params.cpu

//     input:
//       set file(snp_vcf), file(snp_index) from split_snv_vcf

//     output:
//       set file("soft_filtered_snps.vcf"), file("soft_filtered_snps.vcf.idx") into soft_filter_snvs


//     """
//       gatk --java-options "-Xmx48g -Xms48g" \\
//           VariantFiltration \\
//           -R ${reference_uncompressed} \\
//           --variant ${snp_vcf} \\
//           --genotype-filter-expression "DP < ${params.min_depth}" \\
//           --genotype-filter-name "depth" \\
//           --filter-expression "QUAL < ${params.qual}" \\
//           --filter-name "quality" \\
//           --filter-expression "ReadPosRankSum < ${params.readbias}" \\
//           --filter-name "readend" \\
//           --filter-expression "FS > ${params.fisherstrand}" \\
//           --filter-name "fisherstrand" \\
//           --filter-expression "QD < ${params.quality_by_depth}" \\
//           --filter-name "qd" \\
//           --filter-expression "SOR > ${params.strand_odds_ratio}" \\
//           --filter-name "sor" \\
//           -O soft_filtered_snps.vcf
//     """
// }

// /*
// ==============================================
// ~ ~ ~ > *   Combine SNVs and Indels  * < ~ ~ ~
// ==============================================
// */

// process combine_soft_filter_vcfs {

//     cpus params.cpu

//     input:
//       set file("soft_filtered_snps.vcf"), file("soft_filtered_snps.vcf.isx") from soft_filter_snvs
//       set file("soft_filtered_indels.vcf"), file("soft_filtered_indels.vcf.idx") from soft_filter_indels

//     output:
//       set file("WI.${date}.soft-filter.vcf.gz"), file("WI.${date}.soft-filter.vcf.gz.tbi") into soft_filtered_cohort_vcf
//       file("WI.${date}.soft-filter.stats.txt") into soft_filtered_stats_to_mqc


//     """
//       bcftools view -Oz -o soft_filtered_indels.vcf.gz soft_filtered_indels.vcf
//       tabix -p vcf soft_filtered_indels.vcf.gz

//       bcftools view -Oz -o soft_filtered_snps.vcf.gz soft_filtered_snps.vcf
//       tabix -p vcf soft_filtered_snps.vcf.gz

//       bcftools concat \\
//       --threads ${task.cpus-1} \\
//       --allow-overlaps \\
//       soft_filtered_indels.vcf.gz \\
//       soft_filtered_snps.vcf.gz | \\
//       bcftools filter -Oz --threads ${task.cpus-1} --mode + --soft-filter high_missing --include "F_MISSING<=${params.missing}" > WI.${date}.soft-filter.vcf.gz

//       tabix -p vcf WI.${date}.soft-filter.vcf.gz
//       bcftools stats --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt
//     """
// }

// soft_filtered_cohort_vcf
//   .into{soft_vcf_to_strain_list;
//         soft_vcf_to_split_by_strain;
//         soft_vcf_to_other}

// /*
// ==========================================
// ~ ~ ~ > *   Split VCF by sample  * < ~ ~ ~
// ==========================================
// */

// process generate_strain_list {

//     executor 'local'

//     input:
//         set file(softvcf), file(softvcf_index) from soft_vcf_to_strain_list

//     output:
//         file('isotype_list.tsv') into isotype_list

//     """
//         bcftools query -l ${softvcf} > isotype_list.tsv
//     """

// }

// isotype_list
//   .splitText() { it.strip() } 
//   .combine(soft_vcf_to_split_by_strain)
//   .into{isotype_set_vcf; 
//         isotype_set_tsv}

// /*
// ===========================================
// ~ ~ ~ > *   Apply AD soft filter  * < ~ ~ ~
// ===========================================
// */

// process apply_allele_depth_filter {

//     tag { isotype }

//     memory '64 GB'

//     input:
//         set val(isotype), file(softvcf), file(softvcf_vcf) from isotype_set_vcf

//     output:
//         file("${isotype}.AD-filter.vcf.gz") into isotype_AD_soft_vcf
//         file("${isotype}.AD-filter.vcf.gz.tbi") into isotype_AD_soft_vcf_index

//     """
//     bcftools view --samples ${isotype} ${softvcf} |\\
//     bcftools filter -Ov --mode + --soft-filter dv_dp --include "((FORMAT/AD[*:1])/(FORMAT/DP) >= 0.5) || (FORMAT/GT == '0/0') || (TYPE == 'REF')" -Ov -o ${isotype}_temp.vcf

//     gatk --java-options "-Xmx4g -Xms4g" \\
//          VariantFiltration \\
//          -R ${reference_uncompressed} \\
//          --variant ${isotype}_temp.vcf \\
//          --genotype-filter-expression "FILTER != 'PASS'" \\
//          --genotype-filter-name "dv_dp" \\
//          -O temp_soft_filtered.vcf

//     awk 'BEGIN{FS=OFS="\\t"} {gsub("dv_dp",\$7,\$10)} 1' temp_soft_filtered.vcf | \\
//     bcftools view -Oz -o ${isotype}.AD-filter.vcf.gz

//     tabix -p vcf ${isotype}.AD-filter.vcf.gz
//     """

// }

// contigs_soft = Channel.from(CONTIG_LIST)

// contigs_soft
//   .into{contigs_vcf;
//         contigs_index;
//         contigs_impute;
//         }

// contigs_vcf
//   .spread(isotype_AD_soft_vcf)
//   .groupTuple()
//   .set{ to_merge_soft_sm_ad }

// contigs_index
//   .spread(isotype_AD_soft_vcf_index)
//   .groupTuple()
//   .into{ to_merge_soft_sm_ad_index;
//         print_merged_index }

// /*
// ==========================================
// ~ ~ ~ > *   Combine Cohort VCFs  * < ~ ~ ~
// ==========================================
// */

// process merge_sm_soft_vcfs {

//     tag { chrom }

//     memory '64 GB'
//     cpus 4

//     input:
//         set val(chrom), file(softvcf) from to_merge_soft_sm_ad
//         set val(chrom), file(softvcf_vcf) from to_merge_soft_sm_ad_index

//     output:
//         file("WI.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz") into cohort_soft_filter_vcf
//         file("WI.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into cohort_soft_filter_vcf_index

//     """
//         bcftools merge \\
//         -m none \\
//         -r ${chrom} \\
//         -Oz -o WI.${chrom}.COMPLETE.vcf.gz \\
//         ${softvcf}

//         bcftools query -l WI.${chrom}.COMPLETE.vcf.gz | sort > samples.txt

//         bcftools view -S samples.txt WI.${chrom}.COMPLETE.vcf.gz -Oz -o WI.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz

//         tabix -p vcf WI.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz
//     """
// }

// cohort_soft_filter_vcf
//   .toSortedList()
//   .set{ to_concat_soft_filter_vcf }

// cohort_soft_filter_vcf_index
//   .toSortedList()
//   .set{ to_concat_soft_filter_vcf_index }

// /*
// ====================================================
// ~ ~ ~ > *   Concatenate Cohort CHROM VCFs  * < ~ ~ ~
// ====================================================
// */

// process concat_cohort_soft_vcfs {

//     publishDir "${params.out}/variation/", mode: 'copy'

//     memory '64 GB'
//     cpus 20

//     input:
//         file(chromvcf) from to_concat_soft_filter_vcf
//         file(chromvcf_index) from to_concat_soft_filter_vcf_index

//     output:
//         set file("WI.COMPLETE-SOFT-FILTER.vcf.gz"), file("WI.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into cohort_complete_soft_filter_vcf
//         set val("soft"), file("WI.COMPLETE-SOFT-FILTER.vcf.gz"), file("WI.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into soft_vcf_summary
//         set val("soft"), file("WI.COMPLETE-SOFT-FILTER.vcf.gz"), file("WI.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into soft_sample_summary
//         file("WI.COMPLETE-SOFT-FILTER.stats.txt") into complete_soft_filter_vcf_stats

//     """
//       bcftools concat \\
//       --threads ${task.cpus} \\
//       WI.I.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.II.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.III.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.IV.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.V.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.X.COMPLETE-SOFT-FILTER.vcf.gz \\
//       WI.MtDNA.COMPLETE-SOFT-FILTER.vcf.gz |\\
//       bcftools filter -Oz --mode=+x --soft-filter="high_missing" --include 'F_MISSING  <= ${params.missing}' | \\
//       bcftools filter -Oz --mode=+x --soft-filter="high_heterozygosity" --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' |\\
//       bcftools view -Oz -o WI.COMPLETE-SOFT-FILTER.vcf.gz

//       tabix -p vcf WI.COMPLETE-SOFT-FILTER.vcf.gz
//       bcftools stats --verbose WI.COMPLETE-SOFT-FILTER.vcf.gz > WI.COMPLETE-SOFT-FILTER.stats.txt
//     """
// }

// cohort_complete_soft_filter_vcf
//   .into{soft_filtered_vcf_to_hard;
//         soft_filtered_vcf_impute;
//         soft_filter_vcf_strain;
//         soft_filter_vcf_isotype_list;
//         soft_filter_vcf_mod_tracks;
//         soft_filter_vcf_tsv
//       }

// /*
// ================================================
// ~ > *                                      * < ~
// ~ ~ > *                                  * < ~ ~
// ~ ~ ~ > *   Apply Hard Filters on VCF  * < ~ ~ ~
// ~ ~ > *                                  * < ~ ~
// ~ > *                                      * < ~
// ================================================
// */

// process generate_hard_vcf {

//     cpus 20

//     publishDir "${params.out}/variation", mode: 'copy'

//     input:
//         set file(softvcf), file(softvcfindex) from soft_filtered_vcf_to_hard

//     output:
//         set file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.tbi") into hard_vcf
//         set val("hard"), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.tbi") into hard_vcf_summary
//         set val("hard"), file("WI.${date}.hard-filter.vcf.gz"), file("WI.${date}.hard-filter.vcf.gz.tbi") into hard_sample_summary
//         file("WI.${date}.hard-filter.stats.txt") into hard_filter_stats


//     """
//         # Generate hard-filtered VCF
//         function generate_hard_filter {
//             bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${softvcf} |\\
//             bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT != "PASS"' |\\
//             bcftools filter -O u --include 'F_MISSING  <= ${params.missing}' |\\
//             bcftools filter -O u --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' |\\
//             bcftools view -O v --min-af 0.0000000000001 --max-af 0.999999999999 |\\
//             vcffixup - | \\
//             bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
//         }

//         export -f generate_hard_filter

//         parallel --verbose generate_hard_filter {} ::: I II III IV V X MtDNA

//         bcftools concat -O z I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz > WI.${date}.hard-filter.vcf.gz
        
//         tabix -p vcf WI.${date}.hard-filter.vcf.gz

//         bcftools stats --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt

//         # Remove extra files
//         rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
//     """
// }

// /*
// ================================================
// ~ > *                                      * < ~
// ~ ~ > *                                  * < ~ ~
// ~ ~ ~ > *   Apply Hard Filters on VCF  * < ~ ~ ~
// ~ ~ > *                                  * < ~ ~
// ~ > *                                      * < ~
// ================================================
// */


// process imputation {

//     tag {CHROM}

//     cpus params.cores

//     input:
//         set file(softvcf), file(softvcfindex) from soft_filtered_vcf_impute
//         each CHROM from contigs_impute

//     output:
//         file("${CHROM}_imp.vcf.gz") into impute_vcf
//         file("${CHROM}_imp.vcf.gz.csi") into impute_index

//     """

//     java -Xmx50g -jar `which beagle.jar` chrom=${CHROM} window=8000 overlap=3000 impute=true ne=17500 gt=${softvcf} out=${CHROM}_imp

//     bcftools index ${CHROM}_imp.vcf.gz

//     """
// }



// impute_vcf
//   .toSortedList()
//   .set{ imp_v_to_concat }


// impute_index
//   .toSortedList()
//   .set{ imp_i_to_concat }


// process concat_imputed {

//     publishDir "${params.out}/variation/", mode: 'copy'

//     memory '64 GB'
//     cpus 20

//     input:
//         file(vcf) from imp_v_to_concat
//         file(index) from imp_i_to_concat

//     output:
//         set file("WI.${date}.impute.vcf.gz"), file("WI.${date}.impute.vcf.gz.tbi") into cohort_impute

//     """
//       bcftools concat \\
//       --threads ${task.cpus} \\
//       I_imp.vcf.gz \\
//       II_imp.vcf.gz \\
//       III_imp.vcf.gz \\
//       IV_imp.vcf.gz \\
//       V_imp.vcf.gz \\
//       X_imp.vcf.gz \\
//       MtDNA_imp.vcf.gz |\\
//       bcftools view -Oz -o WI.${date}.impute.vcf.gz

//       tabix -p vcf WI.${date}.impute.vcf.gz
//     """
// }

process CREATE_BAM_SAMPLE_SHEET {
    label 'create_bam_sample_sheet'


    input:
    path input
    path bam_folder
    path gvcf_folder

    output:
    path("sample_sheet_bam.txt"),  emit: bam
    path("sample_sheet_gvcf.txt"), emit: gvcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    STRAINS=(`cat ${input}`)
    MISSING="FALSE"
    rm -f sample_sheet_bam.txt
    rm -f sample_sheet_gvcf.txt
    touch sample_sheet_bam.txt
    touch sample_sheet_gvcf.txt
    for S in \${STRAINS[*]}; do
        if [[ ! -e ${gvcf_folder}/\${S}.g.vcf.gz ]]; then
            if [[ ! -e ${bam_folder}/\${S}.bam ]]; then
                echo "Missing \${S}.bam" 1>&2
                MISSING="TRUE"
            fi
            if [[ ! -e ${bam_folder}/\${S}.bam.bai ]]; then
                echo "Missing \${S}.bam" 1>&2
                MISSING="TRUE"
            fi
            if [[ ! -e sample_sheet_bam.txt ]]; then
                echo "\${S}" > sample_sheet_bam.txt
            else
                echo "\${S}" >> sample_sheet_bam.txt
            fi
        else
            echo "\${S}" >> sample_sheet_gvcf.txt
        fi
    done
    if [[ \${MISSING} == "TRUE" ]]; then
        exit 1
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch sample_sheet_bam.txt
    touch sample_sheet_gvcf.txt
    """
}
process LOCAL_SUMMARY {
    executor 'local'
    container null

    input:
        val run
        path "sample_sheet"

    output:
        path "sample_sheet_${date}.tsv", emit: sample_sheet
        path "summary_${date}.txt",      emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    echo '''${log_summary()}''' > summary.txt
    cat ${sample_sheet} > sample_sheet.tsv
    """

    stub:
    """
    touch summary.txt
    touch sample_sheet.tsv
    """
}

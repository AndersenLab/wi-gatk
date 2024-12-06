process LOCAL_REPORT {
    label 'local_report'
    errorStrategy 'retry'
    time { 1.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    path("multiqc.html")
    path(soft_filter_filter)
    path(soft_filter_stats)
    path(hard_filter_stats)
    path("gatk_report.Rmd")
    val timestamp

    output:
    path "gatk_report_${timestamp}.html", emit: html
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat "gatk_report.Rmd" | \\
        sed -e 's/RELEASE_DATE/${timestamp}/g' |
        sed -e 's+variation/WI.{vcf_date}.soft-filter.stats.txt+${soft_filter_stats}+g' |
        sed -e 's+variation/WI.{vcf_date}.hard-filter.stats.txt+${hard_filter_stats}+g' |
        sed -e 's+variation/WI.{vcf_date}.soft-filter.filter_stats.txt+${soft_filter_filter}+g' > gatk_report_${timestamp}.Rmd
    Rscript -e "rmarkdown::render('gatk_report_${timestamp}.Rmd')"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$( Rscript --version | cut -f 4 )
    END_VERSIONS
    """

    stub:
    """
    touch gatk_report_${timestamp}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$( Rscript --version | cut -f 4 )
    END_VERSIONS
    """
}
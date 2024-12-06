process MULTIQC_REPORT {
    label 'multiqc_report'
    errorStrategy 'retry'
    time { 1.hour * task.attempt }
    cpus = { 1 * task.attempt }
    memory = { 4.GB * task.attempt }

    input:
    path("*")

    output:
    path "multiqc_data/*.json", emit: json, optional: true
    path "multiqc.html",        emit: report
    path "versions.yml",        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -k json --filename multiqc.html .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
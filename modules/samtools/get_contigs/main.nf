process SAMTOOLS_GET_CONTIGS {
    label 'samtools_get_contigs'

    input:
    tuple val(meta), path(bam), path(bam_index)
    val partition_size

    output:
    path("contigs.txt"),     emit: contigs
    path("contig_partitions.tsv"), emit: partitions
    path  "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    samtools view -H ${bam} | \\
        awk '{ \\
            if (( \$1 ~ /^@SQ/ )){ \\
                for (I=2;I<=NF;I++){ \\
                    if (( \$I ~ /^SN/ )){ \\
                        split(\$I,A,":"); NAME=A[2]; \\
                    } else if (( \$I ~ /^LN/ )){ \\
                        split(\$I,A,":"); LEN=A[2]; \\
                    } \\
                } \\
                printf "%s\\t%i\\n", NAME, LEN; \\
            }}' > chrom_sizes.txt
    CHROMS=(\$(cut -f 1 chrom_sizes.txt))
    SIZES=(\$(cut -f 2 chrom_sizes.txt))
    for I in \$(seq 0 1 \$(expr \${#CHROMS[*]} - 1 )); do
        if [ \${SIZES[\${I}]} -lt ${partition_size} ]; then
            N=1
        else
            N=\$(expr \${SIZES[\${I}]} / ${partition_size} )
        fi
        STEP=\$(expr \${SIZES[\${I}]} / \$N )
        START=1
        for J in \$(seq 1 1 \$(expr \${N} - 1 )); do
            END=\$(expr \${START} + \${STEP} )
            END2=\$(expr \${END} + 300 )
            echo -e "\${CHROMS[\${I}]}\\t\${SIZES[\${I}]}\\t\${START}\\t\${END2}" >> contig_partitions.tsv
            START=\${END}
        done
        END=\${SIZES[\${I}]}
        echo -e "\${CHROMS[\${I}]}\\t\${SIZES[\${I}]}\\t\${START}\\t\${END}" >> contig_partitions.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch contigs.txt
    touch contig_partitions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
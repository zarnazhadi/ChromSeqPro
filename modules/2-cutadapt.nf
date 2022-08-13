/*
#==============================================
Quality, read and adapter trimming of FASTQ
files
#==============================================
*/

process CUTADAPT {

    publishDir "${params.TRIM_DIR}", mode: 'copy'
    cpus 24

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    tuple val("trimmed_${sample_id}"), file("trimmed_${r1}"), file("trimmed_${r2}")

    script:
    """
    mkdir -p ${params.TRIM_DIR}
    cutadapt -a AGATCGGAAGAG \
             -A AGATCGGAAGAG \
             -o trimmed_${r1} \
             -p trimmed_${r2} \
             ${r1} \
             ${r2} \
             --minimum-length=${params.minimumLength} \
             --overlap=${params.overlap} \
             -q ${params.qscore} \
             -j ${task.cpus}
    """
}

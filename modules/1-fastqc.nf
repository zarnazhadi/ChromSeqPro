/*
#==============================================
Generate FastQC	reports	for each FASTQ file
#==============================================
*/

process FASTQC {

    publishDir "${params.OUT_DIR}", mode: 'copy'

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    path "01_FASTQC/${sample_id}_report"

    script:
    """
    mkdir -p 01_FASTQC
    mkdir 01_FASTQC/${sample_id}_report
    /FastQC/fastqc -o 01_FASTQC/${sample_id}_report -f fastq -q ${r1} ${r2}
    """

}

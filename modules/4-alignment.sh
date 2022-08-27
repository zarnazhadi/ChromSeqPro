process BWA_MEM {

    publishDir "${params.OUT_DIR}", mode: 'copy'
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.sam"

    script:
    """
    bwa mem -M -t ${task.cpus} ${params.reference} \
	${reads[0]} ${reads[1]} > aligned_${sample_id}.sam	

    """
}

/*workflow {

	trimmed_list    = "${params.TRIM_DIR}/*R{1,2}*"

        Channel
        .fromFilePairs(trimmed_list, checkIfExists: true)
        .set { trim_pairs_ch }

        BWA_MEM(trim_pairs_ch)
}
*/

process PICARD {

    publishDir("${params.OUT_DIR}/06_PICARD", mode: 'copy')

    input:
    path(reads)

    output:
    path "*-marked.bam"
    path "*-metrics.txt"

    script:
    if(!params.peakCall)
    """
    mkdir -p 06_PICARD
    
    java -jar /picard.jar MarkDuplicates --INPUT ${reads} --REMOVE_DUPLICATES --REMOVE_SEQUENCING_DUPLICATES \
	--OUTPUT ${reads.baseName}-marked.bam \
	--METRICS_FILE ${reads.baseName}-metrics.txt 
	--TMP_DIR ${params.OUT_DIR}
    """
    else
    """
    java -jar /picard.jar MarkDuplicates --INPUT \
        --OUTPUT ${reads.baseName}-marked.bam \
        --METRICS_FILE ${reads.baseName}-metrics.txt
        --TMP_DIR ${params.OUT_DIR}
    """
}

workflow {

	bam_list = '../bash/aligned/*.sorted.bam'

        Channel
        .fromPath(bam_list, checkIfExists: true)
        .set {bam_ch}

	PICARD(bam_ch)
}

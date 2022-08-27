/*
#==============================================
Generate Samtools summaries for BAM files 
#==============================================
*/

process SORT {
 
    publishDir "${params.OUT_DIR}/05_SAMTOOLS/${bamName}", mode: 'copy'

    input:
    file(bamRead)

    output:
    path "*.bam"

    script:

    bamName = bamRead.toString().split("\\.")[0]

    """
    samtools sort -m 5G -o ${bamName}.sort.bam ${bamRead}
    """
}

process SAM_INDEX {

    publishDir "${params.OUT_DIR}/05_SAMTOOLS/${bamName}", mode: 'copy'

    input:
    file(bamRead)

    output:
    file("*.bai")

    script:

    bamName = bamRead.toString().split("\\.")[0]

    """
    samtools index ${bamRead}
    """
}

process VIEW {

    
    publishDir "${params.OUT_DIR}/05_SAMTOOLS/${bamName}", mode: 'copy'

    input:
    file(bamRead)

    output:
    file("*.bam")

    script:
    
    bamName = bamRead.toString().split("-")[0]

    """
    samtools view -F 1804 -q ${params.mapq_thresh} -b ${bamRead} -o ${bamName}.view.bam
    """

}

process IDXSTAT {


    publishDir "${params.OUT_DIR}/05_SAMTOOLS/${samName}", mode: 'copy'

    input:
    file(samRead)

    output:
    file("*.txt")

    script:

    samName = samRead.toString().split("\\.")[0]

    """
    samtools idxstat ${samRead} > ${samName}.stat.txt
    """

}

process FLAGSTAT {


    publishDir "${params.OUT_DIR}/05_SAMTOOLS/${samName}", mode: 'copy'

    input:
    file(samRead)

    output:
    file("*.txt")

    script:

    samName = samRead.toString().split("\\.")[0]

    """
    samtools flagstat ${samRead} > ${samName}-flagstat.txt
    """

}


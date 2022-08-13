/*
#==============================================
Rename files to SAMPLE_CTYPE_R.fastq.gz - 
requires definition of params.sampleList in the 
command line.
#==============================================
*/

process RENAME {

        publishDir "${params.REG_DIR}", mode: 'copy'

	label "small_job"

        input:
              	tuple val(sampleID),
                        val(experiment),
                        val(cType),
                        file(r1),
                        file(r2)

        output:
               	tuple val("${sampleID}"), file("${sampleID}_${cType}_R1.fastq.gz"), file("${sampleID}_${cType}_R2.fastq.gz")

        script:
        """
	      cp ${r1} ${sampleID}_${cType}_R1.fastq.gz
        cp ${r2} ${sampleID}_${cType}_R2.fastq.gz
	"""
}

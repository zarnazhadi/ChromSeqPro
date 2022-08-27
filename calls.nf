include { MACS2 ; MACS2 as MACS3 } from './modules/7-macs2.nf'
include { REGION_COUNTS; PROMOTER_SIGNAL } from '.modules/7-chrompro.nf'

workflow PEAK_CALL {


	recurrent_list = "../bash/aligned/mem/{EZH2,H3K4me3,H3K27me3}_R-viewed.sort.bam"
        primary_list = "../bash/aligned/mem/EZH2_P-viewed.sort.bam"

        Channel
        .fromPath(recurrent_list, checkIfExists:true)
        .map {it -> [it, file("../bash/aligned/mem/IMPUT_R-viewed.sort.bam", checkIfExists: true)]}
        .set { recurrent_ch }

        Channel
        .fromPath(primary_list, checkIfExists:true)
        .map {it -> [it, file("../bash/aligned/mem/IMPUT_P-viewed.sort.bam", checkIfExists: true)]}
        .set { primary_ch }

        MACS2(recurrent_ch)
        MACS3(primary_ch)

}

workflow CHROM_CALL {

        	println """\
         C H I P S E Q - N F   P I P E L I N E
         ===================================
         base		        : ${params.BASE}
         metadata     	: ${params.METADATA_DIR}
	       results 	      : ${params.RES_DIR}
         """
         .stripIndent()	

	      bam_ch = Channel.fromPath("/nobackup/bs21zh/ZH_prelim/bash/filtered/*.bam", checkIfExists: true)	
	      REGION_COUNTS(bam_ch)

      /*	Channel.fromPath(REGION_COUNTS.out)
    	| branch { path ->
        	imputed: path =~ /control_[RP]-filtered.bam/
        	calculated: true
    	} | set { bams }

    	bams.imputed
    	| cross(bams.calculated) { it.name.find(/[RP]-filtered.bam/) }
    	| view
	
      */
	      PROMOTER_SIGNAL(REGION_COUNTS.out)
}

workflow {

	if (!params.peakCall) {
                CHROM_CALL()
        } else {
                PEAK_CALL()
        }

}


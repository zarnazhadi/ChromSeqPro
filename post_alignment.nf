/*
#==============================================
Import processes
#==============================================
*/

include { SORT; VIEW; FLAGSTAT } from './modules/5-samtools.nf'
include { PICARD               } from './modules/6-picard.nf'

workflow {

	bam_list = './output/*.sam'

        Channel
        .fromPath(bam_list, checkIfExists: true)
        .set {bam_ch}


        VIEW(bam_ch)
        SORT(VIEW.out)
        FLAGSTAT(bam_ch)
        
}



/*
#==============================================
params
#==============================================
*/

params.reference = ""
params.indexType = ""

/*
#==============================================
Import processes
#==============================================
*/

include { BWA_INDEX } from './modules/3-index.nf'

/*
#==============================================
BWA-Indexing
#==============================================
*/

workflow INDEX {

	if (!params.reference || !params.indexType) {
                exit 1, "Missing --reference <*.fasta,*.fa> or --indexType <is,bwtsw>"
        } else {
                ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
                BWA_INDEX(ref_ch)
        }
}

workflow {

  INDEX()
}

workflow.onComplete {

    println     """
                Pipeline execution summary
                ---------------------------
                Completed at: ${workflow.complete}
                Duration    : ${workflow.duration}
                Success     : ${workflow.success}
                workDir     : ${workflow.workDir}
                exit status : ${workflow.exitStatus}
                """
}

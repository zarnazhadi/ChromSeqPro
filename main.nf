/*
#==============================================
params
#==============================================
*/

params.help = ""
params.reference = ""
params.indexType = ""
params.skipQC = false
params.skipIDXSTAT = false

if (params.help){

	log.info"""
	Usage:
    	
	The typical command for running the pipeline is as follows:

	nextflow run main.nf -profile hpc,singularity --reads '/input/*_R{1,2}*.fastq' --outdir '/project/'
    	./nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile standard,docker
    
	Required arguments:
    
	Performance options:
    
	Input File options:
    
	Save options:
    
	References		     If not specified in the configuration file or you wish to overwrite any of the references.
    
	QC Options:
    	
	""".stripIndent()

}

/*
#==============================================
Import processes
#==============================================
*/

include { RENAME 			} from './modules/0-registration'
include { FASTQC ; FASTQC as TRIMMEDQC	} from './modules/1-fastqc.nf'
include { CUTADAPT			} from './modules/2-cutadapt.nf'
//include { BWA_INDEX			} from './modules/3-index.nf'
include { BWA_MEM			} from './modules/4-alignment.nf'
//include { SORT; SAM_INDEX; IDXSTAT; VIEW; FLAGSTAT } from './modules/5-samtoo$
//include { PICARD                        } from './modules/6-picard.nf'

/*
#==============================================
Subworkflows
#==============================================
*/

/*
#==============================================
Registration
#==============================================
*/

/*
#==============================================
FastQC
#==============================================
*/

/*
#==============================================
Trimming
#==============================================
*/

/*
#==============================================
Alignment
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

/*
#==============================================
Samtools summaries
#==============================================
*/

/*
#==============================================
Remove duplicates
#==============================================
*/

/*
#==============================================
Main workflow
#==============================================
*/

workflow {

	if (!params.sampleList) {
        	exit 1, "Please specify sample_list.csv using --sampleList"
	} else {
		Channel.fromPath(params.sampleList) \
            	| splitCsv(header:true, sep:",") \
            	| map { row-> tuple(row.SAMPLEID, row.EXPERIMENT, row.CTYPE, file(row.R1), file(row.R2)) } \
            	| RENAME
		
		if ( params.skipQC ) {
        		CUTADAPT(RENAME.out)
        	} else {
        		FASTQC(RENAME.out)
        		CUTADAPT(RENAME.out)
        		TRIMMEDQC(CUTADAPT.out)
        	}
	} 

	BWA_MEM(CUTADAPT.out)
	
	}

workflow.onComplete {

    println 	"""
        	Pipeline execution summary
        	---------------------------
        	Completed at: ${workflow.complete}
        	Duration    : ${workflow.duration}
        	Success     : ${workflow.success}
        	workDir     : ${workflow.workDir}
        	exit status : ${workflow.exitStatus}
	        """
}

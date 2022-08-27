/*
#==============================================
params
#==============================================
*/

params.skipQC = false

/*
#==============================================
Import processes
#==============================================
*/

include { RENAME                        } from './modules/0-registration'
include { FASTQC ; FASTQC as TRIMMEDQC  } from './modules/1-fastqc.nf'
include { CUTADAPT                      } from './modules/2-cutadapt.nf'

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
}

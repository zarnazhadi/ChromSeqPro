process REGION_COUNTS {

	publishDir "${params.WINDOW_DIR}/", mode: 'copy'

	module 'R/3.6.1'

	input:
	path(INFILE)

	output:
	tuple val("${SAMPLEID}"), path("${SAMPLEID}-genome.windows"), path("${SAMPLEID}-promoters.windows")

	script:

	SAMPLEID=INFILE.name.split("-").first()	

	"""
	mkdir -p ${params.WINDOW_DIR}

        region-counts -ve -b ${params.BLACKLIST} ${params.GENOME} ${params.WINDOWSIZE} ${params.BAM_DIR}/${INFILE} > ${SAMPLEID}-promoters.windows
	region-counts -ve -b ${params.BLACKLIST} ${params.GENOME} ${params.WINDOWSIZE} ${params.BAM_DIR}/${INFILE} > ${SAMPLEID}-genome.windows
	"""
}

process PLOT_COUNTS {

	module 'R/3.6.1'

	input:
	tuple val(SAMPLEID), path(GENFILE), path(PROFILE)

	script:
	"""
	mkdir -p ${params.RES_DIR}/ 

        plot-counts -vt -c "grey25" -o ${params.RES_DIR}/${SAMPLEID}/plots/ -p "${SAMPLEID}-genome-" ${params.GENOME} ${GENFILE} 
        plot-counts -vt -c "orange" -o ${params.RES_DIR}/${SAMPLEID}/plots/ -p "${SAMPLEID}-promoters-" ${params.GENOME} ${PROFILE} 
	"""

}

process PROMOTER_SIGNAL {

	module 'R/3.6.1'

	input:
	tuple val(SAMPLEID), path(GENFILE), path(PROFILE)	

	script:

	CONTROL="control" 
	CTYPE=SAMPLEID.split("_").last()
	"""
	mkdir -p ${params.RES_DIR}/${SAMPLEID} 
	export LAMBDA_GCONTROL=\$(calculate-lambda -t0.5 ${params.GENOME} ${params.WINDOW_DIR}/genome/${CONTROL}_${CTYPE}-genome.windows)
	export LAMBDA_GEXP=\$(calculate-lambda -t0.5 ${params.GENOME} ${GENFILE}) 
	calculate-promoter-signal -v --alpha 0.00001 --lambda-gcontrol \${LAMBDA_GCONTROL} --lambda-gexp \${LAMBDA_GEXP} ${params.GENOME} ${params.WINDOW_DIR}/promoters/${CONTROL}_${CTYPE}-promoters.windows ${PROFILE} > ${params.RES_DIR}/${SAMPLEID}-signal.txt 
	""" 

}

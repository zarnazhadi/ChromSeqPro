process MACS2 {

        publishDir("${params.OUT_DIR}/07_MACS2/", mode: 'copy')

        input:
        tuple path(bamReads),
              path(ctrlReads)

	output:
	path("${bamName}/*")

        script:
	bamName = bamReads.toString().split("-")[0]

        """
	mkdir -p ${params.OUT_DIR}/07_MACS2/${bamName}
	macs2 callpeak -t ${bamReads} -c ${ctrlReads} \
        --outdir ${bamName} \ 
        """
}

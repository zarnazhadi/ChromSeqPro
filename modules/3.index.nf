process BWA_INDEX {
    publishDir "${params.INDEX_DIR}", mode: "copy"

    input:
      path reference

    output:
      tuple val("${reference}"), file("*.{fa,amb,ann,bwt,pac,sa}")

    script:
    if (params.indexType == "is")
        """
	bwa index -a is $reference
        """
    else (params.indexType == "bwtsw")
        """
	bwa index -a bwtsw $reference
        """
}

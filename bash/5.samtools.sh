# OPTIONS FOR GRID ENGINE ==============================================================
#
#$ -l h_rt=10:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE=================================================================

export ALIGNED_DIR="/nobackup/bs21zh/ZH_prelim/bash/aligned"

for bam in ${ALIGNED_DIR}/*.bam
do
	fileid=$(basename ${bam} .bam) # getting file name only
	
	out=$( echo ${ALIGNED_DIR}/${fileid}.sorted.bam) # new file name after alignment
	echo $out
	flag=$( echo ${FLAGSTAT_DIR}/${fileid}.txt)
	echo $flag
	/nobackup/leedsomics_tools/samtools-1.9/samtools sort -m 5G -o $out $bam
	/nobackup/leedsomics_tools/samtools-1.9/samtools flagstat $out > $flag

done

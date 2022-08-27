# OPTIONS FOR GRID ENGINE ==============================================================
#
#$ -l h_rt=10:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE=================================================================

export TRIMMED_DIR="/nobackup/bs21zh/ZH_prelim/bash/trimmed"
export ALIGNED_DIR="/nobackup/bs21zh/ZH_prelim/bash/aligned"

mkdir -p ${ALIGNED_DIR}

cd ../index

for r1 in ${TRIMMED_DIR}/*R1.fastq.gz
do
	r2=${r1/R1/R2}
	fileid=$(basename $r1 _R1.fastq.gz) # getting file name only
	
	out=$( echo ${ALIGNED_DIR}/${fileid}.bam) # new file name after alignment
	echo $out
  
  # BWA
	/nobackup/leedsomics_tools/bwa-0.7.17/bwa mem -P -M -t 8 GRCh38.p13.genome.fa.gz $r1 $r2 | /nobackup/leedsomics_tools/samtools-1.9/samtools sort -m 5G -o $out 

done

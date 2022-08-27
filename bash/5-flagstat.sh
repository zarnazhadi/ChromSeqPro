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
export FLAGSTAT_DIR="/nobackup/bs21zh/ZH_prelim/bash/flag_output"

mkdir -p ${FLAGSTAT_DIR}

for bam in ${ALIGNED_DIR}/*.bam
do
  	fileid=$(basename ${bam} .bam ) # getting file name only
        out=$( echo ${FLAGSTAT_DIR}/${fileid}.txt) # new file name 

        # flagstat
        /nobackup/leedsomics_tools/samtools-1.9/samtools flagstat $bam > $out
done

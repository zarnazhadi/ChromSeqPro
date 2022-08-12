# OPTIONS FOR GRID ENGINE ==============================================================
#
#$ -l h_rt=2:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE=================================================================

mkdir index
cd index

# download ref genome

/nobackup/leedsomics_tools/bwa-0.7.17/bwa index -p reference -a bwtsw gencode.v41.transcripts.fa.gz

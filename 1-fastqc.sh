# OPTIONS FOR GRID ENGINE ===================================================
#
#$ -l h_rt=2:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE====================================================

mkdir -p fastqc_files
/nobackup/leedsomics_tools/FastQC/fastqc \
../data/raw_files/*fastq.gz -threads 4 --outdir fastqc_files

# OPTIONS FOR GRID ENGINE ==============================================================
#
#$ -l h_rt=2:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE=================================================================

mkdir trimmed

for r1 in ../data/raw_files/IMPUT_R_R1.fastq.gz; do
# replace R1 in the file name with R2
    r2=${r1/R1/R2}
    out1=trimmed/$(basename ${r1})
    out2=${out1/R1/R2}
    echo $out1
    echo $out2
    /nobackup/leedsomics_tools/cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o $out1 -p $out2 $r1 $r2 --minimum-length=20 --overlap 3 -q 20 -j 24
done

# OPTIONS FOR GRID ENGINE ==============================================================
#
#$ -l h_rt=2:00:00
#$ -l h_vmem=16G
#$ -pe smp 8
#$ -cwd
#$ -j y
#$ -V
# OPTIONS FOR GRID ENGINE=================================================================

module load java/8u172 

mkdir -p unfiltered

for files in aligned/*.sorted.bam
do
	fileid=$(basename $files .sorted.bam ) # getting file name only
	outname=filtered/$( echo ${fileid}-marked.bam ) # new file name after trim
	metrics=filtered/$( echo ${fileid}-marked.txt )
	echo $outname
	echo $metrics
	# remove duplicates
	java -jar ../picard.jar MarkDuplicates --INPUT $files \
		--OUTPUT $outname --METRICS_FILE $metrics
done

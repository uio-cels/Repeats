#!/bin/sh
#SBATCH --job-name=worker
#SBATCH --account=nn9244k
#SBATCH --time=1:0:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=16

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi


REPO=$1
GENOME=$2
tPSI=$3
CDHIT=$4

module load genometools/1.5.7 
module load repeatmodeler/1.0.8
module load blast/2.2.26
module load perlmodules/5.10_2
module load hmmer/3.0
module load muscle/3.8.31
module load python2/2.7.9
module load blast+/2.2.29
module load repeatmasker/4.0.5
module load bedtools/2.17.0
module load cd-hit/4.6.4

if [ $TASK_ID == 1 ]; then

#-----------------------------------RepeatModeler----------------------------------------#

#	Making sure RepeatModeler chooses to sample from contigs larger than N50 by producing
#	a .csv file of sequence sizes with the N50 value. The file "N50" contains the N50 value,
#	which is used by a python script to print only scaffolds > N50. This makes sense if 
#	the genome FASTA consists of a high amount of small scaffolds, like the genomes in the
#	100 teleost analysis.

	assemblathon_stats.pl -csv $GENOME | grep "N50 scaffold length" | column -t | \
	awk '{print $4}' > N50

	python2 $REPO/N50.py -i $GENOME --N50 N50 > $GENOME.N50

#	Create index for RepeatModeler

	BuildDatabase -name $GENOME.N50 -engine ncbi $GENOME.N50

#	Run RepeatModeler

	RepeatModeler -database $GENOME.N50 -engine ncbi -pa 10

#################################
#	Removal of excess files  

	rm N50				
	rm $GENOME.N50			
#################################
	
elif [ $TASK_ID == 2 ]; then

#---------------------------------TransposonPSI------------------------------------------#
	
	perl $tPSI $GENOME nuc >tPSI.log

elif [ $TASK_ID = 3 ]; then

#-------------------------------LTRharvest/LTRdigest-------------------------------------#
	
#	Create indexes for LTRharvest/LTRdigest
	
	gt suffixerator -db $GENOME -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna

#-------------------------------------LTRharvest-----------------------------------------#

#	Detecting putative LTR retrotransposons with:
#	a minimum length of 100 bp
#	a maximum length of 6000 bp
#	99% similar long terminal repeats
#	An intact TG..CA motif

	gt ltrharvest -index $GENOME -out $GENOME.retrotransposons.out99 \
	-outinner $GENOME.retrotransposons.outinner99 -gff3 $GENOME.retrotransposons.gff99 \
	-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd \
	5 -motif tgca -similar 99 -vic 10 \
		> $GENOME.retrotransposons.result99 &

#	Detecting putative LTR retrotransposons with:
#	a minimum length of 1500 bp
#	a maximum length of 25000 bp
#	85% similar long terminal repeats

	gt ltrharvest -index $GENOME -out $GENOME.retrotransposons.out85 \
	-outinner $GENOME.retrotransposons.outinner85 -gff3 $GENOME.retrotransposons.gff85 \
	-minlenltr 100 -maxlenltr 6000 \
	-mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10 \
		> $GENOME.retrotransposons.result85 &

#	Detecting putative terminal-repeat retrotransposons in miniature (TRIM) with:
#	a minimum length of 280
#	a maximum length of 1500
#	99% similar long terminal repeats

	gt ltrharvest -index $GENOME -out $GENOME.TRIM.outT99 -outinner \
	$GENOME.TRIM.outinnerT99 -gff3 $GENOME.TRIM.gffT99 -minlenltr 70 -maxlenltr 500 \
	-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 \
		> $GENOME.TRIM.resultT99 &

#	Detecting putative terminal-repeat retrotransposons in miniature (TRIM) with:
#	a minimum length of 280
#	a maximum length of 1500
#	85% similar long terminal repeats

	gt ltrharvest -index $GENOME -out $GENOME.TRIM.outT85 -outinner \
	$GENOME.TRIM.outinnerT85 -gff3 $GENOME.TRIM.gffT85 -minlenltr 70 -maxlenltr 500 \
	-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -vic 10 \
		> $GENOME.TRIM.resultT85 &

	wait

#-------------------------------------LTRdigest------------------------------------------#

#	Sorting the GFF files for LTRdigest

	gt gff3 -sort $GENOME.retrotransposons.gff99 > $GENOME.retrotransposons.gff99.sort &
	gt gff3 -sort $GENOME.retrotransposons.gff85 > $GENOME.retrotransposons.gff85.sort &
	gt gff3 -sort $GENOME.TRIM.gffT85 > $GENOME.TRIM.gffT85.sort &
	gt gff3 -sort $GENOME.TRIM.gffT99 > $GENOME.TRIM.gffT99.sort &

	wait

#	Running LTRdigest, detecting elements with a primer binding site and retrotransposon
# 	specific enzymes for each element found by LTRharvest.

	gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
	$GENOME.retrotransposons.gff99.sort \
	$GENOME > $GENOME.retrotransposons.gff99.dgt &

	gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
	$GENOME.retrotransposons.gff85.sort \
	$GENOME > $GENOME.retrotransposons.gff85.dgt &

	gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
	$GENOME.TRIM.gffT85.sort \
	$GENOME > $GENOME.TRIM.gffT85.dgt &

	gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
	$GENOME.TRIM.gffT99.sort \
	$GENOME > $GENOME.TRIM.gffT99.dgt &

	wait

fi

#!/bin/sh
#SBATCH --job-name=ltrharvest/ltrdigest
#SBATCH --account=nn9244k
#SBATCH --time=0:10:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=1

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi
# Mail at beginning/end/on suspension
#SBATCH --mail-user=williambrynildsen@gmail.com
#SBATCH --mail-type=ALL

GENOME=$1

#-------------------------------------Indexing-------------------------------------------#

#Create indexes for LTRharvest/LTRdigest
gt suffixerator -db $GENOME -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna

#-------------------------------------LTRharvest-----------------------------------------#

##Detecting putative LTR retrotransposons with:
#	a minimum length of 100 bp
#	a maximum length of 6000 bp
#	99% similar long terminal repeats
#	An intact TG..CA motif

gt ltrharvest -index $GENOME -out $GENOME.retrotransposons.out99 \
-outinner $GENOME.retrotransposons.outinner99 -gff3 $GENOME.retrotransposons.gff99 \
-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 \
-motif tgca -similar 99 -vic 10 \
	> $GENOME.retrotransposons.result99

##Detecting putative LTR retrotransposons with:
#	a minimum length of 1500 bp
#	a maximum length of 25000 bp
#	85% similar long terminal repeats

gt ltrharvest -index $GENOME -out $GENOME.retrotransposons.out85 \
-outinner $GENOME.retrotransposons.outinner85 -gff3 $GENOME.retrotransposons.gff85 \
-minlenltr 100 -maxlenltr 6000 \
-mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10 \
	> $GENOME.retrotransposons.result85

##Detecting putative terminal-repeat retrotransposons in miniature (TRIM) with:
#	a minimum length of 280
#	a maximum length of 1500
#	99% similar long terminal repeats

gt ltrharvest -index $GENOME -out $GENOME.TRIM.outT99 -outinner \
$GENOME.TRIM.outinnerT99 -gff3 $GENOME.TRIM.gffT99 -minlenltr 70 -maxlenltr 500 \
-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 \
	> $GENOME.TRIM.resultT99

##Detecting putative terminal-repeat retrotransposons in miniature (TRIM) with:
#	a minimum length of 280
#	a maximum length of 1500
#	85% similar long terminal repeats

gt ltrharvest -index $GENOME -out $GENOME.TRIM.outT85 -outinner \
$GENOME.TRIM.outinnerT85 -gff3 $GENOME.TRIM.gffT85 -minlenltr 70 -maxlenltr 500 \
-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -vic 10 \
	> $GENOME.TRIM.resultT85

#-------------------------------------LTRdigest------------------------------------------#

#Preparing the GFF files for LTRdigest

gt gff3 -sort $GENOME.retrotransposons.gff99 > $GENOME.retrotransposons.gff99.sort
gt gff3 -sort $GENOME.retrotransposons.gff85 > $GENOME.retrotransposons.gff85.sort
gt gff3 -sort $GENOME.TRIM.gffT85 > $GENOME.TRIM.gffT85.sort
gt gff3 -sort $GENOME.TRIM.gffT99 > $GENOME.TRIM.gffT99.sort 

##Running LTRdigest, detecting elements with a primer binding site and retrotransposon
# specific enzymes for each harvested group

gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
$GENOME.retrotransposons.gff99.sort \
$GENOME > $GENOME.retrotransposons.gff99.dgt

gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
$GENOME.retrotransposons.gff85.sort \
$GENOME > $GENOME.retrotransposons.gff85.dgt 

gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
$GENOME.TRIM.gffT85.sort \
$GENOME > $GENOME.TRIM.gffT85.dgt 

gt ltrdigest -trnas $REPO/eukaryotic-tRNAs.fa -hmms $REPO/gydb/*hmm -- \
$GENOME.TRIM.gffT99.sort \
$GENOME > $GENOME.TRIM.gffT99.dgt 

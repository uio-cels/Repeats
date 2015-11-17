#!/bin/sh
#SBATCH --job-name=pipeline
#SBATCH --account=nn9244k
#SBATCH --time=0:30:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=1

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi
# Mail at beginning/end/on suspension
#SBATCH --mail-user=williambrynildsen@gmail.com
#SBATCH --mail-type=ALL

#Genome headers must not contain any special characters (only letters and/or numbers).

export REPO=/work/users/willibr/testing/repository
export GENOME=/work/users/willibr/testing/c.elegans/yeast.fa
export tPSI=/projects/cees/bin/TransposonPSI/transposonPSI.pl
export CDHIT=/projects/cees/bin/cd-hit/cd-hit-v4.6.4-2015-0603/cd-hit-est

#$REPO/eukaryotic-tRNAs.fa
#$REPO/TRIM99.CRL_Step1.pl
#$REPO/TRIM85.CRL_Step1.pl
#$REPOretro99.CRL_Step1.pl
#$REPO/retro85.CRL_Step1.pl
#$REPO/CRL_Step2.pl
#$REPO/CRL_Step3.pl
#$REPO/filter_protein_match.lua
#$REPO/uniprot_sprot.fasta

module load genometools/1.5.7 
module load repeatmodeler/1.0.8
module load blast/2.2.26
module load perlmodules/5.10_2
module load hmmer/3.0
module load muscle/3.8.31
module load python2/2.7.9
module load usearch/7.0.1090
module load blast+/2.2.29
module load repeatmasker/4.0.5
module load bedtools/2.17.0
module load cd-hit/4.6.4

#-------------------------------------Indexing-------------------------------------------#

#Create indexes for LTRharvest/LTRdigest
gt suffixerator -db $GENOME -indexname $GENOME -tis -suf -lcp -des -ssp -sds -dna
#Create index for RepeatModeler
BuildDatabase -name $GENOME  -engine ncbi $GENOME 

#-----------------------------------RepeatModeler----------------------------------------#

RepeatModeler -database $GENOME -engine ncbi -pa 10

#Create indexes for BLASTX

makeblastdb -in $REPO/uniprot_sprot.fasta -dbtype prot
makeblastdb -in $REPO/RepeatPeps.lib -dbtype prot

#------------------------------------TransposonPSI---------------------------------------#

perl $tPSI $GENOME nuc >tPSI.log
 
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

#------------------------------Filtration and extraction---------------------------------#
##Running custom scripts provided by Ning Jiang (she will publish them soon) to filter
# the LTRdigest results

perl $REPO/retro99.CRL_Step1.pl -gff $GENOME.retrotransposons.gff99.dgt
perl $REPO/retro85.CRL_Step1.pl -gff $GENOME.retrotransposons.gff85.dgt 
perl $REPO/TRIM85.CRL_Step1.pl -gff $GENOME.TRIM.gffT85.dgt 
perl $REPO/TRIM99.CRL_Step1.pl -gff $GENOME.TRIM.gffT99.dgt 

#Running a lua script to select the LTR retrotransposons with specific enzymes

gt select -rule_files $REPO/filter_protein_match.lua -- \
< $GENOME.retrotransposons.gff99.dgt \
	> $GENOME.retrotransposons.gff99.dgt.withdomains

gt select -rule_files $REPO/filter_protein_match.lua -- \
< $GENOME.retrotransposons.gff85.dgt \
	> $GENOME.retrotransposons.gff85.dgt.withdomains

##Continuing running scripts provided by Ning Jiang.
#
#First, creating folders
mkdir $GENOME.TRIM85.fasta_files
mkdir $GENOME.TRIM99.fasta_files
mkdir $GENOME.retro99.fasta_files
mkdir $GENOME.retro85.fasta_files

#Copying the last script into each folder (mandatory)
cp $REPO/CRL_Step3.pl $GENOME.TRIM85.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.TRIM99.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.retro99.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.retro85.fasta_files

#Copying other files as well, including the genome file (not good for space)

cp retro85.CRL_Step1_Passed_Elements.txt $GENOME.retro85.fasta_files/
cp $GENOME.retrotransposons.out85 $GENOME.retro85.fasta_files/
cp $GENOME.retrotransposons.result85 $GENOME.retro85.fasta_files/
cp $GENOME $GENOME.retro85.fasta_files/

cp retro99.CRL_Step1_Passed_Elements.txt $GENOME.retro99.fasta_files/
cp $GENOME.retrotransposons.out99 $GENOME.retro99.fasta_files/
cp $GENOME.retrotransposons.result99 $GENOME.retro99.fasta_files/
cp $GENOME $GENOME.retro99.fasta_files/

cp TRIM99.CRL_Step1_Passed_Elements.txt $GENOME.TRIM99.fasta_files/
cp $GENOME.TRIM.outT99 $GENOME.TRIM99.fasta_files/
cp $GENOME.TRIM.resultT99 $GENOME.TRIM99.fasta_files/
cp $GENOME $GENOME.TRIM99.fasta_files/

cp TRIM85.CRL_Step1_Passed_Elements.txt $GENOME.TRIM85.fasta_files/
cp $GENOME.TRIM.outT85 $GENOME.TRIM85.fasta_files/
cp $GENOME.TRIM.resultT85 $GENOME.TRIM85.fasta_files/
cp $GENOME $GENOME.TRIM85.fasta_files/

##Running part II of the custom script, and moving output FASTA files into the right
# folder

perl $REPO/CRL_Step2.pl --step1 retro85.CRL_Step1_Passed_Elements.txt --repeatfile \
$GENOME.retrotransposons.out85 --resultfile $GENOME.retrotransposons.result85 \
--sequencefile $GENOME --removed_repeats $GENOME.retro85.CRL_Step2_Passed_Elements.fasta

mv Repeat_* $GENOME.retro85.fasta_files/

perl $REPO/CRL_Step2.pl --step1 retro99.CRL_Step1_Passed_Elements.txt --repeatfile \
$GENOME.retrotransposons.out99 --resultfile $GENOME.retrotransposons.result99 \
--sequencefile $GENOME --removed_repeats $GENOME.retro99.CRL_Step2_Passed_Elements.fasta
 
mv Repeat_* $GENOME.retro99.fasta_files/
 
perl $REPO/CRL_Step2.pl --step1 TRIM99.CRL_Step1_Passed_Elements.txt --repeatfile \
$GENOME.TRIM.outT99 --resultfile $GENOME.TRIM.resultT99  --sequencefile $GENOME \
--removed_repeats $GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta
 
mv Repeat_* $GENOME.TRIM99.fasta_files/

perl $REPO/CRL_Step2.pl --step1 TRIM85.CRL_Step1_Passed_Elements.txt --repeatfile \
$GENOME.TRIM.outT85 --resultfile $GENOME.TRIM.resultT85  --sequencefile $GENOME \
--removed_repeats $GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta
 
mv Repeat_* $GENOME.TRIM85.fasta_files/

##Here, part III of the custom scripts need to wait until all Repeat_* files are moved,
# so I add a little nap. But I need a better solution.

sleep 60

cd $GENOME.retro85.fasta_files/

perl CRL_Step3.pl --directory . --step2  \
$GENOME.retro85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25


cd  ../$GENOME.TRIM85.fasta_files/

perl CRL_Step3.pl --directory . --step2  \
$GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25

cd ../$GENOME.TRIM99.fasta_files/

perl CRL_Step3.pl --directory . --step2  \
$GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25

cd $GENOME.retro99.fasta_files/

perl CRL_Step3.pl --directory . --step2  \
$GENOME.retro99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25


#Changing the headers of the genome FASTA file before extraction

cd ..

python2 $REPO/change_headers_to_seqN.py -i $GENOME | sed 's/ //g' \
	> $GENOME.changed_headers

#Getting the full positions of the elements out of the GFF

grep "ID=repeat_region" $GENOME.retrotransposons.gff85.dgt.withdomains \
	> $GENOME.retrotransposons.gff85.dgt.withdomains.full
	
grep "ID=repeat_region" $GENOME.retrotransposons.gff99.dgt.withdomains \
	> $GENOME.retrotransposons.gff99.dgt.withdomains.full

#Extracting retrotransposon sequences with proteins
bedtools getfasta -fi $GENOME.changed_headers -bed \
$GENOME.retrotransposons.gff85.dgt.withdomains.full -fo \
$GENOME.retrotransposons.gff85.dgt.withdomains.full.fasta

bedtools getfasta -fi $GENOME.changed_headers -bed \
$GENOME.retrotransposons.gff99.dgt.withdomains.full -fo \
$GENOME.retrotransposons.gff99.dgt.withdomains.full.fasta

#In order to extract sequences from the TransposonPSI file, files had to be altered
cut -f 1,4,5,9 $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | sed 's/ID=.*Target.//g' | \
sed 's/;.*//g' | column -t | awk '{print $4}' \
	> tPSI.classes

paste tPSI.classes $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | column -s $'\t' -t \
	> tmp4

awk '{print $2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' tmp4 \
	> tmp5

sed 's/ /\t/g' tmp5 \
	> tmp6

#The actual extraction done using bedtools (getfasta)

bedtools getfasta -name -fi $GENOME -bed tmp6 -fo $GENOME.tPSI.fasta

#The TransposonPSI library is classified according to protein homology
python2 $REPO/reprint.tPSI.lib.py -i $GENOME.tPSI.fasta | fold -w 60 \
	> $GENOME.tPSI.classified.fasta

#--------------------------------Merging of libraries------------------------------------#
#

cp $GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta $GENOME.retro99.fasta
cp $GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta $GENOME.retro85.fasta
cp $GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta $GENOME.TRIM99.fasta
cp $GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta $GENOME.TRIM85.fasta

#The python script renames the elements according to their group

python2 $REPO/reprint.ltrharvest.lib.py -i $GENOME.retro99.fasta --name retro99 \
--type Unknown > $GENOME.retro99.renamed.fasta

python2 $REPO/reprint.ltrharvest.lib.py -i $GENOME.retro85.fasta --name retro85 \
--type Unknown > $GENOME.retro85.renamed.fasta

python2 $REPO/reprint.ltrharvest.lib.py -i $GENOME.TRIM99.fasta --name TRIM99 \
--type TRIM > $GENOME.TRIM99.renamed.fasta

python2 $REPO/reprint.ltrharvest.lib.py -i $GENOME.TRIM85.fasta --name TRIM85 \
--type TRIM > $GENOME.TRIM85.renamed.fasta

python2 $REPO/reprint.ltrharvest.lib.py -i \
$GENOME.retrotransposons.gff99.dgt.withdomains.full.fasta --name retro99withdomains \
--type Unknown > $GENOME.retrotransposons.gff99.dgt.withdomains.full.renamed.fasta

python2 $REPO/reprint.ltrharvest.lib.py -i \
$GENOME.retrotransposons.gff85.dgt.withdomains.full.fasta --name retro85withdomains \
--type Unknown > $GENOME.retrotransposons.gff85.dgt.withdomains.full.renamed.fasta

#Merging of libraries

cat $GENOME.retro99.renamed.fasta \
$GENOME.retro85.renamed.fasta \
$GENOME.TRIM99.renamed.fasta \
$GENOME.TRIM85.renamed.fasta \
$GENOME.retrotransposons.gff99.dgt.withdomains.full.renamed.fasta \
$GENOME.retrotransposons.gff85.dgt.withdomains.full.renamed.fasta \
> $GENOME.retrotransposons.TRIMs.fasta
cd RM*
sed 's/ .*//g' consensi.fa.classified > consensi.fa.classified.stripped
cd ..
cat $GENOME.retrotransposons.TRIMs.fasta \
RM*/consensi.fa.classified.stripped \
$GENOME.tPSI.classified.fasta | fold -w 60 > $GENOME.repeats.fasta

#---------------------------------Clustering---------------------------------------------#

$CDHIT -i $GENOME.repeats.fasta -o $GENOME.repeats.srt.nr

#---------------------------Removing false positives-------------------------------------#

#Removing false transposons by checking against SwissProt-UniProt

blastx -query $GENOME.repeats.srt.nr \
-db $REPO/uniprot_sprot.fasta -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50

blastx -query  $GENOME.repeats.srt.nr \
-db $REPO/RepeatPeps.lib -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50

sort -k1,1 -k12,12nr -k11,11n \
$GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring

sort -k1,1 -k12,12nr -k11,11n \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring

cut -f 1 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring \
> $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column

cut -f 1 $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring \
> $GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column

echo "Num entries only hit to RepeatPeps" 
comm -13  $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column 
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column \
	> $GENOME.repeats.only_repeatpeps
	
cat $GENOME.repeats.only_repeatpeps |wc -l

echo "Num entries only hit to UniProt"

comm -23 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column \
	> $GENOME.repeats.only_uniprot
	
cat $GENOME.repeats.only_uniprot |wc -l

echo "Num entries hit to both"

comm -12 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
$GENOME.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column | wc -l

python2 $REPO/reprint.filtered.lib.py -i $GENOME.repeats.srt.nr -l \
$GENOME.repeats.only_uniprot | fold -w 60 \
	> $GENOME.repeats.srt.nr.no_uniprot

#-----------------------------------------RepeatMasker-----------------------------------#

##Creating a total repeat library

cat $GENOME.repeats.srt.nr.no_uniprot $REPO/repbase.update.04102015.teleostei.lib \
	> $GENOME.total.repeat.library

#Making folders for each library

mkdir $GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output
mkdir $GENOME.total.repeat.library.repeatmasker_output
mkdir repbase.update.04102015.teleostei.lib.repeatmasker_output

RepeatMasker -lib $GENOME.repeats.srt.nr.no_uniprot -a -s -dir \
$GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output/ $GENOME

RepeatMasker -lib $GENOME.total.repeat.library -a -s -dir \
$GENOME.total.repeat.library.repeatmasker_output/ $GENOME

RepeatMasker -lib $REPO/repbase.update.04102015.teleostei.lib -a -s -dir \
repbase.update.04102015.teleostei.lib.repeatmasker_output/ $GENOME

#summary.landscape.slurm

/cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
$GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output/*.out \
	> $GENOME.repmask.denovo.summary
	
/cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
$GENOME.total.repeat.library.repeatmasker_output/*.out \
	> $GENOME.repmask.total.summary
	
/cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
repbase.update.04102015.teleostei.lib.repeatmasker_output/*.out\
	> $GENOME.repmask.repbase.summary


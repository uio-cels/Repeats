#!/bin/sh
#SBATCH --job-name=pipeline
#SBATCH --account=nn9244k
#SBATCH --time=168:0:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=16

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi

module purge

TMP=$LOCALTMP
TMPDIR=/work/users/willibr/tmp

mkdir -p repeats

species=${1%}
origdir=$PWD

mkdir -p repeats/${species}
echo ${species}

mkdir ${TMPDIR}/${species}_repeats
cd ${TMPDIR}/${species}_repeats

cp $origdir/repeats_worker_script.slurm ${TMPDIR}/${species}_repeats
cp $origdir/$1 ${TMPDIR}/${species}_repeats

#------------------------------Pointing to paths-----------------------------------------#

#	The repository contains databases and scripts used during the pipeline.

export REPO=/projects/cees/in_progress/william/pipeline/repository

#---Content of the repository:

#	Databases:
#	eukaryotic-tRNAs.fa !
#	uniprot_sprot17nov2015.fasta
#	RepeatPeps.lib !
#	TEseeker.DNAtransposons.fasta !

#	Python scripts (only tested in python2/2.7.9):
#	change_headers_to_seqN.py (retrived from github)
#	reprint.tPSI.lib.py (retrived from github)
#	reprint.ltrharvest.lib.py (retrived from github)
#	reprint.ltrharvest.lib.new.py !
#	reprint.filtered.lib.py (retrived from github)
#	sequence_cleaner.py !
#	N50.py !
#	reprint.genome.py !

#	LTRharvest/LTRdigest filtering (courtesy of Ning Jiang):
#	TRIM99.CRL_Step1.pl
#	TRIM85.CRL_Step1.pl
#	retro99.CRL_Step1.pl
#	retro85.CRL_Step1.pl
#	CRL_Step2.pl 
#	CRL_Step3.pl 
#	CRL_Step4.pl !
#	CRL_Step5.pl ! 
#	ltr_library.pl !

#	RepeatMasker parsing:
#	get.rm.master.sh !
#	get.rm.stats.sh !
#	get.rm.stats.rename.sh !
#	get.rm.overview.sh !

#	Misc:
#	assemblathon_stats.pl !	
#	filter_protein_match.lua !

#	The only parameter given the pipeline is a FASTA file of the genome.

export GENOME=$1

#	Point to a modified version of TransposonPSI with dual core usage.

export tPSI=/projects/cees/bin/TransposonPSI/dualcore_transposonPSI.pl

#	Point to CD-HIT -- Can probably change this with using the module system

export CDHIT=/projects/cees/bin/cd-hit/cd-hit-v4.6.4-2015-0603/cd-hit-est

#-----------------------------Loading modules--------------------------------------------#

module load genometools/1.5.7 
module load repeatmodeler/1.0.8
module load blast/2.2.26
module load perlmodules/5.10_2
module load hmmer/3.1b1
module load muscle/3.8.31
module load python2/2.7.9
module load blast+/2.2.29
module load repeatmasker/4.0.6
module load bedtools/2.17.0
module load cd-hit/4.6.4

# #------------------------------Genome index table----------------------------------------#
# #	NB! Genome FASTA headers must not contain any special characters (only letters and/or 
# #	numbers).
# 
#	Renaming headers (to SEQ0, SEQ1, SEQ2 etc.)
python2 $REPO/reprint.genome.py -i $GENOME > $GENOME.seq

#	Printing headers of the original genome file
grep ">" $GENOME > halfoftable1

#	Printing headers of the renamed genome file
grep ">" $GENOME.seq > halfoftable2

#	Creating a file with two columns containing original header name and new header name
paste halfoftable1 halfoftable2 > $GENOME.table

#	Pointing the main variable, $GENOME to the version with renamed headers.

GENOME=$GENOME.seq

#---------------------------------Harvesting repeats-------------------------------------#

#	Starting LTRharvest/LTRdigest, RepeatModeler and TransposonPSI

arrayrun 1-3 repeats_worker_script.slurm $REPO $GENOME $tPSI $CDHIT 

wait

#------------------------------Filtration and extraction---------------------------------#

#	Running custom scripts provided by Ning Jiang (she will publish them soon) to filter
# 	the LTRdigest results.

#	First, the putative retrotransposons are only kept if they contain a retrotransposon
#	specific domain.

gt select -rule_files $REPO/filter_protein_match.lua -- \
< $GENOME.retrotransposons.gff85.dgt \
	> $GENOME.retrotransposons.gff85.dgt.withdomains

gt select -rule_files $REPO/filter_protein_match.lua -- \
< $GENOME.retrotransposons.gff99.dgt \
	> $GENOME.retrotransposons.gff99.dgt.withdomains

#	Step1: The output of LTRdigest (.dgt files) are processed so that only elements with
#	a polypurine tract (PPT) and/or a primer binding site (PBS) are chosen, but only if 
#	at least 50% of the PPT/PBS sequence are located in the internal regions of the 
#	putative LTR retrotransposon and only if the distance between the long terminal 
#	repeats (LTRs) and PPT's/PBS's are maximum 20 bp.

perl $REPO/retro99.CRL_Step1.pl -gff $GENOME.retrotransposons.gff99.dgt.withdomains
perl $REPO/retro85.CRL_Step1.pl -gff $GENOME.retrotransposons.gff85.dgt.withdomains
perl $REPO/TRIM85.CRL_Step1.pl -gff $GENOME.TRIM.gffT85.dgt 
perl $REPO/TRIM99.CRL_Step1.pl -gff $GENOME.TRIM.gffT99.dgt 

#	OUTPUT: ....CRL_Step1_Passed_Elements.txt

#	Creating folders to place a lot of files

mkdir $GENOME.TRIM85.fasta_files
mkdir $GENOME.TRIM99.fasta_files
mkdir $GENOME.retro99.fasta_files
mkdir $GENOME.retro85.fasta_files

# 	Copying the third script into each folder (mandatory). CRL_Step3 needs to be run in
#	the same folder as the genome file and some other files.

cp $REPO/CRL_Step3.pl $GENOME.TRIM85.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.TRIM99.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.retro99.fasta_files
cp $REPO/CRL_Step3.pl $GENOME.retro85.fasta_files

#	Retro85

cp retro85.CRL_Step1_Passed_Elements.txt $GENOME.retro85.fasta_files/
cp $GENOME.retrotransposons.out85 $GENOME.retro85.fasta_files/
cp $GENOME.retrotransposons.result85 $GENOME.retro85.fasta_files/
cp $GENOME $GENOME.retro85.fasta_files/

#	Retro99

cp retro99.CRL_Step1_Passed_Elements.txt $GENOME.retro99.fasta_files/
cp $GENOME.retrotransposons.out99 $GENOME.retro99.fasta_files/
cp $GENOME.retrotransposons.result99 $GENOME.retro99.fasta_files/
cp $GENOME $GENOME.retro99.fasta_files/

#	TRIM99

cp TRIM99.CRL_Step1_Passed_Elements.txt $GENOME.TRIM99.fasta_files/
cp $GENOME.TRIM.outT99 $GENOME.TRIM99.fasta_files/
cp $GENOME.TRIM.resultT99 $GENOME.TRIM99.fasta_files/
cp $GENOME $GENOME.TRIM99.fasta_files/

#	TRIM85

cp TRIM85.CRL_Step1_Passed_Elements.txt $GENOME.TRIM85.fasta_files/
cp $GENOME.TRIM.outT85 $GENOME.TRIM85.fasta_files/
cp $GENOME.TRIM.resultT85 $GENOME.TRIM85.fasta_files/
cp $GENOME $GENOME.TRIM85.fasta_files/

#	The elements are filtered so that no element contain sequence gaps of more
#	than 50N's. 
#	
#	False elements might have LTR-like sequence outside their "LTRs", to check this, 
#	50 bp of flanking sequence are collected from the remaining elements.

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

#	OUTPUT: $GENOME....CRL_Step2_Passed_Elements.fasta

#	Files are copied.

cp $GENOME.retro85.CRL_Step2_Passed_Elements.fasta $GENOME.retro85.fasta_files/
cp $GENOME.retro99.CRL_Step2_Passed_Elements.fasta $GENOME.retro99.fasta_files/
cp $GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta $GENOME.TRIM85.fasta_files/
cp $GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta $GENOME.TRIM99.fasta_files/

#	Flanking sequences are aligned using MUSCLE, and elements with flanking regions that
#	are more than 60% similar are considered false positives.

cd $GENOME.retro85.fasta_files/

perl CRL_Step3.pl --directory . --step2 \
$GENOME.retro85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 

cd ../$GENOME.TRIM85.fasta_files/

perl CRL_Step3.pl --directory . --step2 \
$GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 

cd ../$GENOME.TRIM99.fasta_files/

perl CRL_Step3.pl --directory . --step2 \
$GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 

cd ../$GENOME.retro99.fasta_files/

perl CRL_Step3.pl --directory . --step2 \
$GENOME.retro99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25

#	OUTPUT: CRL_Step3_Passed_Elements.fasta

#	Detection of nested transposons. 
#
#	LTR sequences are used to mask internal 
#	regions of the elements (matches indicate nested LTRs). In addition, internal 
#	sequences are compared with a database of known non-LTR transposons (using 
#	RepeatMasker). 

cd ..

mkdir retro85
mkdir retro99
mkdir TRIM99
mkdir TRIM85

cd retro85

perl $REPO/ltr_library.pl --resultfile ../$GENOME.retrotransposons.result85 \
--step3 ../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
--sequencefile ../$GENOME

cp lLTR_Only.lib $GENOME.retro85.lLTR_Only.lib

cd ../retro99

perl $REPO/ltr_library.pl --resultfile ../$GENOME.retrotransposons.result99 \
--step3 ../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
--sequencefile ../$GENOME

cp lLTR_Only.lib $GENOME.retro99.lLTR_Only.lib

cd ../TRIM85

perl $REPO/ltr_library.pl --resultfile ../$GENOME.TRIM.resultT85 \
--step3 ../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
--sequencefile ../$GENOME

cp lLTR_Only.lib $GENOME.TRIM85.lLTR_Only.lib

cd ../TRIM99

perl $REPO/ltr_library.pl --resultfile ../$GENOME.TRIM.resultT99 \
--step3 ../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
--sequencefile ../$GENOME

cp lLTR_Only.lib $GENOME.TRIM99.lLTR_Only.lib

cd ..

cat retro85/$GENOME.retro85.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
	> retro85/$GENOME.repeats_to_mask_LTR85.fasta

cat retro99/$GENOME.retro99.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
	> retro99/$GENOME.repeats_to_mask_LTR99.fasta

cat TRIM85/$GENOME.TRIM85.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
	> TRIM85/$GENOME.repeats_to_mask_TRIM85.fasta
	
cat TRIM99/$GENOME.TRIM99.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
	> TRIM99/$GENOME.repeats_to_mask_TRIM99.fasta

RepeatMasker -lib retro85/$GENOME.repeats_to_mask_LTR85.fasta -pa 16 -nolow -no_is \
-dir . $GENOME.retrotransposons.outinner85 

RepeatMasker -lib retro99/$GENOME.repeats_to_mask_LTR99.fasta -pa 16 -nolow -no_is \
-dir . $GENOME.retrotransposons.outinner99

RepeatMasker -lib TRIM85/$GENOME.repeats_to_mask_TRIM85.fasta -pa 16 -nolow -no_is \
-dir . $GENOME.TRIM.outinnerT85

RepeatMasker -lib TRIM99/$GENOME.repeats_to_mask_TRIM99.fasta -pa 16 -nolow -no_is \
-dir . $GENOME.TRIM.outinnerT99

#	OUTPUT: ...outinnerXX.masked

#	Building exemplars
#
#	Representative sequences are chosen by an all vs. all BLASTN search. First sequences
#	are extracted.
	
MASKED1=$GENOME.retrotransposons.outinner85.masked
MASKED2=$GENOME.retrotransposons.outinner99.masked
MASKED3=$GENOME.TRIM.outinnerT85.masked
MASKED3=$GENOME.TRIM.outinnerT99.masked

cd retro85

if [ -f "$MASKED1" ];
then
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.retrotransposons.result85 --innerfile \
	../$GENOME.retrotransposons.outinner85.masked --sequencefile ../$GENOME
else
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.retrotransposons.result85 --innerfile \
	../$GENOME.retrotransposons.outinner85 --sequencefile ../$GENOME
fi

cd ../retro99

if [ -f "$MASKED2" ];
then
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.retrotransposons.result99 --innerfile \
	../$GENOME.retrotransposons.outinner99.masked --sequencefile ../$GENOME
else
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.retrotransposons.result99 --innerfile \
	../$GENOME.retrotransposons.outinner99 --sequencefile ../$GENOME
fi

cd ../TRIM99

if [ -f "$MASKED3" ];
then
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.TRIM.resultT99 --innerfile \
	../$GENOME.TRIM.outinnerT99.masked --sequencefile ../$GENOME
else
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.TRIM.resultT99 --innerfile \
	../$GENOME.TRIM.outinnerT99 --sequencefile ../$GENOME
fi

cd ../TRIM85

if [ -f "$MASKED4" ];
then
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.TRIM.resultT85 --innerfile \
	../$GENOME.TRIM.outinnerT85.masked --sequencefile ../$GENOME
else
	perl $REPO/CRL_Step4.pl --step3 \
	../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
	--resultfile ../$GENOME.TRIM.resultT85 --innerfile \
	../$GENOME.TRIM.outinnerT85 --sequencefile ../$GENOME
fi

cd ..

module load blast+/2.2.29

makeblastdb -in retro85/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in retro85/Inner_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in retro99/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in retro99/Inner_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in TRIM85/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in TRIM85/Inner_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in TRIM99/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
makeblastdb -in TRIM99/Inner_Seq_For_BLAST.fasta -dbtype nucl

blastn -query retro85/lLTRs_Seq_For_BLAST.fasta -db retro85/lLTRs_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out retro85/lLTRs_Seq_For_BLAST.fasta.out

blastn -query retro85/Inner_Seq_For_BLAST.fasta -db retro85/Inner_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out retro85/Inner_Seq_For_BLAST.fasta.out

blastn -query retro99/lLTRs_Seq_For_BLAST.fasta -db retro99/lLTRs_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out retro99/lLTRs_Seq_For_BLAST.fasta.out

blastn -query retro99/Inner_Seq_For_BLAST.fasta -db retro99/Inner_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out retro99/Inner_Seq_For_BLAST.fasta.out

blastn -query TRIM85/lLTRs_Seq_For_BLAST.fasta -db TRIM85/lLTRs_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out TRIM85/lLTRs_Seq_For_BLAST.fasta.out

blastn -query TRIM85/Inner_Seq_For_BLAST.fasta -db TRIM85/Inner_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out TRIM85/Inner_Seq_For_BLAST.fasta.out

blastn -query TRIM99/lLTRs_Seq_For_BLAST.fasta -db TRIM99/lLTRs_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out TRIM99/lLTRs_Seq_For_BLAST.fasta.out

blastn -query TRIM99/Inner_Seq_For_BLAST.fasta -db TRIM99/Inner_Seq_For_BLAST.fasta \
-evalue 1e-10 -num_descriptions 1000 -out TRIM99/Inner_Seq_For_BLAST.fasta.out

#	OUTPUT: lLTR_Seq_For_BLAST.fasta
#			Inner_Seq_For_BLAST.fasta

#	Selecting exemplars. In the end no elements should share 80% bp's over 90% of their 
#	length.

perl $REPO/CRL_Step5.pl --LTR_blast retro85/lLTRs_Seq_For_BLAST.fasta.out \
--inner_blast retro85/Inner_Seq_For_BLAST.fasta.out \
--step3 $GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
--final $GENOME.retro85.lib

perl $REPO/CRL_Step5.pl --LTR_blast retro99/lLTRs_Seq_For_BLAST.fasta.out \
--inner_blast retro99/Inner_Seq_For_BLAST.fasta.out \
--step3 $GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
--final $GENOME.retro99.lib

perl $REPO/CRL_Step5.pl --LTR_blast TRIM85/lLTRs_Seq_For_BLAST.fasta.out \
--inner_blast TRIM85/Inner_Seq_For_BLAST.fasta.out \
--step3 $GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
--final $GENOME.TRIM85.lib

perl $REPO/CRL_Step5.pl --LTR_blast TRIM99/lLTRs_Seq_For_BLAST.fasta.out \
--inner_blast TRIM99/Inner_Seq_For_BLAST.fasta.out \
--step3 $GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
--final $GENOME.TRIM99.lib

#	OUTPUT: ...lib

#	Renaming headers in each .lib file according to how they were generated.

python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.retro85.lib \
--name retrotransposon85 > $GENOME.retro85.renamed.lib

python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.retro99.lib \
--name retrotransposon99 > $GENOME.retro99.renamed.lib 

python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.TRIM85.lib \
--name TRIM85 > $GENOME.TRIM85.renamed.lib 

python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.TRIM99.lib \
--name TRIM99 > $GENOME.TRIM99.renamed.lib

#	Merging each group

cat $GENOME.TRIM99.renamed.lib $GENOME.TRIM85.renamed.lib $GENOME.retro99.renamed.lib \
$GENOME.retro85.renamed.lib \
	> $GENOME.retrotransposons.TRIMs.lib

#	Trying to classify the consensus sequences using RepeatClassifier

RepeatClassifier -consensi $GENOME.retrotransposons.TRIMs.lib -engine ncbi

#	Output is $GENOME.retrotransposons.TRIMs.lib.classified.
#	Elements classified as simple repeats, buffers or unknowns are reclassified into the
#	LTR superfamily using sed. 

sed 's/#Simple_repeat/#LTR\/Unknown/g' \
$GENOME.retrotransposons.TRIMs.lib.classified | \
sed 's/#buffer/#LTR\/Unknown/g' | \
sed 's/#Unknown/#LTR\/Unknown/g'| \
sed 's/#DNA.*/#LTR\/Unknown/g' | \
sed 's/LINE.*/#LTR\/Unknown/g' \
	> $GENOME.retrotransposons.TRIMs.lib.classified.reclassified

#-------------------------------tPSI FASTA extraction------------------------------------#

In order to extract sequences from the TransposonPSI file, files had to be altered
cut -f 1,4,5,9 $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | sed 's/ID=.*Target.//g' | \
sed 's/;.*//g' | column -t | awk '{print $4}' \
	> tPSI.classes

paste tPSI.classes $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | column -s $'\t' -t \
	> tmp4

awk '{print $2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' tmp4 \
	> tmp5

sed 's/ /\t/g' tmp5 \
	> tmp6

#	The actual extraction done using bedtools (getfasta)

bedtools getfasta -name -fi $GENOME -bed tmp6 -fo $GENOME.tPSI.fasta

#	The TransposonPSI library is classified according to protein homology
python2 $REPO/reprint.tPSI.lib.py -i $GENOME.tPSI.fasta | fold -w 60 \
	> $GENOME.tPSI.classified.fasta

#################################
#	Removal of excess files  

rm tmp* tPSI.classes *tPSI.fasta *TPSI.allHits.chains.gff3  tPSI.log
rm *TPSI.allHits.chains.bestPerLocus *TPSI.allHits.chains *TPSI.allHits
################################

---------------------------------Clustering---------------------------------------------#
Clustering of tPSI elements are done using CD-HIT-EST with 80% cutoff.

$CDHIT -i $GENOME.tPSI.classified.fasta -o $GENOME.tPSI.classified.fasta.nr.80.n5 \
-c 0.8 -T 16 -n 5

#--------------------------------Merging of libraries------------------------------------#

cd RM*
sed 's/ .*//g' consensi.fa.classified > consensi.fa.classified.stripped
cd ..
cat $GENOME.retrotransposons.TRIMs.lib.classified.reclassified \
RM*/consensi.fa.classified.stripped \
$GENOME.tPSI.classified.fasta.nr.80.n5 | fold -w 60 > $GENOME.repeats.fasta

#---------------------------Removing false positives-------------------------------------#

#	Removing false transposons by checking against SwissProt-UniProt. TIME: ~1 hour 

blastx -query $GENOME.repeats.fasta \
-db $REPO/uniprot_sprot17nov2015.fasta -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50

blastx -query  $GENOME.repeats.fasta \
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
comm -13 $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column \
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

python2 $REPO/reprint.filtered.lib.py -i $GENOME.repeats.fasta -l \
$GENOME.repeats.only_uniprot | fold -w 60 \
	> $GENOME.repeats.srt.nr.no_uniprot

#--------------------------BLASTX for classification-------------------------------------#

blastx -query $GENOME.repeats.srt.nr.no_uniprot \
-db $REPO/RepeatPeps.lib -num_threads 10 -evalue 1e-5 -max_target_seqs 1 -outfmt \
'6 bitscore length pident qseqid sacc' | awk '{print $5," ",$4," ",$1," ",$2," ",$3}' \
	> $GENOME.ClassificationScheme

python2 $REPO/classifier.py -i \
$GENOME.repeats.srt.nr.no_uniprot -c $GENOME.ClassificationScheme \
	> $GENOME.classified_repeats.fasta

awk '{print $2}' $GENOME.classified_repeats.fasta | grep "#" > list_to_fastagrep

fastagrep -v -f list_to_fastagrep $GENOME.repeats.srt.nr.no_uniprot \
	> $GENOME.non_classified_repeats.fasta

cat $GENOME.non_classified_repeats.fasta $GENOME.classified_repeats.fasta \
	> $GENOME.denovolib

#-----------------------------HMMR classification of "Unknowns"--------------------------#

#	Retrieving the FASTA sequences of unknown repeats from RepeatModeler

fastagrep -p "#Unknown" $GENOME.denovolib | fastagrep -p ">rnd" \
	> $GENOME.repmod.unknowns.lib

#	Running nhmmscan to search for homology in Dfam.

gt splitfasta -numfiles 10 $GENOME.repmod.unknowns.lib

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.1 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.1 1> $GENOME.DFAM.out.1 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.2 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.2 1> $GENOME.DFAM.out.2 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.3 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.3 1> $GENOME.DFAM.out.3 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.4 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.4 1> $GENOME.DFAM.out.4 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.5 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.5 1> $GENOME.DFAM.out.5 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.6 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.6 1> $GENOME.DFAM.out.6 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.7 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.7 1> $GENOME.DFAM.out.7 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.8 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.8 1> $GENOME.DFAM.out.8 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.9 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.9 1> $GENOME.DFAM.out.9 

nhmmscan --noali --cut_tc --dfamtblout $GENOME.repmod.unknowns.dfam_results.10 \
--cpu=16 $REPO/Dfam.hmm $GENOME.repmod.unknowns.lib.10 1> $GENOME.DFAM.out.10

cat $GENOME.repmod.unknowns.dfam_results* > $GENOME.repmod.unknowns.all.dfam_results

#	Running a set of python scripts to produce annotated versions of the libraries

python2 $REPO/dfammer.py -i $GENOME.repmod.unknowns.all.dfam_results \
	> $GENOME.repmod.unknowns.dfam_tophits

python2 $REPO/dfammer2.py -i $GENOME.repmod.unknowns.dfam_tophits \
	> $GENOME.repmod.unknowns.annotated.list

awk '{print $2}' $GENOME.repmod.unknowns.annotated.list \
	> $GENOME.repmod.unknowns.annotated.list.to_fastagrep

python2 $REPO/classifier.py -i $GENOME.repmod.unknowns.lib -c \
$GENOME.repmod.unknowns.annotated.list \
	> $GENOME.repmod.unknowns.dfam_classified.lib

#	Extracting all sequences except the unknowns not classified by HMM model homology

fastagrep -v -f $GENOME.repmod.unknowns.annotated.list.to_fastagrep $GENOME.denovolib \
	> $GENOME.denovolib.no_classifiable_unknowns

#	Removing bottom empty line

sed -i '/^$/d' $GENOME.repmod.unknowns.dfam_classified.lib 

cat $GENOME.repmod.unknowns.dfam_classified.lib \
$GENOME.denovolib.no_classifiable_unknowns \
	> $GENOME.denovolib.classified

--------------------------------RepeatMasker------------------------------------------#

#	Creating a total repeat library

cat $GENOME.denovolib.classified $REPO/RepBase_for_RepeatMasker.20150807.eukaryotes \
> $GENOME.total.repeat.library

#	Making folders for each library

mkdir $GENOME.repmask.denovo
mkdir $GENOME.repmask.total
mkdir $GENOME.repmask.repbase

#	RepeatMasking and creating summaries

RepeatMasker -lib $GENOME.denovolib.classified -gff -a -s -pa 16 -dir \
$GENOME.repmask.denovo/ $GENOME

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.denovo/*.out \
	> $GENOME.repmask.denovo.summary

RepeatMasker -lib $GENOME.total.repeat.library -gff -a -s -pa 16 -dir \
$GENOME.repmask.total/ $GENOME 

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.total/*.out \
	> $GENOME.repmask.total.summary

RepeatMasker -species eukaryota -gff -a -s -pa 16 -dir \
$GENOME.repmask.repbase/ $GENOME 

/cluster/software/VERSIONS/repeatmasker-4.0.6/util/buildSummary.pl -species teleostei \
$GENOME.repmask.repbase/*.out \
	> $GENOME.repmask.repbase.summary

#--------------------------------RepeatLandscape-----------------------------------------#

#	The createRepeatLandscape.pl script requires genome size estimates to create pie 
#	charts. 

#	Getting genome size

assemblathon_stats.pl -csv $GENOME | grep "Total size of scaffolds" | awk '{print $5}' \
> genomesize

sh $REPO/landscaping.sh $GENOME

#---------------------------------Result parsing-----------------------------------------#

cd $GENOME.repmask.total/
sh $REPO/get.rm.master.sh $GENOME
cd ..

cd $GENOME.repmask.denovo/
sh $REPO/get.rm.master.sh $GENOME
cd ..

cd $GENOME.repmask.repbase
sh $REPO/get.rm.master.sh $GENOME
cd ..

################################
#	Removal of excess files  

rm list_of_CACTAS
rm $GENOME.denovolib.classified.2
rm $GENOME.csv
rm $GENOME.fai
rm $GENOME.sds $GENOME.md5 $GENOME.esq $GENOME.des $GENOME.csv $GENOME.suf $GENOME.prj
rm $GENOME.llv $GENOME.lcp $GENOME.N50* formatdb.log error.log *.masked
rm $GENOME.TRIM* $GENOME.retro99*  $GENOME.retro85* *.dgt *.sort *Passed*
rm -r retro85 TRIM85 retro99 TRIM99 *fasta_files
rm $GENOME.retrotransposons.result99 $GENOME.retrotransposons.outinner99*
rm $GENOME.retrotransposons.out99 *retrotransposons.gff99
rm $GENOME.retrotransposons.result85 $GENOME.retrotransposons.outinner85*
rm $GENOME.retrotransposons.out85 $GENOME.retrotransposons.gff85 *withdomains*
rm $GENOME.retrotransposons.TRIMs.lib $GENOME.retrotransposons.TRIMs.lib.masked
rm *DFAM* $GENOME.repmod.unknowns.lib $GENOME.repmod.unknowns.all.dfam_results
rm $GENOME.repmod.unknowns.dfam_tophits $GENOME.repmod.unknowns.annotated.list
rm $GENOME.repmod.unknowns.annotated.list.to_fastagrep 
rm $GENOME.repmod.unknowns.dfam_classified.lib $GENOME.denovolib.no_classifiable_unknowns
rm $GENOME.tPSI.classified.fasta *80* $GENOME.repeats.fasta 
rm $GENOME.repeats.srt.nr.no_uniprot *unknowns* $GENOME.denovolib
rm list_of_CACTAS
rm $GENOME.denovolib.classified.2
rm $GENOME.ClassificationScheme $GENOME.classified_repeats.fasta list_to_fastagrep
rm $GENOME.non_classified_repeats.fasta
rm *blastx*
rm *repeats.only_uniprot
rm *repeats.only_repeatpeps
rm $GENOME.retrotransposons.TRIMs.lib.classified
rm unaligned.fa
rm halfoftable1					
rm halfoftable2	
################################

#------------------------------------End of pipeline-------------------------------------#

rsync -ravz * ${origdir}/repeats/${species} 

#!/bin/sh
#SBATCH --job-name=pipeline
#SBATCH --account=nn9244k
#SBATCH --time=24:0:0
#SBATCH --mem-per-cpu=3000M 
#SBATCH --cpus-per-task=10

if [ -n "$SLURM_JOB_ID" ]; then
    # running in a slurm job
    source /cluster/bin/jobsetup
fi
# Mail at beginning/end/on suspension
#SBATCH --mail-user=williambrynildsen@gmail.com
#SBATCH --mail-type=ALL

#------------------------------Pointing to paths-----------------------------------------#

#The repository contains databases and scripts used during the pipeline.
export REPO=/work/users/willibr/testing/repository

#---Content of the repository:

#---Databases
#	eukaryotic-tRNAs.fa
#	repbase.update.04102015.teleostei.lib
#	uniprot_sprot17nov2015.fasta
#	RepeatPeps.lib
#	TEseeker.DNAtransposons.fasta

#---Python scripts (only tested in python2/2.7.9)
#	change_headers_to_seqN.py
#	reprint.tPSI.lib.py
#	reprint.ltrharvest.lib.py
#	reprint.ltrharvest.lib.new.py
#	reprint.filtered.lib.py
#	sequence_cleaner.py
#	N50.py
#	reprint.genome.py

#---LTRharvest/LTRdigest filtering (courtesy of Ning Jiang)
#	TRIM99.CRL_Step1.pl
#	TRIM85.CRL_Step1.pl
#	retro99.CRL_Step1.pl
#	retro85.CRL_Step1.pl
#	CRL_Step2.pl
#	CRL_Step3.pl
#	CRL_Step4.pl
#	CRL_Step5.pl
#	ltr_library.pl

#---RepeatMasker parsing
#	get.rm.master.sh
#	get.rm.stats.sh
#	get.rm.stats.rename.sh
#	get.rm.overview.sh

#The only parameter given the pipeline is the genome.
export GENOME=$1

#Point to TransposonPSI
export tPSI=/projects/cees/bin/TransposonPSI/transposonPSI.pl

#Point to CD-HIT
export CDHIT=/projects/cees/bin/cd-hit/cd-hit-v4.6.4-2015-0603/cd-hit-est

#-----------------------------Loading modules-----------------------------------------------#

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

#NB! Genome FASTA headers must not contain any special characters (only letters and/or 
#numbers).

#grep ">" $ORIGINALGENOME > halfoftable1
#python2 $REPO/reprint.genome.py -i $ORIGINALGENOME > reprinted.genome
#GENOME=reprinted.genome
#grep ">" $GENOME > halfoftable2
#paste halfoftable1 halfoftable2 > $GENOME



#Starting long jobs in separate scripts.
#sh slave_LTRharvest.slurm $GENOME &
#sh slave_tPSI.slurm $GENOME &
# sh slave_RepeatModeler.slurm $GENOME &
# wait
# 
# #------------------------------Filtration and extraction---------------------------------#
# ##Running custom scripts provided by Ning Jiang (she will publish them soon) to filter
# # the LTRdigest results
# 
# perl $REPO/retro99.CRL_Step1.pl -gff $GENOME.retrotransposons.gff99.dgt
# perl $REPO/retro85.CRL_Step1.pl -gff $GENOME.retrotransposons.gff85.dgt 
# perl $REPO/TRIM85.CRL_Step1.pl -gff $GENOME.TRIM.gffT85.dgt 
# perl $REPO/TRIM99.CRL_Step1.pl -gff $GENOME.TRIM.gffT99.dgt 
# # 
# 
# ##Continuing running scripts provided by Ning Jiang.
# 
# #First, creating folders
# mkdir $GENOME.TRIM85.fasta_files
# mkdir $GENOME.TRIM99.fasta_files
# mkdir $GENOME.retro99.fasta_files
# mkdir $GENOME.retro85.fasta_files
# 
# #Copying the last script into each folder (mandatory)
# cp $REPO/CRL_Step3.pl $GENOME.TRIM85.fasta_files
# cp $REPO/CRL_Step3.pl $GENOME.TRIM99.fasta_files
# cp $REPO/CRL_Step3.pl $GENOME.retro99.fasta_files
# cp $REPO/CRL_Step3.pl $GENOME.retro85.fasta_files
# 
# #Copying other files as well, including the genome file (not good for space)
# 
# cp retro85.CRL_Step1_Passed_Elements.txt $GENOME.retro85.fasta_files/
# cp $GENOME.retrotransposons.out85 $GENOME.retro85.fasta_files/
# cp $GENOME.retrotransposons.result85 $GENOME.retro85.fasta_files/
# cp $GENOME $GENOME.retro85.fasta_files/
# 
# cp retro99.CRL_Step1_Passed_Elements.txt $GENOME.retro99.fasta_files/
# cp $GENOME.retrotransposons.out99 $GENOME.retro99.fasta_files/
# cp $GENOME.retrotransposons.result99 $GENOME.retro99.fasta_files/
# cp $GENOME $GENOME.retro99.fasta_files/
# 
# cp TRIM99.CRL_Step1_Passed_Elements.txt $GENOME.TRIM99.fasta_files/
# cp $GENOME.TRIM.outT99 $GENOME.TRIM99.fasta_files/
# cp $GENOME.TRIM.resultT99 $GENOME.TRIM99.fasta_files/
# cp $GENOME $GENOME.TRIM99.fasta_files/
# 
# cp TRIM85.CRL_Step1_Passed_Elements.txt $GENOME.TRIM85.fasta_files/
# cp $GENOME.TRIM.outT85 $GENOME.TRIM85.fasta_files/
# cp $GENOME.TRIM.resultT85 $GENOME.TRIM85.fasta_files/
# cp $GENOME $GENOME.TRIM85.fasta_files/
# 
# ##Running part II of the custom script, and moving output FASTA files into the right
# # folder
# 
# perl $REPO/CRL_Step2.pl --step1 retro85.CRL_Step1_Passed_Elements.txt --repeatfile \
# $GENOME.retrotransposons.out85 --resultfile $GENOME.retrotransposons.result85 \
# --sequencefile $GENOME --removed_repeats $GENOME.retro85.CRL_Step2_Passed_Elements.fasta
# 
# mv Repeat_* $GENOME.retro85.fasta_files/
# 
# perl $REPO/CRL_Step2.pl --step1 retro99.CRL_Step1_Passed_Elements.txt --repeatfile \
# $GENOME.retrotransposons.out99 --resultfile $GENOME.retrotransposons.result99 \
# --sequencefile $GENOME --removed_repeats $GENOME.retro99.CRL_Step2_Passed_Elements.fasta
#  
# mv Repeat_* $GENOME.retro99.fasta_files/
#  
# perl $REPO/CRL_Step2.pl --step1 TRIM99.CRL_Step1_Passed_Elements.txt --repeatfile \
# $GENOME.TRIM.outT99 --resultfile $GENOME.TRIM.resultT99  --sequencefile $GENOME \
# --removed_repeats $GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta
#  
# mv Repeat_* $GENOME.TRIM99.fasta_files/
# 
# perl $REPO/CRL_Step2.pl --step1 TRIM85.CRL_Step1_Passed_Elements.txt --repeatfile \
# $GENOME.TRIM.outT85 --resultfile $GENOME.TRIM.resultT85  --sequencefile $GENOME \
# --removed_repeats $GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta
#  
# mv Repeat_* $GENOME.TRIM85.fasta_files/
# 
# #Copying files before running part III
# 
# cp $GENOME.retro85.CRL_Step2_Passed_Elements.fasta $GENOME.retro85.fasta_files/
# cp $GENOME.retro99.CRL_Step2_Passed_Elements.fasta $GENOME.retro99.fasta_files/
# cp $GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta $GENOME.TRIM85.fasta_files/
# cp $GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta $GENOME.TRIM99.fasta_files/
# 
# cd $GENOME.retro85.fasta_files/
# 
# perl CRL_Step3.pl --directory . --step2 \
# $GENOME.retro85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 
# 
# cd ../$GENOME.TRIM85.fasta_files/
# 
# perl CRL_Step3.pl --directory . --step2 \
# $GENOME.TRIM85.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 
# 
# cd ../$GENOME.TRIM99.fasta_files/
# 
# perl CRL_Step3.pl --directory . --step2 \
# $GENOME.TRIM99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25 
# 
# cd ../$GENOME.retro99.fasta_files/
# 
# perl CRL_Step3.pl --directory . --step2 \
# $GENOME.retro99.CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
# 
# #Preparing for part IV
# 
# cd ..
# 
# mkdir retro85
# mkdir retro99
# mkdir TRIM99
# mkdir TRIM85
# 
# cd retro85
# 
# perl $REPO/ltr_library.pl --resultfile ../$GENOME.retrotransposons.result85 \
# --step3 ../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --sequencefile ../$GENOME
# 
# cp lLTR_Only.lib $GENOME.retro85.lLTR_Only.lib
# 
# cd ../retro99
# 
# perl $REPO/ltr_library.pl --resultfile ../$GENOME.retrotransposons.result99 \
# --step3 ../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --sequencefile ../$GENOME
# 
# cp lLTR_Only.lib $GENOME.retro99.lLTR_Only.lib
# 
# cd ../TRIM85
# 
# perl $REPO/ltr_library.pl --resultfile ../$GENOME.TRIM.resultT85 \
# --step3 ../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --sequencefile ../$GENOME
# 
# cp lLTR_Only.lib $GENOME.TRIM85.lLTR_Only.lib
# 
# cd ../TRIM99
# 
# perl $REPO/ltr_library.pl --resultfile ../$GENOME.TRIM.resultT99 \
# --step3 ../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --sequencefile ../$GENOME
# 
# cp lLTR_Only.lib $GENOME.TRIM99.lLTR_Only.lib
# 
# cd ..
# 
# cat retro85/$GENOME.retro85.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
# 	> retro85/$GENOME.repeats_to_mask_LTR85.fasta
# 
# cat retro99/$GENOME.retro99.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
# 	> retro99/$GENOME.repeats_to_mask_LTR99.fasta
# 
# cat TRIM85/$GENOME.TRIM85.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
# 	> TRIM85/$GENOME.repeats_to_mask_TRIM85.fasta
# 	
# cat TRIM99/$GENOME.TRIM99.lLTR_Only.lib $REPO/TEseeker.DNAtransposons.fasta \
# 	> TRIM99/$GENOME.repeats_to_mask_TRIM99.fasta
# 	
# #Repeatmasking a bit
# 
# RepeatMasker -lib retro85/$GENOME.repeats_to_mask_LTR85.fasta -pa 16 -nolow -no_is \
# -dir . $GENOME.retrotransposons.outinner85 
# 
# RepeatMasker -lib retro99/$GENOME.repeats_to_mask_LTR99.fasta -pa 16 -nolow -no_is \
# -dir . $GENOME.retrotransposons.outinner99
# 
# RepeatMasker -lib TRIM85/$GENOME.repeats_to_mask_TRIM85.fasta -pa 16 -nolow -no_is \
# -dir . $GENOME.TRIM.outinnerT85
# 
# RepeatMasker -lib TRIM99/$GENOME.repeats_to_mask_TRIM99.fasta -pa 16 -nolow -no_is \
# -dir . $GENOME.TRIM.outinnerT99
# 
# #Step 4
# MASKED1=$GENOME.retrotransposons.outinner85.masked
# MASKED2=$GENOME.retrotransposons.outinner99.masked
# MASKED3=$GENOME.TRIM.outinnerT85.masked
# MASKED3=$GENOME.TRIM.outinnerT99.masked
# 
# cd retro85
# 
# if [ -f "$MASKED1" ];
# then
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.retrotransposons.result85 --innerfile \
# 	../$GENOME.retrotransposons.outinner85.masked --sequencefile ../$GENOME
# else
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.retrotransposons.result85 --innerfile \
# 	../$GENOME.retrotransposons.outinner85 --sequencefile ../$GENOME
# fi
# 
# cd ../retro99
# 
# if [ -f "$MASKED2" ];
# then
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.retrotransposons.result99 --innerfile \
# 	../$GENOME.retrotransposons.outinner99.masked --sequencefile ../$GENOME
# else
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.retrotransposons.result99 --innerfile \
# 	../$GENOME.retrotransposons.outinner99 --sequencefile ../$GENOME
# fi
# 
# cd ../TRIM99
# 
# if [ -f "$MASKED3" ];
# then
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.TRIM.resultT99 --innerfile \
# 	../$GENOME.TRIM.outinnerT99.masked --sequencefile ../$GENOME
# else
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.TRIM.resultT99 --innerfile \
# 	../$GENOME.TRIM.outinnerT99 --sequencefile ../$GENOME
# fi
# 
# cd ../TRIM85
# 
# if [ -f "$MASKED4" ];
# then
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.TRIM.resultT85 --innerfile \
# 	../$GENOME.TRIM.outinnerT85.masked --sequencefile ../$GENOME
# else
# 	perl $REPO/CRL_Step4.pl --step3 \
# 	../$GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# 	--resultfile ../$GENOME.TRIM.resultT85 --innerfile \
# 	../$GENOME.TRIM.outinnerT85 --sequencefile ../$GENOME
# fi
# 
# cd ..
# 
# module load blast+/2.2.29
# 
# makeblastdb -in retro85/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in retro85/Inner_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in retro99/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in retro99/Inner_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in TRIM85/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in TRIM85/Inner_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in TRIM99/lLTRs_Seq_For_BLAST.fasta -dbtype nucl
# makeblastdb -in TRIM99/Inner_Seq_For_BLAST.fasta -dbtype nucl
# 
# blastn -query retro85/lLTRs_Seq_For_BLAST.fasta -db retro85/lLTRs_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out retro85/lLTRs_Seq_For_BLAST.fasta.out
# 
# blastn -query retro85/Inner_Seq_For_BLAST.fasta -db retro85/Inner_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out retro85/Inner_Seq_For_BLAST.fasta.out
# 
# blastn -query retro99/lLTRs_Seq_For_BLAST.fasta -db retro99/lLTRs_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out retro99/lLTRs_Seq_For_BLAST.fasta.out
# 
# blastn -query retro99/Inner_Seq_For_BLAST.fasta -db retro99/Inner_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out retro99/Inner_Seq_For_BLAST.fasta.out
# 
# blastn -query TRIM85/lLTRs_Seq_For_BLAST.fasta -db TRIM85/lLTRs_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out TRIM85/lLTRs_Seq_For_BLAST.fasta.out
# 
# blastn -query TRIM85/Inner_Seq_For_BLAST.fasta -db TRIM85/Inner_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out TRIM85/Inner_Seq_For_BLAST.fasta.out
# 
# blastn -query TRIM99/lLTRs_Seq_For_BLAST.fasta -db TRIM99/lLTRs_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out TRIM99/lLTRs_Seq_For_BLAST.fasta.out
# 
# blastn -query TRIM99/Inner_Seq_For_BLAST.fasta -db TRIM99/Inner_Seq_For_BLAST.fasta \
# -evalue 1e-10 -num_descriptions 1000 -out TRIM99/Inner_Seq_For_BLAST.fasta.out
# 
# 
# #STEP 5
# 
# perl $REPO/CRL_Step5.pl --LTR_blast retro85/lLTRs_Seq_For_BLAST.fasta.out \
# --inner_blast retro85/Inner_Seq_For_BLAST.fasta.out \
# --step3 $GENOME.retro85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --final $GENOME.retro85.lib
# 
# perl $REPO/CRL_Step5.pl --LTR_blast retro99/lLTRs_Seq_For_BLAST.fasta.out \
# --inner_blast retro99/Inner_Seq_For_BLAST.fasta.out \
# --step3 $GENOME.retro99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --final $GENOME.retro99.lib
# 
# perl $REPO/CRL_Step5.pl --LTR_blast TRIM85/lLTRs_Seq_For_BLAST.fasta.out \
# --inner_blast TRIM85/Inner_Seq_For_BLAST.fasta.out \
# --step3 $GENOME.TRIM85.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --final $GENOME.TRIM85.lib
# 
# perl $REPO/CRL_Step5.pl --LTR_blast TRIM99/lLTRs_Seq_For_BLAST.fasta.out \
# --inner_blast TRIM99/Inner_Seq_For_BLAST.fasta.out \
# --step3 $GENOME.TRIM99.fasta_files/CRL_Step3_Passed_Elements.fasta \
# --final $GENOME.TRIM99.lib
# 
# #Renaming headers in each library:
# 
# python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.retro85.lib \
# --name retrotransposon85 > $GENOME.retro85.renamed.lib
# 
# python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.retro99.lib \
# --name retrotransposon99 > $GENOME.retro99.renamed.lib 
# 
# python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.TRIM85.lib \
# --name TRIM85 > $GENOME.TRIM85.renamed.lib 
# 
# python2 $REPO/reprint.ltrharvest.lib.new.py -i $GENOME.TRIM99.lib \
# --name TRIM99 > $GENOME.TRIM99.renamed.lib
# 
# #Merging each group
# 
# cat $GENOME.TRIM99.renamed.lib $GENOME.TRIM85.renamed.lib $GENOME.retro99.renamed.lib \
# $GENOME.retro85.renamed.lib \
# 	> $GENOME.retrotransposons.TRIMs.lib
# 
# #Trying to classify the consensus sequences using RepeatClassifier
# 
# RepeatClassifier -consensi $GENOME.retrotransposons.TRIMs.lib -engine ncbi
# 
# #Output is $GENOME.retrotransposons.TRIMs.lib.classified.
# #Elements classified as simple repeats, buffers or unknowns are reclassified into the
# #LTR superfamily using sed. 
# 
# sed 's/#Simple_repeat/#LTR\/Unknown/g' \
# $GENOME.retrotransposons.TRIMs.lib.classified | \
# sed 's/#buffer/#LTR\/Unknown/g' | \
# sed 's/#Unknown/#LTR\/Unknown/g'| \
# sed 's/#DNA.*/#LTR\/Unknown/g' | \
# sed 's/LINE.*/#LTR\/Unknown/g' \
# 	> $GENOME.retrotransposons.TRIMs.lib.classified.reclassified
# 
# 
# #-------------------------------tPSI FASTA extraction------------------------------------#
# 
# #In order to extract sequences from the TransposonPSI file, files had to be altered
# cut -f 1,4,5,9 $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | sed 's/ID=.*Target.//g' | \
# sed 's/;.*//g' | column -t | awk '{print $4}' \
# 	> tPSI.classes
# 
# paste tPSI.classes $GENOME.TPSI.allHits.chains.bestPerLocus.gff3 | column -s $'\t' -t \
# 	> tmp4
# 
# awk '{print $2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' tmp4 \
# 	> tmp5
# 
# sed 's/ /\t/g' tmp5 \
# 	> tmp6
# 
# #The actual extraction done using bedtools (getfasta)
# 
# bedtools getfasta -name -fi $GENOME -bed tmp6 -fo $GENOME.tPSI.fasta
# 
# #The TransposonPSI library is classified according to protein homology
# python2 $REPO/reprint.tPSI.lib.py -i $GENOME.tPSI.fasta | fold -w 60 \
# 	> $GENOME.tPSI.classified.fasta
# 
# #--------------------------------Merging of libraries------------------------------------#
# 
cd RM*
sed 's/ .*//g' consensi.fa.classified > consensi.fa.classified.stripped
cd ..
cat $GENOME.retrotransposons.TRIMs.lib.classified.reclassified \
RM*/consensi.fa.classified.stripped \
$GENOME.tPSI.classified.fasta | fold -w 60 > $GENOME.repeats.fasta

#---------------------------------Clustering---------------------------------------------#
#Clustering of remaining elements are done using CD-HIT-EST with standard settings,
#meaning that elements 90% similar will be clustered together. This clustering 
#algorithm is not perfect for long sequences. 80% is preferred to 90%. Takes ~10 min.

$CDHIT -i $GENOME.repeats.fasta -o $GENOME.repeats.srt.nr.80 -c 0.8 -T 10

#---------------------------Removing false positives-------------------------------------#

#Removing false transposons by checking against SwissProt-UniProt. TIME: ~1 hour 

blastx -query $GENOME.repeats.srt.nr.80 \
-db $REPO/uniprot_sprot17nov2015.fasta -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out $GENOME.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50

blastx -query  $GENOME.repeats.srt.nr.80 \
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

python2 $REPO/reprint.filtered.lib.py -i $GENOME.repeats.srt.nr -l \
$GENOME.repeats.only_uniprot | fold -w 60 \
	> $GENOME.repeats.srt.nr.no_uniprot
# 
#------------------------------Removing singletons---------------------------------------#

# makeblastdb -in $GENOME -dbtype nucl
# blastn -query $GENOME.repeats.srt.nr.no_uniprot \
# -db $GENOME -num_threads 1 -evalue 1e-5 -max_target_seqs 10 \
# -outfmt 6 -out denovolib.vs.$GENOME.blastn.1e-5.max50
# 
# cat denovolib.vs.$GENOME.blastn.1e-5.max50 | cut -f 1 | uniq -c | sort -n | column -t | \
# grep "1     " | cut -f 2 | awk '{ print $2}' > $GENOME.list_of_singletons
# 
# python2 $REPO/reprint.filtered.lib.py -i $GENOME.repeats.srt.nr.no_uniprot \
# -l $GENOME.list_of_singletons > $GENOME.repeats.srt.nr.no_uniprot.no_singletons

#-----------------------Removing simple repeats designated as TEs------------------------#

#First splitting the library in elements found by tPSI and elements detected by other
#programs.

# fastagrep -i -p "TransposonPSI" $GENOME.repeats.srt.nr.no_uniprot \
# 	> $GENOME.denovolib_tPSI
# 	
# fastagrep -i -v -p "TransposonPSI" $GENOME.repeats.srt.nr.no_uniprot \
# 	> $GENOME.denovolib_non_tPSI
# 
# 
# #The transposons designated as CACTA by tPSI are mostly simple repeats. These are 
# #masked using TRF and then all sequences that were masked are filtered out. 
# 
# fastagrep -i -p "/CACTA" $GENOME.denovolib_tPSI \
# 	> $GENOME.denovolib_tPSI.CACTA_elements_only
# 
# fastagrep -i -v -p "/CACTA" $GENOME.denovolib_tPSI \
# 	> $GENOME.denovolib_tPSI.no_CACTA
# 
# trf $GENOME.denovolib_tPSI.CACTA_elements_only \
# 2 7 7 80 10 50 500 -h -m
# 
# #output: $GENOME.denovolib_tPSI.CACTA_elements_only.2.7.7.80.10.50.500.mask
# 
# #Only keeping sequences that was not masked at all (cutoff = 1%)
# 
# python2 $REPO/sequence_cleaner.py \
# $GENOME.denovolib_tPSI.CACTA_elements_only.2.7.7.80.10.50.500.mask 0 1
# 
# #output: $clear_GENOME.denovolib_tPSI.CACTA_elements_only.2.7.7.80.10.50.500.mask 
# 
# #Merging library:
# 
# cat $GENOME.denovolib_tPSI.no_CACTA \
# clear_denovolib_tPSI.CACTA_elements_only.2.7.7.80.10.50.500.mask \
# $GENOME.denovolib_non_tPSI \
# 	> $GENOME.final_denovo_lib
# 
# #A final filter of 80 bp as not many TEs are smaller than that.
# 
# python2 $REPO/sequence_cleaner.py $GENOME.final_denovo_lib 80 0 

#output clear_$GENOME.final_denovo_lib


# #--------------------------------RepeatMasker------------------------------------------#
# 
##Creating a total repeat library

# cat clear_$GENOME.final_denovo_lib $REPO/repbase.update.04102015.teleostei.lib \
# > $GENOME.total.repeat.library
# 
# #Making folders for each library
# 
# mkdir $GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output
# mkdir $GENOME.total.repeat.library.repeatmasker_output
# mkdir repbase.update.04102015.teleostei.lib.repeatmasker_output
# 
# RepeatMasker -lib clear_$GENOME.final_denovo_lib -gff -a -s -dir \
# clear_$GENOME.final_denovo_lib.repeatmasker_output/ $GENOME
# 
# RepeatMasker -lib $GENOME.total.repeat.library -gff -a -s -dir \
# $GENOME.total.repeat.library.repeatmasker_output/ $GENOME
# 
# RepeatMasker -lib $REPO/repbase.update.04102015.teleostei.lib -gff -a -s -dir \
# repbase.update.04102015.teleostei.lib.repeatmasker_output/ $GENOME
# 
# #Creating summaries
# 
# /cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
# $GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output/*.out \
# 	> $GENOME.repmask.denovo.summary
# 	
# /cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
# $GENOME.total.repeat.library.repeatmasker_output/*.out \
# 	> $GENOME.repmask.total.summary
# 	
# /cluster/software/repeatmasker/util/buildSummary.pl -species teleostei \
# repbase.update.04102015.teleostei.lib.repeatmasker_output/*.out\
# 	> $GENOME.repmask.repbase.summary

#---------------------------------Result parsing-----------------------------------------#

# cd $GENOME.repeats.srt.nr.no_uniprot.repeatmasker_output
# get.rm.master.sh $GENOME
# cd ..
# 
# cd $GENOME.total.repeat.library.repeatmasker_output
# get.rm.master.sh $GENOME
# cd ..
# 
# cd repbase.update.04102015.teleostei.lib.repeatmasker_output
# get.rm.master.sh $GENOME
# cd ..


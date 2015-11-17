#Repeats -- Work in progress

In short the pipeline does the following:

* LTRharvest are run for collection of LTR retrotransposons and TRIMs
* LTRdigest refines the output of LTRharvest by looking for primer binding sites and retrotransposon specific enzymes. FASTA sequences are extracted from the GFF file produced.
* Repeats are also collected using RepeatModeler.
* More divergent LTR retrotransposons and DNA transposons are collected using TransposonPSI. FASTA sequences are extracted from the GFF and classified according to protein homology.
* All sequences are clustered using CD-HIT-EST
* The sequences are subject to BLASTX searches against SwissProt-Uniprot and RepeatMaskers repeat peptide library. If hits occur in the uniprot database, but not the repeat peptide database the sequence is discarded.
* RepeatMasker is run three times; only with the _de novo_ library, only with the RepBase library and once with a merged variant.
* Finally, summaries of the RepeatMasker results are produced.

Flowchart:
![alt text](https://github.com/uio-cees/willibr-TE-pipeline/blob/master/repeatspipeline.png)


The pipeline is dependent on these programs:

[genometools/1.5.7](http://www.genometools.org)

[repeatmasker/4.0.5](http://www.repeatmasker.org/RepeatModeler.html)

[blast/2.2.26](http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/release/2.2.26/)

[hmmer/3.0](http://hmmer.janelia.org/download.html)

[muscle/3.8.31](http://www.drive5.com/muscle/downloads.htm)

[python2/2.7.9](https://www.python.org/downloads/release/python-279/)

[usearch/7.0.1090](http://www.drive5.com/usearch/download.html)

[cd-hit/4.6.4](http://weizhongli-lab.org/cd-hit/)

[blast+/2.2.29](http://www.ncbi.nlm.nih.gov/books/NBK279671/)

[repeatmasker/4.0.5](http://www.repeatmasker.org/RMDownload.html)

[bedtools/2.17.0](https://code.google.com/p/bedtools/downloads/list)


--

You will need these files to be in the directory along with the genome of study.

>eukaryotic-tRNAs.fa
>
>retro99.custom_script1.pl
>
>retro85.custom_script1.pl
>
>TRIM85.custom_script1.pl
>
>TRIM99.custom_script1.pl
>
>custom_script2.pl
>
>filter\_protein_match.lua
>
>custom_script3.pl
>
>change\_headers\_to_seqN.py
>
>reprint.tPSI.lib.py
>
>reprint.ltrharvest.lib.py
>
>reprint.filtered.lib.py

>uniprot_sprot.fasta

>RepeatPeps.lib

>repbase.update.lib

####Load necessary modules

```perl
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
```

##Indexing genome

####For LTRharvest/LTRdigest
```perl
gt suffixerator -db scaffolds -indexname scaffolds -tis -suf -lcp -des -ssp -sds -dna
```

####For RepeatModeler
```perl
BuildDatabase -name scaffolds  -engine ncbi scaffolds 
```

##RepeatModeler
```perl
RepeatModeler -database scaffolds -engine ncbi -pa 10
```

##TransposonPSI
```perl
perl /projects/cees/bin/TransposonPSI/transposonPSI.pl scaffolds nuc >tPSI.log
```


##LTRharvest
```perl
gt ltrharvest -index scaffolds -out scaffolds.retrotransposons.out99 \
-outinner scaffolds.retrotransposons.outinner99 -gff3 scaffolds.retrotransposons.gff99 -minlenltr 100 -maxlenltr 6000 \
-mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 > scaffolds.retrotransposons.result99

gt ltrharvest -index scaffolds -out scaffolds.retrotransposons.out85 \
-outinner scaffolds.retrotransposons.outinner85 -gff3 scaffolds.retrotransposons.gff85 -minlenltr 100 -maxlenltr 6000 \
-mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10  > scaffolds.retrotransposons.result85

gt ltrharvest -index scaffolds -out scaffolds.TRIM.outT99 -outinner \
scaffolds.TRIM.outinnerT99 -gff3 scaffolds.TRIM.gffT99 -minlenltr 70 -maxlenltr 500 \
-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 > scaffolds.TRIM.resultT99

gt ltrharvest -index scaffolds -out scaffolds.TRIM.outT85 -outinner \
scaffolds.TRIM.outinnerT85 -gff3 scaffolds.TRIM.gffT85 -minlenltr 70 -maxlenltr 500 \
-mindistltr 280 -maxdistltr 1500 -mintsd 5 -maxtsd 5 -vic 10 > scaffolds.TRIM.resultT85
```

####Sorting for LTRdigest
```perl
gt gff3 -sort scaffolds.retrotransposons.gff99 > scaffolds.retrotransposons.gff99.sort
gt gff3 -sort scaffolds.retrotransposons.gff85 > scaffolds.retrotransposons.gff85.sort
gt gff3 -sort scaffolds.TRIM.gffT85 > scaffolds.TRIM.gffT85.sort
gt gff3 -sort scaffolds.TRIM.gffT99 > scaffolds.TRIM.gffT99.sort 
```
##LTRdigest
```perl
gt ltrdigest -trnas /work/users/willibr/testing/folder/eukaryotic-tRNAs.fa -hmms /work/users/willibr/testing/repository/gydb/*hmm -- scaffolds.retrotransposons.gff99.sort \
scaffolds> scaffolds.retrotransposons.gff99.dgt
```
```perl
gt ltrdigest -trnas /work/users/willibr/testing/folder/eukaryotic-tRNAs.fa -hmms /work/users/willibr/testing/repository/gydb/*hmm -- scaffolds.retrotransposons.gff85.sort \
scaffolds> scaffolds.retrotransposons.gff85.dgt 
```
```perl
gt ltrdigest -trnas /work/users/willibr/testing/folder/eukaryotic-tRNAs.fa -hmms /work/users/willibr/testing/repository/gydb/*hmm -- scaffolds.TRIM.gffT85.sort \
scaffolds> scaffolds.TRIM.gffT85.dgt 
```
```perl
gt ltrdigest -trnas /work/users/willibr/testing/folder/eukaryotic-tRNAs.fa -hmms /work/users/willibr/testing/repository/gydb/*hmm -- scaffolds.TRIM.gffT99.sort \
scaffolds> scaffolds.TRIM.gffT99.dgt 
```

####Filtering LTRdigest results
```perl
perl retro99.custom_script1.pl -gff scaffolds.retrotransposons.gff99.dgt
perl retro85.custom_script1.pl -gff scaffolds.retrotransposons.gff85.dgt 
perl TRIM85.custom_script1.pl -gff scaffolds.TRIM.gffT85.dgt 
perl TRIM99.custom_script1.pl -gff scaffolds.TRIM.gffT99.dgt 
```

####Seleting retrotransposons with enzymes
```perl
gt select -rule_files /work/users/willibr/testing/folder/filter_protein_match.lua -- < scaffolds.retrotransposons.gff99.dgt > scaffolds.retrotransposons.gff99.dgt.withdomains
gt select -rule_files /work/users/willibr/testing/folder/filter_protein_match.lua -- < scaffolds.retrotransposons.gff85.dgt > scaffolds.retrotransposons.gff85.dgt.withdomains
```
####Making folders
```perl
mkdir scaffolds.TRIM85.fasta_files
mkdir scaffolds.TRIM99.fasta_files
mkdir scaffolds.retro99.fasta_files
mkdir scaffolds.retro85.fasta_files
```
####Copying files to folders
```perl
cp /work/users/willibr/testing/custom_script3.pl scaffolds.TRIM85.fasta_files
cp /work/users/willibr/testing/custom_script3.pl scaffolds.TRIM99.fasta_files
cp /work/users/willibr/testing/custom_script3.pl scaffolds.retro99.fasta_files
cp /work/users/willibr/testing/custom_script3.pl scaffolds.retro85.fasta_files

cp retro85.custom_script1_Passed_Elements.txt scaffolds.retro85.fasta_files/
cp scaffolds.retrotransposons.out85 scaffolds.retro85.fasta_files/
cp scaffolds.retrotransposons.result85 scaffolds.retro85.fasta_files/
cp scaffolds scaffolds.retro85.fasta_files/

cp retro99.custom_script1_Passed_Elements.txt scaffolds.retro99.fasta_files/
cp scaffolds.retrotransposons.out99 scaffolds.retro99.fasta_files/
cp scaffolds.retrotransposons.result99 scaffolds.retro99.fasta_files/
cp scaffolds scaffolds.retro99.fasta_files/

cp TRIM99.custom_script1_Passed_Elements.txt scaffolds.TRIM99.fasta_files/
cp scaffolds.TRIM.outT99 scaffolds.TRIM99.fasta_files/
cp scaffolds.TRIM.resultT99 scaffolds.TRIM99.fasta_files/
cp scaffolds scaffolds.TRIM99.fasta_files/

cp TRIM85.custom_script1_Passed_Elements.txt scaffolds.TRIM85.fasta_files/
cp scaffolds.TRIM.outT85 scaffolds.TRIM85.fasta_files/
cp scaffolds.TRIM.resultT85 scaffolds.TRIM85.fasta_files/
cp scaffolds scaffolds.TRIM85.fasta_files/
```

####Extract flanks for later alignment

```perl
perl custom_script2.pl --step1 retro85.custom_script1_Passed_Elements.txt --repeatfile \
scaffolds.retrotransposons.out85 --resultfile scaffolds.retrotransposons.result85  --sequencefile scaffolds \
--removed_repeats scaffolds.retro85.custom_script2_Passed_Elements.fasta

wait

mv Repeat_* scaffolds.retro85.fasta_files/

wait

perl custom_script2.pl --step1 retro99.custom_script1_Passed_Elements.txt --repeatfile \
scaffolds.retrotransposons.out99 --resultfile scaffolds.retrotransposons.result99  --sequencefile scaffolds \
--removed_repeats scaffolds.retro99.custom_script2_Passed_Elements.fasta

wait 

mv Repeat* scaffolds.retro99.fasta_files/

wait

perl custom_script2.pl --step1 TRIM99.custom_script1_Passed_Elements.txt --repeatfile \
scaffolds.TRIM.outT99 --resultfile scaffolds.TRIM.resultT99  --sequencefile scaffolds \
--removed_repeats scaffolds.TRIM99.custom_script2_Passed_Elements.fasta

wait

mv Repeat* scaffolds.TRIM99.fasta_files/

wait

perl custom_script2.pl --step1 TRIM85.custom_script1_Passed_Elements.txt --repeatfile \
scaffolds.TRIM.outT85 --resultfile scaffolds.TRIM.resultT85  --sequencefile scaffolds \
--removed_repeats scaffolds.TRIM85.custom_script2_Passed_Elements.fasta

wait

mv Repeat* scaffolds.TRIM85.fasta_files/
```


####Aligning flanks

```perl
cd scaffolds.retro85.fasta_files/

wait

perl custom_script3.pl --directory . --step2  \
../scaffolds.retro85.custom_script2_Passed_Elements.fasta --pidentity 60 --seq_c 25

cd  ../scaffolds.TRIM85.fasta_files/

wait

perl custom_script3.pl --directory . --step2  \
../scaffolds.TRIM85.custom_script2_Passed_Elements.fasta --pidentity 60 --seq_c 25

cd ../scaffolds.TRIM99.fasta_files/

wait 

perl custom_script3.pl --directory . --step2  \
../scaffolds.TRIM99.custom_script2_Passed_Elements.fasta --pidentity 60 --seq_c 25

cd scaffolds.retro99.fasta_files/

wait

perl custom_script3.pl --directory . --step2  \
../scaffolds.retro99.custom_script2_Passed_Elements.fasta --pidentity 60 --seq_c 25
```

###Extract sequence of retrotransposons with domains

```perl
python2 change_headers_to_seqN.py -i scaffolds | sed 's/ //g' > scaffolds.changed_headers

grep "ID=repeat_region" scaffolds.retrotransposons.gff85.dgt.withdomains > scaffolds.retrotransposons.gff85.dgt.withdomains.full
grep "ID=repeat_region" scaffolds.retrotransposons.gff99.dgt.withdomains > scaffolds.retrotransposons.gff99.dgt.withdomains.full

wait

bedtools getfasta -fi scaffolds.changed_headers -bed scaffolds.retrotransposons.gff85.dgt.withdomains.full -fo scaffolds.retrotransposons.gff85.dgt.withdomains.full.fasta
bedtools getfasta -fi scaffolds.changed_headers -bed scaffolds.retrotransposons.gff99.dgt.withdomains.full -fo scaffolds.retrotransposons.gff99.dgt.withdomains.full.fasta
```

###Extract sequences from TransposonPSI results
```perl
cut -f 1,4,5,9 scaffolds.TPSI.allHits.chains.bestPerLocus.gff3 > tmp 

wait

sed 's/ID=.*Target.//g' tmp > tmp2

wait

sed 's/;.*//g' tmp2 | column -t > tmp3

awk '{print $4}' tmp3 > tPSI.classes

wait 

paste tPSI.classes scaffolds.TPSI.allHits.chains.bestPerLocus.gff3 | column -s $'\t' -t > tmp4

wait 

awk '{print $2,$3,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' tmp4 > tmp5

wait

sed 's/ /\t/g' tmp5 > tmp6

wait

bedtools getfasta -name -fi scaffolds -bed tmp6 -fo scaffolds.tPSI.fasta

```

###Classify the tPSI library in a RepeatMasker fashion

```perl
python2 reprint.tPSI.lib.py -i scaffolds.tPSI.fasta | fold -w 60 > scaffolds.tPSI.classified.fasta
```

###Merge libs, renaming headers in the process
```perl
cp scaffolds.retro99.fasta_files/custom_script3_Passed_Elements.fasta scaffolds.retro99.fasta
cp scaffolds.retro85.fasta_files/custom_script3_Passed_Elements.fasta scaffolds.retro85.fasta
cp scaffolds.TRIM99.fasta_files/custom_script3_Passed_Elements.fasta scaffolds.TRIM99.fasta
cp scaffolds.TRIM85.fasta_files/custom_script3_Passed_Elements.fasta scaffolds.TRIM85.fasta

wait 

python2 reprint.ltrharvest.lib.py -i scaffolds.retro99.fasta --name retro99 --type Unknown > scaffolds.retro99.renamed.fasta
python2 reprint.ltrharvest.lib.py -i scaffolds.retro85.fasta --name retro85 --type Unknown > scaffolds.retro85.renamed.fasta
python2 reprint.ltrharvest.lib.py -i scaffolds.TRIM99.fasta --name TRIM99 --type TRIM > scaffolds.TRIM99.renamed.fasta
python2 reprint.ltrharvest.lib.py -i scaffolds.TRIM85.fasta --name TRIM85 --type TRIM > scaffolds.TRIM85.renamed.fasta
python2 reprint.ltrharvest.lib.py -i scaffolds.retrotransposons.gff99.dgt.withdomains.full.fasta --name retro99withdomains --type Unknown > scaffolds.retrotransposons.gff99.dgt.withdomains.full.renamed.fasta
python2 reprint.ltrharvest.lib.py -i scaffolds.retrotransposons.gff85.dgt.withdomains.full.fasta --name retro85withdomains --type Unknown > scaffolds.retrotransposons.gff85.dgt.withdomains.full.renamed.fasta

cat scaffolds.retro99.renamed.fasta \
scaffolds.retro85.renamed.fasta \
scaffolds.TRIM99.renamed.fasta \
scaffolds.TRIM85.renamed.fasta \
scaffolds.retrotransposons.gff99.dgt.withdomains.full.renamed.fasta \
scaffolds.retrotransposons.gff85.dgt.withdomains.full.renamed.fasta \
> scaffolds.retrotransposons.TRIMs.fasta

wait 

cat scaffolds.retrotransposons.TRIMs.fasta \
RM*/consensi.fa.classified \
scaffolds.tPSI.classified.fasta \
> scaffolds.repeats.fasta
```

###Cluster the sequences
```perl
usearch -sortbylength scaffolds.repeats.fasta \
--output scaffolds.repeats.srt
```

NB:Do this instead of using usearch:

```
cd-hit-est -i repeats.srt -o scaffolds.repeats.srt.nr -c 0.80 -n 5 -T 10
```

wait

usearch -cluster_fast scaffolds.repeats.srt \
--id 0.8 -centroids scaffolds.repeats.srt.consensus \
--uc scaffolds.repeats.srt.clusters.uc --consout scaffolds.repeats.srt.nr.bad_headers

wait 

sed 's/centroid=//g' scaffolds.repeats.srt.nr.bad_headers > scaffolds.repeats.srt.nr.bad_headers2
sed 's/;.*//g' scaffolds.repeats.srt.nr.bad_headers2 > scaffolds.repeats.srt.nr

```

####Remove false positives (non transposon genes)
```perl
makeblastdb -in scaffolds.repeats.srt.nr -dbtype nucl 

wait 

blastx -query scaffolds.repeats.srt.nr \
-db uniprot_sprot.fasta -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50

blastx -query  scaffolds.repeats.srt.nr \
-db RepeatPeps.lib -num_threads 10 -evalue 1e-5 -max_target_seqs 50 \
-outfmt 6 -out scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50

sort -k1,1 -k12,12nr -k11,11n \
scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring

sort -k1,1 -k12,12nr -k11,11n \
scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50 | sort -u -k1,1 \
--merge > scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring

cut -f 1 scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring \
> scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column

cut -f 1 scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.highest_scoring \
> scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column

echo "Num entries only hit to RepeatPeps" 
comm -13  scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column > scaffolds.repeats.only_repeatpeps
cat scaffolds.repeats.only_repeatpeps |wc -l
echo "Num entries only hit to UniProt"
comm -23 scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column > scaffolds.repeats.only_uniprot
cat scaffolds.repeats.only_uniprot |wc -l
echo "Num entries hit to both"
comm -12 scaffolds.repeats.srt.nr.uniprot-sprot.blastx.1e-5.max50.first_column scaffolds.repeats.srt.nr.repeatpeps.blastx.1e-5.max50.first_column |wc -l

wait

sed 's/\[/(/g' scaffolds.repeatlib.only_uniprot > scaffolds.repeatlib.only_uniprot.no_brackets1
sed 's/]/)/g' scaffolds.repeatlib.only_uniprot.no_brackets1 > scaffolds.repeatlib.only_uniprot.no_brackets2
sed 's/\[/(/g' scaffolds.repeatlib.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring > scaffolds.repeatlib.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring.no_brackets1
sed 's/]/)/g' scaffolds.repeatlib.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring.no_brackets1 > scaffolds.repeatlib.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring.no_brackets2

wait 

grep -f  scaffolds.repeatlib.only_uniprot.no_brackets2 scaffolds.repeatlib.srt.nr.uniprot-sprot.blastx.1e-5.max50.highest_scoring.no_brackets2 > scaffolds.repeatlib.only_uniprot_blast_hit.no_brackets

wait

sed  's/(/\[/g' scaffolds.repeatlib.only_uniprot_blast_hit.no_brackets > scaffolds.repeatlib.only_uniprot_blast_hit1
sed  's/)/]/g' scaffolds.repeatlib.only_uniprot_blast_hit1 > scaffolds.repeatlib.only_uniprot_blast_hit

wait

python2 reprint.filtered.lib.py -i scaffolds.repeats.srt.nr -l scaffolds.repeats.only_uniprot | fold -w 60 > scaffolds.repeats.srt.nr.no_uniprot

```

##RepeatMasker

```perl
sed 's/ .*//g' scaffolds.repeats.srt.nr.no_uniprot > scaffolds.repeats.srt.nr.no_uniprot.stripped

wait 

cat scaffolds.repeats.srt.nr.no_uniprot.stripped repbase.update.lib > scaffolds.total.repeat.library

wait 
RepeatMasker -lib scaffolds.repeats.srt.nr.no_uniprot.stripped -a -s -pa 10 -dir scaffolds.repmask.denovo/ scaffolds
RepeatMasker -lib scaffolds.total.repeat.library -a -s -pa 10 -dir scaffolds.repmask.total/ scaffolds
RepeatMasker -lib repbase.update.lib -a -s -pa 10 -dir scaffolds.repmask.repbase/ scaffolds
```
#Create summaries

```perl
/cluster/software/repeatmasker/util/buildSummary.pl -species eukaryota scaffolds.repmask.denovo/scaffolds.out > scaffolds.repmask.denovo.summary
/cluster/software/repeatmasker/util/buildSummary.pl -species eukaryota scaffolds.repmask.total/scaffolds.out > scaffolds.repmask.total.summary
/cluster/software/repeatmasker/util/buildSummary.pl -species eukaryota scaffolds.repmask.repbase/scaffolds.out > scaffolds.repmask.repbase.summary
```

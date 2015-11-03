# Repeat library creation

##1. Masking simple repeats

```perl 
module load repeatmasker/4.0.5
RepeatMasker -noint -pa 16 -species eukaryota -dir . genome.fa
```
Important output files: genome.fa.masked and genome.fa.out 
							
genome.fa.masked is a FASTA file of the genome with simple repeats masked as N. This also produces a GFF file of detected simple repeats that can be used as a simple repeat track in a genome browser.

The reason for masking simple repeats first is that transposable element detection software often gets confused by simple repeats and designate simple repeat regions as transposable elements. 



##2. LTRharvest and LTRdigest

LTRharvest searches the genome for putative LTR retrotransposons. LTRharvest needs index files in order for it to work. The suffixerator command in GenomeTools produces such index files.

```perl
module load genometools/1.5.7
gt suffixerator -db genome.fa.masked -indexname genome -tis -suf -lcp -des -ssp -sds -dna
```
Time: 10 minutes for a 600 mb FASTA file.

###Running LTRharvest
```perl
module load genometools/1.5.7
gt ltrharvest -index genome -gff3 genome.fa.gff -out genome.fa.ltr.fa
```
This took 3 minutes. Number of putative elements found is checked by running:

```perl
grep "ID=repeat_region" genome.fa.gff | wc l
```

###Running LTRdigest

Download HMM (hidden markov models) of known retrotransposon enzymes from: __http://gydb.org/index.php/Collection_HMM__

Sort the genome.fa.gff file:

```
gt gff3 -sort genome.fa.gff > genome.fa.sorted.gff
```

Then run LTRdigest using -hmms option

```perl
module load genometools/1.5.7
module load hmmer/3.0

gt -j 10 ltrdigest -hmms ../gydb-hmmr-complete/*hmm -outfileprefix gadMor2_ltrdigest ../gadMor2.masked.ltrharvest.sort.gff3 ../suffix_gadMor2_masked/gadMor2masked > gadMor2_ltrdigest.output.gff

```
gadMor2_ltrdigest.output.gff is a GFF file containing all LTR retrotransposons with some containing LTR retrotransposon proteins. In able to only get the elements containing proteins we need to use the 'select' command using GenomeTools.

Sort the gadMor2_ltrdigest.output.gff file:

```perl
gt gff3 -sort gadMor2_ltrdigest.output.gff > gadMor2_ltrdigest.output.sort.gff
```
And select for the relevant LTR retrotransposons.

```
gt select -rule_files filter_protein_match.lua -- < gadMor2_ltrdigest.output.sort.gff > gadMor2_ltrdigest.output.withdomains.gff
```

The rule file (from **https://github.com/satta/ltrsift/tree/master/filters**)

```
name        = "Protein Domain Filter"
author      = "Sascha Kastens"
version     = "1.0"
email       = "mail@skastens.de"
short_descr = "Filters out candidates without protein domains"
description = "Filters out a candidate if it does not contain at " ..
              "least one node of type 'protein_match'."

function filter(gn)
  gfi = gt.feature_node_iterator_new(gn)
  node = gfi:next()
  while not (node == nil) do
    if (node:get_type() == "protein_match") then
      return false
    end
    node = gfi:next()
  end
  return true
end
```

I did not manage to retrieve fasta files from gadMor2_ltrdigest.output.withdomains.gff using the entire gadMor2 assembly due to LTRharvests way of indexing the genome. I made a version of the gadMor2 genome only containing the assembled chromosomes (LG01-LG23). Due to LTRharvest indexing and getting errors running the 'extractfeat' command in GenomeTools I instead ran 'gt select' to get all LTR retrotransposon found in the assembled chromosomes.

```perl
gt select -seqid seq0 gadMor2_ltrdigest.output.withdomains.gff > LG01.ltrs.withdomains.gff
gt select -seqid seq1 gadMor2_ltrdigest.output.withdomains.gff > LG02.ltrs.withdomains.gff
gt select -seqid seq2 gadMor2_ltrdigest.output.withdomains.gff > LG03.ltrs.withdomains.gff
gt select -seqid seq3 gadMor2_ltrdigest.output.withdomains.gff > LG04.ltrs.withdomains.gff
gt select -seqid seq4 gadMor2_ltrdigest.output.withdomains.gff > LG05.ltrs.withdomains.gff
gt select -seqid seq5 gadMor2_ltrdigest.output.withdomains.gff > LG06.ltrs.withdomains.gff
gt select -seqid seq6 gadMor2_ltrdigest.output.withdomains.gff > LG07.ltrs.withdomains.gff
gt select -seqid seq7 gadMor2_ltrdigest.output.withdomains.gff > LG08.ltrs.withdomains.gff
gt select -seqid seq8 gadMor2_ltrdigest.output.withdomains.gff > LG09.ltrs.withdomains.gff
gt select -seqid seq9 gadMor2_ltrdigest.output.withdomains.gff > LG10.ltrs.withdomains.gff
gt select -seqid seq10 gadMor2_ltrdigest.output.withdomains.gff > LG11.ltrs.withdomains.gff
gt select -seqid seq11 gadMor2_ltrdigest.output.withdomains.gff > LG12.ltrs.withdomains.gff
gt select -seqid seq12 gadMor2_ltrdigest.output.withdomains.gff > LG13.ltrs.withdomains.gff
gt select -seqid seq13 gadMor2_ltrdigest.output.withdomains.gff > LG14.ltrs.withdomains.gff
gt select -seqid seq14 gadMor2_ltrdigest.output.withdomains.gff > LG15.ltrs.withdomains.gff
gt select -seqid seq15 gadMor2_ltrdigest.output.withdomains.gff > LG16.ltrs.withdomains.gff
gt select -seqid seq16 gadMor2_ltrdigest.output.withdomains.gff > LG17.ltrs.withdomains.gff
gt select -seqid seq17 gadMor2_ltrdigest.output.withdomains.gff > LG18.ltrs.withdomains.gff
gt select -seqid seq18 gadMor2_ltrdigest.output.withdomains.gff > LG19.ltrs.withdomains.gff
gt select -seqid seq19 gadMor2_ltrdigest.output.withdomains.gff > LG20.ltrs.withdomains.gff
gt select -seqid seq20 gadMor2_ltrdigest.output.withdomains.gff > LG21.ltrs.withdomains.gff
gt select -seqid seq21 gadMor2_ltrdigest.output.withdomains.gff > LG22.ltrs.withdomains.gff
gt select -seqid seq22 gadMor2_ltrdigest.output.withdomains.gff > LG23.ltrs.withdomains.gff
```
Then all the individual retrotransposons were selected from each GFF and a new GFF was made:

```perl 
grep "ID=repeat_region" LG* >LGALL.ltrs.withdomains
```
A little bit of cleaning of the output:
```perl
sed 's/.*://g' LGALL.ltrs.withdomains > LGALL.ltrs.withdomains.seq
```
Since LTRdigest is calling LG01 seq0 and so on I had to rename all the 'seq' into there corresponding LG:

```perl
sed 's/seq22/LG23/g' LGALL.ltrs.withdomains.seq> tmp2
sed 's/seq21/LG22/g' tmp2 > tmp
sed 's/seq20/LG21/g' tmp > tmp2
sed 's/seq19/LG20/g' tmp2 > tmp
sed 's/seq17/LG18/g' tmp > tmp2
sed 's/seq16/LG17/g' tmp2 > tmp
sed 's/seq15/LG16/g' tmp > tmp2
sed 's/seq14/LG15/g' tmp2 > tmp
sed 's/seq13/LG14/g' tmp > tmp2
sed 's/seq12/LG13/g' tmp2 > tmp
sed 's/seq11/LG12/g' tmp > tmp2
sed 's/seq10/LG11/g' tmp2 > tmp
sed 's/seq9/LG10/g' tmp > tmp2
sed 's/seq8/LG09/g' tmp2 > tmp
sed 's/seq7/LG08/g' tmp > tmp2
sed 's/seq6/LG07/g' tmp2 > tmp
sed 's/seq5/LG06/g' tmp > tmp2
sed 's/seq4/LG05/g' tmp2 > tmp
sed 's/seq3/LG04/g' tmp > tmp2
sed 's/seq2/LG03/g' tmp2 > tmp
sed 's/seq1/LG02/g' tmp > tmp2
sed 's/seq0/LG01/g' tmp2 > tmp
sed 's/LG028/LG19/g' tmp > retrotransposons.withdomains.gff
```

Then bedtools getfasta was run to retrieve FASTA sequences:

```perl
bedtools getfasta -fi gadMor2.assembly.no_underscores.fasta.sr_masked.LGs -bed retrotransposons.withdomains.gff -fo retrotransposons.withdomains.fasta
```
The output _retrotransposons.withdomains.fasta_ is the final library of LTR retrotransposons to be used further.

##3. TransposonPSI

TransposonPSI uses PSIBLAST to search the genome against a library of transposon specific enzymes. NOTE: There will still be simple repeats designated as transposons, so it might be better to run it on a masked genome.

Running TransposonPSI:

```perl
module load blast/2.2.26
module load perlmodules
/projects/cees/bin/TransposonPSI/transposonPSI.pl genome.fa nuc
```
Time: 15 hours on a 600 mb genome

The output file of interest is _genome.fa.TPSI.allHits.chains.bestPerLocus.gff3_. Bedtools was used to extract FASTA sequences of found elements:

Here I loose a lot of potentially useful information, as the transposons detected in each sequence is lost using the getfasta command. Should find a way to transfer that information to the headers.

```perl
bedtools getfasta -fi genome.fa -bed genome.fa.TPSI.allHits.chains.bestPerLocus.gff3 -fo tPSI.fasta
```

_tPSI.fasta_ is the final output.

##4. RepeatModeler
###Build database
```perl
module load repeatmodeler/1.0.8
BuildDatabase -name repmoddatabase -engine ncbi genome.fa
```

The .out file is like this:


###Run RepeatModeler
```perl
module load repeatmodeler/1.0.8
RepeatModeler -database repmoddatabase
```
Time: 12 hours on a 600 mb genome.

A lot of files are produces, but the file of interest is _consensi.fa.classified_, containing all repeats discovered by RepeatModeler (including simple repeats).

##5. Classification of libraries using RepeatClassifier and TEclass.

The RepeatModeler library is classified by RepeatModeler, but the TransposonPSI and LTRharvest/LTRdigest libraries are not. 

First, run RepeatClassifier (a script part of RepeatModeler) on a merged version of the transposonPSI and LTRharvest/LTRdigest libraries:

```perl
cat tPSI.fasta retrotransposons.withdomains.fasta > tPSI.LTRdigest.lib.unclassifed
```
##RepeatClassifier
Run RepeatClassifier on the merged library:

```perl
module load repeatmodeler/1.0.8
RepeatClassifier -engine ncbi -consensi tPSI.LTRdigest.lib.unclassifed
```
Time: 8 hours. 

Output file: tPSI.LTRdigest.lib.unclassifed.classified.

##TEclass
Run TEclass on *tPSI.LTRdigest.lib.unclassifed.classified* and on *consensi.fa.classified*

```perl
/projects/cees/bin/TEclass/test_consensi_2.1.pl tPSI.LTRdigest.lib.unclassifed.classified
/projects/cees/bin/TEclass/test_consensi_2.1.pl consensi.fa.classified
```

Output file: tPSI.LTRdigest.lib.unclassifed.classified.lib, consensi.fa.classified.lib

Now a lot of elements in the library has two classifications. In order to avoid this, run the python script classification_chooser.py. The script merges RepeatClassifer and TEclass classifications, choosing RepeatClassifiers classification if there are disagreements. And prints to standard output. Send the output to a file. The output includes both classification, so strip all headers of '|TEclass*' using sed.

```python

#!/usr/bin/env python
"""
FINISHED, William Brynildsen, 6.10.2015
"""
from Bio import SeqIO
import sys, argparse, re

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=
	'This parser merges RepeatModeler and TEclass classifications, choosing RepeatModelers classification if there are disagreements')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	args = parser.parse_args()
		
for record in SeqIO.parse(args.input, 'fasta'):
	if "Unknown" in record.description:
		if "DNA" in record.description:
			print ">" + re.sub("Unknown", "DNA", record.description) + "\n" + record.seq
		elif "SINE" in record.description:
			print ">" + re.sub("Unknown", "SINE", record.description) + "\n" + record.seq
		elif "LINE" in record.description:
			print ">" + re.sub("Unknown", "LINE", record.description) + "\n" + record.seq
		elif "nonLTR" in record.description:
			print ">" + re.sub("Unknown", "nonLTR", record.description) + "\n" + record.seq
		elif "LTR" in record.description:
			print ">" + re.sub("Unknown", "LTR", record.description) + "\n" + record.seq
		elif "Retro" in record.description:
			print ">" + re.sub("Unknown", "Retro", record.description) + "\n" + record.seq
		else:
			print ">" + record.description + "\n" + record.seq
	else:
		print ">" + record.description + "\n" + record.seq
```
```perl
module load python2/2.7.9
python2 classification_chooser.py -i tPSI.LTRdigest.lib.unclassifed.classified.lib \
> tPSI.LTRdigest.lib.unclassifed.classified.chosen.lib

python2 classification_chooser.py -i consensi.fa.classified.lib > consensi.fa.classified.chosen.lib

#Removing parts of headers:
sed '/ .*//g' tPSI.LTRdigest.lib.unclassifed.classified.chosen.lib > tPSI.LTRdigest.lib.unclassifed.classified.chosen.stripped.lib

sed '/ .*//g' consensi.fa.classified.chosen.lib > consensi.fa.classified.chosen.stripped.lib
```

##6. Filtering for non-TE genes

First, concatenate the three libraries of RepeatModeler, LTRharvest/LTRdigest and TransposonPSI.

```perl
cat tPSI.LTRdigest.lib.unclassifed.classified.chosen.stripped.lib consensi.fa.classified.chosen.stripped.lib \
> AllRepeats.lib
```
###BLASTX against UniProt and RepeatPeps.lib
AllRepeats.lib are blasted against UniProt and RepeatPeps. If a sequence is found in UniProt, but not in RepeatPeps it is excluded from the library.

###BLASTX:

```perl
module load blast+/2.2.29

blastx -query ../AllRepeats.lib \
-db /cluster/software/repeatmasker/Libraries/RepeatPeps.lib -num_threads 10 -evalue 1e-10 -max_target_seqs 50 -outfmt 6 \
    -out AllRepeats.repeatpeps.blastx.1e-10.max50

blastx -query ../AllRepeats.lib \
  -db /projects/cees/in_progress/cod2/data/swiss_prot/uniprot_sprot.fasta -num_threads 10 -evalue 1e-10 -max_target_seqs 50 -outfmt 6 \
  -out AllRepeats.uniprot-sport.blastx.1e-10.max50
```
###Sorting:

```perl
sort -k1,1 -k12,12nr -k11,11n AllRepeats.repeatpeps.blastx.1e-10.max50 | \
sort -u -k1,1 --merge > AllRepeats.lib.repeatpeps.blastx.1e-10.max50.highest_scoring

sort -k1,1 -k12,12nr -k11,11n AllRepeats.uniprot-sport.blastx.1e-10.max50 | \
sort -u -k1,1 --merge > AllRepeats.lib.uniprot-sport.blastx.1e-10.max50.highest_scoring

cut -f 1 AllRepeats.lib.uniprot-sport.blastx.1e-10.max50.highest_scoring \
>  AllRepeats.lib.uniprot-sport.blastx.1e-10.max500.first_column
cut -f 1 AllRepeats.lib.repeatpeps.blastx.1e-10.max50.highest_scoring \
 > AllRepeats.lib.repeatpeps.blastx.1e-10.max50.first_column

echo "Num entries only hit to RepeatPeps" 
comm -13 AllRepeats.lib.uniprot-sport.blastx.1e-10.max500.first_column \
AllRepeats.lib.repeatpeps.blastx.1e-10.max50.first_column > AllRepeats.only_repeatpeps

cat AllRepeats.only_repeatpeps |wc -l
echo "Num entries only hit to UniProt"
comm -23 AllRepeats.lib.uniprot-sport.blastx.1e-10.max500.first_column \
AllRepeats.lib.repeatpeps.blastx.1e-10.max50.first_column > AllRepeats.only_uniprot

cat AllRepeats.only_uniprot |wc -l
echo "Num entries hit to both"
comm -12 AllRepeats.lib.uniprot-sport.blastx.1e-10.max500.first_column \
AllRepeats.lib.repeatpeps.blastx.1e-10.max50.first_column |wc -l
```
###Removing sequences only found in UniProt
```perl
module load python2/2.7.9
python2 reprint.filtered.lib.py -i AllRepeats.lib -l AllRepeats.only_uniprot > AllRepeats.chosen.stripped.filtered.lib
```

###reprint.filtered.lib.py:

```python
#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
        
        
        parser = argparse.ArgumentParser(description=
        'This parser will take a list and a library and reprint the library without the elements from the list')
        parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
        parser.add_argument('-l', '--list', action='store', help='', type=argparse.FileType('r'), default = '-')
        args = parser.parse_args()
                
        lst = []
        for i in args.list:
                lst.append(i.rstrip())
        for record in SeqIO.parse(args.input, 'fasta'):
                if record.id not in lst:
                        print ">" + record.id + "\n" + record.seq
```

##7. Clustering using USEARCH

###Filtering out simple repeats to avoid issues while clustering

```perl
python2 scripts/no_simple_repeats.py -i AllRepeats.chosen.stripped.filtered.lib > AllRepeats.chosen.stripped.filtered.no_sr.lib
```

The python script:

```python
#!/usr/bin/env python

"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
        
        
        parser = argparse.ArgumentParser(description=
        'This parser will take a list and a library and reprint the library without the elements from the list')
        parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')

        args = parser.parse_args()
                
        for record in SeqIO.parse(args.input, 'fasta'):
                if "Simple_repeat" not in record.id: 
			if "buffer" not in record.id:
                        	print ">" + record.id + "\n" + record.seq
```

###Checking the length of sequences

```perl
python2 scripts/lenght.of.sequences.py -i AllRepeats.chosen.stripped.filtered.no_sr.lib | sort -h | uniq -c |less
```

###Sorting for USEARCH

```
module load usearch

usearch -sortbylength ../AllRepeats.chosen.stripped.filtered.no_sr.lib --output AllRepeats.chosen.stripped.filtered.no_sr.usearch.sorted.lib --log usearch.log
```

###Running USEARCH

```
usearch -cluster_fast AllRepeats.chosen.stripped.filtered.no_sr.usearch.sorted.lib --id 0.8 --centroids my_centroids.fa --uc result.uc -consout final.nr.consensus.fa -msaout aligned.fasta --log usearch2.log
```

The output file final.nr.consensus.fa is the consensus sequences. Remove header text using:

```perl
sed 's/centroid=//g' final.nr.consensus.fa > final.nr.consensus.fa.stripped1
sed 's/;seqs=.*//g' final.nr.consensus.fa.stripped1 > final.nr.consensus.fa.stripped2
```

##8. Mask the genome using the final library

```
RepeatMasker -pa 16 -dir . -lib ../usearch/final.nr.consensus.fa.stripped2 ../gadMor2.assembly.no_underscores.fasta.masked 1>repmask_output 2>repmask_error
```



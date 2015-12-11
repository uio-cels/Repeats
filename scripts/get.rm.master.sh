#-----------------------Making small tables of the results-------------------------------#
export REPO=/work/users/willibr/testing/repository
GENOME=$1

sh $REPO/get.rm.stats.sh *.tbl > $GENOME.small_tbl

#   SINEs:            18366      1743314 bp    0.27 %
#   LINEs:            40122     12219712 bp    1.90 %
#   LTR elements:     57133     11610087 bp    1.80 %
#DNA transposons     234181     23106671 bp    3.59 %
#Unclassified:        13403      1094652 bp    0.17 %
#Simple repeats:     634033     57502326 bp    8.93 %
#Total interspersed repeats:    49774436 bp    7.73 %
#bases masked:  112098025 bp ( 17.42 %)

#	Cleaning up
sh $REPO/get.rm.stats.rename.sh $GENOME.small_tbl | column -t > $GENOME.final_table

#SINEs         18366      1743314   0.27
#LINEs         40122      12219712  1.90
#LTR           57133      11610087  1.80
#DNA           234181     23106671  3.59
#Unclassified  13403      1094652   0.17
#Simple        634033     57502326  8.93
#Transposons   49774436   7.73
#Total         112098025  17.42

sh $REPO/get.rm.overview.sh *.tbl > $GENOME.overview 


#	I want the "Simple", "Transposons" and "Total" in a separate table, as the columns of
#	"Transposons" and "Total" differ from the rest.

sh $REPO/get.rm.stats.rename.sh $GENOME.overview \
| column -t | awk '{print $1, $3, $4}' | grep "Simple" | column -t > simple

sh $REPO/get.rm.stats.rename.sh $GENOME.overview \
| column -t | grep "Transposons" | column -t > Transposons

sh $REPO/get.rm.stats.rename.sh $GENOME.overview \
| column -t | grep "Total" | column -t > Total

cat simple Transposons Total | column -t > $GENOME.final_overview

#--------------------------------Making GFF for IGV--------------------------------------#

#	NB! With the -gff option RepeatMasker produces a GFF2 file that does not work with
#	these commands, so make a GFF3 file using the RepeatMasker utility script 
#	mOutToGFF3.pl. 

/cluster/software/repeatmasker/util/rmOutToGFF3.pl $GENOME.out > $GENOME.gff3

#	Making files for IGV
#	Removing header
grep "[0-9]" $GENOME.out > $GENOME.out.no_header

#	Removing all "##sequence_region ..." headers in order to merge files
grep -v "##" $GENOME.gff3 > $GENOME.gff3.no_headers

#	Making a new GFF with columns from the .out and the .gff3 file, notice that columns
# 	$24 and $11 becomes element#class. The grep -v "*" chooses all elements that doesn't 
#	overlap with another element with a higher SW score and that are less than 80% of 
#	the length of the other element (reported by RepeatMasker).



paste $GENOME.out.no_header $GENOME.gff3.no_headers | grep -v "*" | \
awk '{print $5,"\t",$17,"\t",$18,"\t",$19,"\t",$20,"\t",$21,"\t",$22,"\t", \
$23,"\t",$24"#"$11";""Position=["$12":"$13":"$14"]"}' | sed 's/ //g' \
        > $GENOME.classified.no_overlaps.gff

#In IGV, first sort, second load and press OK to index.

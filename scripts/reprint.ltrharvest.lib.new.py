#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'Reprints FASTA files and changing the headers')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	parser.add_argument('-n', '--name')
	parser.add_argument('-t', '--type')
	args = parser.parse_args()
		
	counter = 0
	for record in SeqIO.parse(args.input, 'fasta'):
		print ">" + "%s-%s" % (args.name, counter) + "\n" + record.seq 
		counter += 1
	


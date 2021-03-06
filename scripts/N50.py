#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=
	'This parser will take a file containing a number (i.e. the N50 value), and only reprint sequences ≥ that number')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	parser.add_argument('-n', '--N50')
	args = parser.parse_args()
	
	f = open("N50")
	N50 = f.readline()
	f.close()
	for record in SeqIO.parse(args.input, 'fasta'):
		if len(record.seq) >= int(N50):
			print ">" + record.id + "\n" + record.seq

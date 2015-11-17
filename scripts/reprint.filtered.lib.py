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
	


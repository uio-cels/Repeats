#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	parser.add_argument('-c', '--classy', action='store', help='', type=argparse.FileType('r'), default = '-')
	args = parser.parse_args()
	
	d = args.classy
	classy = d.readlines()
	
	temp = ""
	
	for record in SeqIO.parse(args.input, 'fasta'):
		for line in classy:				
			if record.id in line and record.id not in temp:	
				temp = temp + ">" + line.rstrip() + "\n" + record.seq + "\n"
	print temp
	

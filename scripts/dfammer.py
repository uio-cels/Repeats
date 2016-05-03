#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'Input: nhmmscan output with -dfamtblout flag. Output: temporary file for dfammer2.py')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	
	args = parser.parse_args()
	tmp = ""
	
	for line in args.input:
		linesplit = line.split()
		if line.startswith("#") == False and linesplit[2] not in tmp and "Unknown" in line:
			tmp = tmp + line
	print tmp.rstrip()
	
	
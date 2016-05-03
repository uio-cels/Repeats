#!/usr/bin/env python
"""
FINISHED, William Brynildsen
"""

import sys, argparse

if __name__ == '__main__':
	
	
	parser = argparse.ArgumentParser(description=
	'Input: output of dfammer.py. Ouput: list of classified RepeatModeler repeats, not \
	classified by RepeatClassifier or by BLASTX homology to RepeatPeps.lib')
	parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')
	
	args = parser.parse_args()
	tmp = ""
	
	for line in args.input:
		linesplit = line.split()
		if "DNA" in line:
			print linesplit[2].split("#")[0] + "#" + "DNA" + "\t" + linesplit[2]
		if "SINE" in line:
			print linesplit[2].split("#")[0] + "#" + "SINE" + "\t" + linesplit[2]
		if "LINE" in line:
			print linesplit[2].split("#")[0] + "#" + "LINE" + "\t" + linesplit[2]
		if "non-LTR" not in line and "LTR" in line: 
			print linesplit[2].split("#")[0] + "#" + "LTR" + "\t" + linesplit[2]
		
	
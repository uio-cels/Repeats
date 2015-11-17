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
	args = parser.parse_args()
		
	counter = 0
	for record in SeqIO.parse(args.input, 'fasta'):
		if "cacta" in record.id:
			print ">" + "TransposonPSI-%s#DNA/CACTA" % (counter) + "\n" + record.seq 
			counter += 1
		if "DDE" in record.id:
			print ">" + "TransposonPSI-%s#DNA/DDE" % (counter) + "\n" + record.seq 
			counter += 1
		if "gypsy" in record.id:
			print ">" + "TransposonPSI-%s#LTR/Gypsy" % (counter) + "\n" + record.seq 
			counter += 1
		if "hAT" in record.id:
			print ">" + "TransposonPSI-%s#DNA/hAT" % (counter) + "\n" + record.seq 
			counter += 1
		if "helitron" in record.id:
			print ">" + "TransposonPSI-%s#RC/Helitron" % (counter) + "\n" + record.seq 
			counter += 1
		if "ISa" in record.id:
			print ">" + "TransposonPSI-%s#DNA/ISa" % (counter) + "\n" + record.seq 
			counter += 1
		if "ISb" in record.id:
			print ">" + "TransposonPSI-%s#DNA/ISb" % (counter) + "\n" + record.seq 
			counter += 1
		if "ISC1316" in record.id:
			print ">" + "TransposonPSI-%s#DNA/ISC1316" % (counter) + "\n" + record.seq 
			counter += 1
		if "LINE" in record.id:
			print ">" + "TransposonPSI-%s#LINE" % (counter) + "\n" + record.seq 
			counter += 1
		if "ltr_Roo" in record.id:
			print ">" + "TransposonPSI-%s#LTR/Pao" % (counter) + "\n" + record.seq 
			counter += 1
		if "mariner_ant1" in record.id:
			print ">" + "TransposonPSI-%s#DNA/TcMar-Ant1" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "mariner":
			print ">" + "TransposonPSI-%s#DNA/TcMar" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "MuDR_A_B":
			print ">" + "TransposonPSI-%s#DNA/MuLE-MuDR" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "P_element":
			print ">" + "TransposonPSI-%s#DNA/P" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "piggybac":
			print ">" + "TransposonPSI-%s#DNA/PiggyBac" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "TY1_Copia":
			print ">" + "TransposonPSI-%s#LTR/Copia" % (counter) + "\n" + record.seq 
			counter += 1
		if record.id == "Crypton":
			print ">" + "TransposonPSI-%s#DNA/Crypton" % (counter) + "\n" + record.seq 
			counter += 1
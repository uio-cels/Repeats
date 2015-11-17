"""
NOT FINISHED, William Brynildsen
"""

from Bio import SeqIO
import sys, argparse

if __name__ == '__main__':
        
        
        parser = argparse.ArgumentParser(description=
        'Input is a multiple fasta file, and output is a multiple fasta file with changed headers according to LTRdigest indexing')
        parser.add_argument('-i', '--input', action='store', help='', type=argparse.FileType('r'), default = '-')

        args = parser.parse_args()
        counter = 0
        for record in SeqIO.parse(args.input, 'fasta'): 
    		print ">" + "seq",counter, "\n" + record.seq
    		counter += 1
    		

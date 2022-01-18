#!/usr/bin/python
DESCRIPTION = '''
Get the frequency of each character (e.g. number of A, T, G and C's) in a fasta file. 
These numbers are required for calculationg the D2S statistic.
'''
import sys
import argparse
import logging
from D2S_tools import *
from itertools import groupby

## Pass arguments.
def main():
	# Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('--fasta', metavar='input.fasta', type=lambda x: read_file_check_compression(x), required=True, help='Input fasta file, can be gziped')
	parser.add_argument('--freq', metavar='input.fasta.chrFreq', type=lambda x: write_file_check_compression(x), required=True, help='Output character frequecny file, can be gziped')
	parser.add_argument('--debug', action='store_true', required=False, help='Print DEBUG info (default: %(default)s)')
	args = parser.parse_args()
	
	# Set up basic debugger
	if args.debug:
		logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s', stream=sys.stdout, level=logging.DEBUG)
	else:
		logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s', stream=sys.stdout, level=logging.ERROR)
	logger = logging.getLogger(__name__)
	
	logger.debug('%s', args) ## DEBUG
	
	char_set = {'A':0, 'C':0, 'G':0, 'T':0}
	charcter_freq_from_fasta(args.fasta, args.freq, char_set, logger)
	
	


def charcter_freq_from_fasta(in_fasta_fh, out_freq_fh, char_set, logger):
	'''
	Count the number of times each character in char_set are observed in the fasta file.
	Write the frequecny of each character (char_count / total_chars)
	
	Will also write "NUM_SEQUENCES" and "NUM_CHARACTERS" lines which details the 
	number of sequences and character (bases) in the fasta file.
	
	NUM_SEQUENCES and NUM_CHARACTERS is used by downstream scripts to calculate 
	the number of possible Kmers in a dataset. 
	'''
	
	# Count the number of sequences. Used later on for calculating the number of Kmers.
	totalSeqCount = 0
	
	for header, seq in fasta_iter(in_fasta_fh):
		totalSeqCount += 1
		seq.upper() # Convert seq to upper case.
		
		for char in char_set.keys():
			char_set[char] += seq.count(char)
	
	# Get total number of characters across all sequences. 
	totalChars = float(sum(char_set.itervalues()))
	
	for char in sorted(char_set.keys()):
		out_freq_fh.write(char + '\t' + str(char_set[char] / totalChars) + '\n')
	out_freq_fh.write("NUM_SEQUENCES" + '\t' + str(totalSeqCount) + '\n')
	out_freq_fh.write("NUM_CHARACTERS" + '\t' + str(totalChars) + '\n')



def fasta_iter(fh):
    	"""
    	Given a fasta file. yield tuples of header, sequence
	Clears description from seq name.
	
	From: https://www.biostars.org/p/710/
	Updated: 11/09/2018
	Version: 0.2
	"""
    	# ditch the boolean (x[0]) and just keep the header or sequence since
    	# we know they alternate.
    	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    	for header in faiter:
        	# drop the ">" and description
        	header = header.next()[1:].strip().split(' ')[0]
        	# join all sequence lines to one.
        	seq = "".join(s.strip() for s in faiter.next())
        	yield header, seq
	
	
	


if __name__ == '__main__':
	main()

#!/usr/bin/python
DESCRIPTION = '''
Pass a sorted kmer count file (kmer_seq<\\t>kmer_count) and prepare it for calculation of D2S.

This script converts the kmer sequence into a numeric representation. This is useful when you want
to compare two different kmer sets. 

NOTE - This acript assumes:
	- Thie kmers in the given file are sorted lexicographically (like the output for jellyfish).
	- The kmers are compoased of only [ATGC].

Input kmer file: kmer_seq<\\t>kmer_count
Output kmer file: kmer_value<\\t>kmer_seq<\\t>kmer_count
'''
import sys
from D2S_tools import *
import argparse
import logging

## Pass arguments.
def main():
	# Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-i', '--in_kmers', metavar='input.txt', type=lambda x: read_file_check_compression(x), default=sys.stdin, required=False, help='Input kmer count file (sorter), can be gziped (default: stdin)')
	parser.add_argument('-o', '--out_kmers', metavar='output.txt', type=lambda x: write_file_check_compression(x), required=True, help='Output kmer file, can be gziped')
	parser.add_argument('--debug', action='store_true', required=False, help='Print DEBUG info (default: %(default)s)')
	args = parser.parse_args()
	
	# Set up basic debugger
	if args.debug:
		logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s', stream=sys.stdout, level=logging.DEBUG)
	else:
		logging.basicConfig(format='#%(levelname)s :: %(asctime)s :: %(message)s', stream=sys.stdout, level=logging.ERROR)
	logger = logging.getLogger(__name__)
	
	
	logger.debug('%s', args) ## DEBUG
	
	char_mapping = {'A':'0', 'C':'1','G':'2','T': '3'} # Chracter mapping for ATGC only
	Kmers_2_NumbericRepresentation(args.in_kmers, args.out_kmers, char_mapping, logger=logger)
	
	
	
	
	
	
def Kmers_2_NumbericRepresentation(input_kmer_file, output_kmer_file, char_mapping, logger, sep='\t'):
	'''
	Convert kmer_seq into a numberic representation.
	Check that kmmer_value is > last kmer_value, thus file is sorted correctly. 
	'''
	
	last_kmer_value = -1
	
	for kmer_seq, kmer_count in pass_column_file(input_kmer_file, sep):
		kmer_seq.upper() # Convert to upper case.
		
		logger.debug('%s\t%s', kmer_seq, kmer_count) ## DEBUG
		
		kmer_value = int(''.join([char_mapping[base] for base in list(kmer_seq)]))
		
		logger.debug('kmer value: %s', kmer_value) ## DEBUG
		
		if kmer_value <= last_kmer_value: # Break if kmer has value <= last kmer
			sys.exit("Kmers not sorted")
		
		last_kmer_value = kmer_value
		
		output_kmer_file.write(str(kmer_value) + '\t' + kmer_seq + '\t' + kmer_count + '\n')



if __name__ == '__main__':
	main()

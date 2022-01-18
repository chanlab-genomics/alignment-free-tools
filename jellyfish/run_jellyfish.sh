#!/bin/bash

# Example extraction of 21-mers from the FASTA file sequence.fa
# and preparation of input files for distance calculations

file=sequence.fa
k=21

echo "## File:" $file
jellyfish count -m $k -s $s -t $NCPUS -o $file.$k.jf $file
jellyfish dump -ct $file.$k.jf | sort -k1,1 | python2 Kmers_2_NumericRepresentation.py -o $file.${k}mer.nkc.gz
python2 Composition_of_InputSeqs.py --fasta $file --freq $file.CharFreq
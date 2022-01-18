'''
Functions used by multiple scripts.
'''
import os
import sys
import gzip


def read_file_check_compression(arg):
    '''
    Check passed file name exists and opens using gzip when needed. 
    '''
    if not os.path.exists(arg):
        sys.exit('ERROR: The file %s does not exist!' % arg)
    else:
        # open with gzip if it has the *.gz extension
        if arg.endswith(".gz"):
            return gzip.open(arg, 'rb')
        else:
            return open(arg, 'r')


def check_file_exists(arg):
    '''
    Check file exists. Can open it later. 
    '''
    if not os.path.exists(arg):
        sys.exit('ERROR: The file %s does not exist!' % arg)
    else:
        return arg


def write_file_check_compression(arg):
    '''
    Opens file for writing, uses gzip when needed. 
    '''
    # open with gzip if it has the *.gz extension
    if arg.endswith(".gz"):
        return gzip.open(arg, 'wb')
    else:
        return open(arg, 'w')


def pass_column_file(fh, sep='\t'):
    '''
    Takes a file with columns seperated by 'sep' and yields each line:
    Will:
            - Strip new line characters (\n)
            - Split line using seperator (sep)
            - Ignores blank and comment lines

    Yields cleaned split line as a list.
    '''
    for line in fh:

        if isinstance(line, bytes):
            line = line.decode('utf-8')

        line = line.strip()

        if not line or line.startswith('#'):  # Ignore blank and comment lines
            continue

        yield line.split(sep)


def Next_Kmer(Kmer_iter):
    '''
    Iterate and return the next Kmer in the set. 
    Catchs the StopIteration returned by the pass_column_file generator and 
    returns a False value. Easiest way to handle instances where one file 
    finished and the other still has Kmers left. 
    '''
    try:
        try:
            value, seq, count = Kmer_iter.next()
        except AttributeError:
            value, seq, count = next(Kmer_iter)

        return int(value), seq, int(count), True
    except StopIteration:
        return None, None, None, False


def iterate_Kmer_sets(KmerSet1_fh, KmerSet2_fh, logger, Both_KmerSets=True, KmerSet1_Only=True, KmerSet2_Only=True):
    '''
    Iterates over two sorterd kmer files and returns kmers where:
            (1) Kmer is in KmerSet1 and KmerSet2 (Both_KmerSets)
            (1) Kmer is in KmerSet1 NOT KmerSet2 (KmerSet1_Only)
            (1) Kmer is in KmerSet2 NOT KmerSet1 (KmerSet2_Only)

    Assumes files are sorted lexicographically and the file was created by Kmers_2_NumbericRepresentation.py.
    '''
    report_interval = 10000000

    KmerSet1_fh.seek(0)
    KmerSet2_fh.seek(0)

    KmerSet1_notDone = True
    KmerSet2_notDone = True

    KmerSet1 = pass_column_file(KmerSet1_fh)
    KmerSet2 = pass_column_file(KmerSet2_fh)

    KmerSet1_value, KmerSet1_seq, KmerSet1_count, KmerSet1_notDone = Next_Kmer(
        KmerSet1)
    KmerSet2_value, KmerSet2_seq, KmerSet2_count, KmerSet2_notDone = Next_Kmer(
        KmerSet2)

    i = 0  # LOOP - Keep count
    while KmerSet1_notDone or KmerSet2_notDone:
        #logger.debug('Iteration %s: %s\t%s\t%s\t%s\t%s\t%s', i, KmerSet1_value, KmerSet1_seq, KmerSet1_count, KmerSet2_value, KmerSet2_seq, KmerSet2_count)

        # If Kmer is in BOTH datasets
        if KmerSet1_value == KmerSet2_value and KmerSet1_notDone and KmerSet2_notDone:
            logger.debug('%s\t%s\t%s\t%s\t%s\t%s', KmerSet1_value, KmerSet1_seq,
                         KmerSet1_count, KmerSet2_value, KmerSet2_seq, KmerSet2_count)

            if Both_KmerSets:
                yield KmerSet1_seq, KmerSet1_count, KmerSet2_seq, KmerSet2_count

            KmerSet1_value, KmerSet1_seq, KmerSet1_count, KmerSet1_notDone = Next_Kmer(
                KmerSet1)
            KmerSet2_value, KmerSet2_seq, KmerSet2_count, KmerSet2_notDone = Next_Kmer(
                KmerSet2)

        # If Kmer is in Dataset 1 ONLY
        elif KmerSet1_value < KmerSet2_value and KmerSet1_notDone:
            logger.debug('%s\t%s\t%s\t%s\t%s\t%s', KmerSet1_value,
                         KmerSet1_seq, KmerSet1_count, None, None, None)

            if KmerSet1_Only:
                yield KmerSet1_seq, KmerSet1_count, None, None

            KmerSet1_value, KmerSet1_seq, KmerSet1_count, KmerSet1_notDone = Next_Kmer(
                KmerSet1)

        # If Kmer is in Dataset 2 ONLY
        elif KmerSet1_value > KmerSet2_value and KmerSet2_notDone:
            logger.debug('%s\t%s\t%s\t%s\t%s\t%s', None, None, None,
                         KmerSet2_value, KmerSet2_seq, KmerSet2_count)

            if KmerSet2_Only:
                yield None, None, KmerSet2_seq, KmerSet2_count

            KmerSet2_value, KmerSet2_seq, KmerSet2_count, KmerSet2_notDone = Next_Kmer(
                KmerSet2)

        # If Dataset 2 is done BUT Dataset 1 still has some left
        elif KmerSet1_notDone and not KmerSet2_notDone:
            logger.debug('%s\t%s\t%s\t%s\t%s\t%s', KmerSet1_value,
                         KmerSet1_seq, KmerSet1_count, None, None, None)

            if KmerSet1_Only:
                yield KmerSet1_seq, KmerSet1_count, None, None

            KmerSet1_value, KmerSet1_seq, KmerSet1_count, KmerSet1_notDone = Next_Kmer(
                KmerSet1)

        # If Dataset 1 is done BUT Dataset 2 still has some left
        elif KmerSet2_notDone and not KmerSet1_notDone:
            logger.debug('%s\t%s\t%s\t%s\t%s\t%s', None, None, None,
                         KmerSet2_value, KmerSet2_seq, KmerSet2_count)

            if KmerSet2_Only:
                yield None, None, KmerSet2_seq, KmerSet2_count

            KmerSet2_value, KmerSet2_seq, KmerSet2_count, KmerSet2_notDone = Next_Kmer(
                KmerSet2)

        # LOOP - Print count
        i += 1
        if i % report_interval == 0:
            logger.info('K-mer loop round %s', i)  # INFO

        ## LOOP - Safety
        # if i > 1000:
        #	break

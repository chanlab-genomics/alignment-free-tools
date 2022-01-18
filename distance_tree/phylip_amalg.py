#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import sys
import shutil
import tarfile
import argparse
import pandas as pd
from glob import glob
from pprint import pprint
from contextlib import contextmanager
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple

CORRUPT_FILES: int = 0


@contextmanager
def change_dir(destination):
    try:
        cwd = os.getcwd()
        os.chdir(destination)
        yield
    finally:
        os.chdir(cwd)


def get_names(target_dir: str) -> list:
    """
    Gets a list of all the gene names from the a target directory which 
    holds all the comparison files.

    Parameters:
        A directory that holds all the result file named with the following
        convention:
            "id_name_1"-"id_name_2"."ext" (DON'T include paraenthesis in
            the actual filename, they are only used here for the purpose 
            of clarity and demonstration)

    Return:
        Returns a sorted list of all the unique gene names found.
    """

    # Get a list of all the files within the target directory
    target_files = glob(os.path.join(target_dir, "**"))

    # Just get the base directory of these target files
    target_files = map(os.path.basename, target_files)

    # Remove any filename extensions
    target_files = [file_.rsplit(
        '.', maxsplit=1)[0] if '.' in file_ else file_ for file_ in target_files]

    gene_names = []

    for target_file in target_files:
        gene_names.extend(target_file.split('-'))

    return sorted(list(set(gene_names)))


def get_gene_ids_from_path(result_path: str) -> List[str]:
    """
    Extarct the gene ids from a given result path.

    Return:
        Returns the extracted gene ids as a list of strings.
    """

    result_path = os.path.basename(result_path)

    result_path = result_path.rsplit('.', maxsplit=1)[
        0] if '.' in result_path else result_path

    return result_path.split('-')


def setup_df(name_list):
    """
    Creates a blank dataframe with the appropriate row and columns names.
    """

    blank_df = pd.DataFrame(index=name_list, columns=name_list, dtype=float)
    blank_df.fillna(0.0, inplace=True)

    return blank_df


def extract_result(result_path: str):
    """
    Extracts a single result from the specified result path.

    Parameters:
        result_path:
            A path to a single file containing a distance value.

    Returns:
        Returns the distance value (as a float) from the specified result path.
    """

    global CORRUPT_FILES

    with open(result_path, 'r') as result_file:
        result_str = result_file.read()

    # The result should be the the very end column after performing a split
    # on the semi-colon
    value = result_str.rsplit(';', maxsplit=1)[-1]

    try:
        return float(value)
    except ValueError:
        CORRUPT_FILES += 1
        #######################################################################
        # NOTE: Might want to change in the future!
        #######################################################################
        return 0


def poplate_single_result(phylip_df: pd.DataFrame, result_path: str):
    """
    Populates the dataframe with a single value from the results folder.

    Parameters:
        phylip_df:
            A dataframe that will eventually contain all the distance results.

        result_path:
            A path to a single file containing a distance value.
    """

    # Get the result value (as a float) from the specified result path
    result_value = extract_result(result_path)

    # Retrieve the gene ids from the file name
    gene_id_1, gene_id_2 = get_gene_ids_from_path(result_path)

    phylip_df.loc[gene_id_1][gene_id_2] = result_value
    phylip_df.loc[gene_id_2][gene_id_1] = result_value

    return


def populate_all_results(phylip_df: pd.DataFrame, result_dir: str):
    """
    Populates the dataframe with all results from a given result directory.

    Parameters:
        phylip_df:
            A dataframe that will eventually contain all the distance results.

        result_dir:
            A directory containing all the distance results.
    """

    # Get all the result files from the the result directory.
    target_files = glob(os.path.join(result_dir, "**"))

    for target_file in target_files:
        poplate_single_result(phylip_df, target_file)

    if CORRUPT_FILES > 0:
        print("[WARN] %d corrupted file/s found in %s (skipped)." % (CORRUPT_FILES, result_dir), file=sys.stderr)

    return


def print_phylip(phylip_df: pd.DataFrame, output_path: str):
    """
    Prints the input dataframe in standard PHYLIP format to the an output file.

    Parameter:
        phylip_df:
            A dataframe to be printed to the output directory.

        output_path:
            The output file path for the PHYLIP matrix.
    """

    # Get the number of columns from this dataframe
    rows, _ = phylip_df.shape

    # Left justify all of the indexes by 10
    phylip_df.rename(index=lambda index_: index_.ljust(10), inplace=True)

    with open(output_path, 'w', newline='') as output_file:

        # Write the number of rows/cols in the first line
        print('\t' + str(rows), end='\n', flush=True, file=output_file)

        # Write the remaining matrix, omit the column (header) names
        phylip_df.to_csv(output_file, sep='\t', float_format="%.8f", header=False,
                         index=True, index_label=False)

    return


def create_matrix(data_folder, output_file):

    zipped = data_folder.endswith(".tz.gz")
    zipped_data_folder = None

    if zipped:

        zipped_data_folder = data_folder
        data_folder = data_folder[:-len(".tz.gz")]

        dirname = os.path.dirname(zipped_data_folder)
        filename = os.path.basename(zipped_data_folder)

        with change_dir(dirname):
            tar = tarfile.open(filename, "r:gz")
            tar.extractall()
            tar.close()

    name_list = get_names(data_folder)
    blank_df = setup_df(name_list)
    populate_all_results(blank_df, data_folder)
    print_phylip(blank_df, output_file)

    if zipped:
        shutil.rmtree(data_folder)

    return


def main():

    parser = argparse.ArgumentParser(
        description="Creates a distance matrix from individual distance files.")

    parser.add_argument('--data', type=str, required=True,
                        help='A path to a directory or tarball that has the individual distances.')
    parser.add_argument('--matrix', type=str, required=True,
                        help='A path to a text file to dump the contents of the matrix.')

    args = parser.parse_args()

    create_matrix(args.data, args.matrix)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import sys
import csv
import time
import shutil
import socket
import itertools
import argparse

from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple
from glob import glob
from pprint import pprint
from string import Formatter
from datetime import timedelta

"""
Example Usage:
    (Unix)
python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny_red_40 --data_output_path /90days/s4430291/Genomes_for_AFphylogeny_red_40_1_D2S --temp True
python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny --data_output_path /90days/s4430291/Genomes_for_AFphylogeny_D2S --temp F --submit F
"""

if sys.platform.startswith('win32'):
    ROOT_DIR = os.path.join(
        'D:\\', '2020_SS', 'BioInfo', 'jackknifing')

elif sys.platform.startswith('linux'):
    ROOT_DIR = os.path.join(
        '/', 'home', 's4430291', 'chanlab-genomics', 'jackknifing')

    if not os.path.exists(ROOT_DIR):
        ROOT_DIR = os.path.join(
            '/', 'home', '564', 'mc7636', 'chanlab-genomics', 'jackknifing')

# Job time in minutes to run each python script
JOB_TIME = 9
JOB_MEM = "15GB"
# JOB_TIME = 0.8
# JOB_MEM = "12GB"
JOB_NODES = 1
NCPUS = 4

PYTHON_VERSION = "2.7"

if 'gadi' in socket.gethostname().lower():
    JOB_TEMPLATE = """#!/bin/bash
    #PBS -N {file_name}
    #PBS -j oe
    #PBS -o {stdout_file}
    #PBS -l ncpus={ncpus},mem={job_mem}
    #PBS -l walltime={job_time}
    #PBS -l jobfs=20MB

    #CHANGE THIS TO YOUR UQ-FACULTY-SCHOOL group name. 
    #USE the groups command to find out your exact group name. 
    #PBS -P d85 
    #PBS -q normalbw
    #PBS -l wd

    export OMP_NUM_THREADS={ncpus}

    DATE=$(date +"%d/%m/%Y %H:%M")
    echo "time started  "$DATE
    echo ------------------------------------------------------
    echo -n 'Job is running on the following nodes '; cat $PBS_NODEFILE
    echo ------------------------------------------------------
    echo PBS: qsub is running on $PBS_O_HOST
    echo PBS: originating queue is $PBS_O_QUEUE
    echo PBS: executing queue is $PBS_QUEUE
    echo PBS: working directory is $PBS_O_WORKDIR
    echo PBS: execution mode is $PBS_ENVIRONMENT
    echo PBS: job identifier is $PBS_JOBID
    echo PBS: job name is $PBS_JOBNAME
    echo PBS: node file is $PBS_NODEFILE
    echo PBS: current home directory is $PBS_O_HOME
    echo PBS: PATH = $PBS_O_PATH
    echo ------------------------------------------------------
    export TIMEFORMAT="%E sec"

    cd $PBS_O_WORKDIR
    pwd

    module load python
    {d2s_cmd}

    DATE=$(date +"%d/%m/%Y %H:%M")
    echo "time finished "$DATE
    """
else:
    JOB_TEMPLATE = """#!/bin/bash
    #PBS -N {file_name}
    #PBS -j oe
    #PBS -o {stdout_file}
    #PBS -l select=1:ncpus={ncpus}:mem={job_mem}
    #PBS -l walltime={job_time}

    #CHANGE THIS TO YOUR UQ-FACULTY-SCHOOL group name. 
    #USE the groups command to find out your exact group name. 
    #PBS -A NCMAS-d85

    export OMP_NUM_THREADS={ncpus}

    DATE=$(date +"%d/%m/%Y %H:%M")
    echo "time started  "$DATE
    echo ------------------------------------------------------
    echo -n 'Job is running on the following nodes '; cat $PBS_NODEFILE
    echo ------------------------------------------------------
    echo PBS: qsub is running on $PBS_O_HOST
    echo PBS: originating queue is $PBS_O_QUEUE
    echo PBS: executing queue is $PBS_QUEUE
    echo PBS: working directory is $PBS_O_WORKDIR
    echo PBS: execution mode is $PBS_ENVIRONMENT
    echo PBS: job identifier is $PBS_JOBID
    echo PBS: job name is $PBS_JOBNAME
    echo PBS: node file is $PBS_NODEFILE
    echo PBS: current home directory is $PBS_O_HOME
    echo PBS: PATH = $PBS_O_PATH
    echo ------------------------------------------------------
    export TIMEFORMAT="%E sec"

    cd $PBS_O_WORKDIR
    pwd

    module load python
    {d2s_cmd}

    DATE=$(date +"%d/%m/%Y %H:%M")
    echo "time finished "$DATE
    """

JOB_TEMPLATE = '\n'.join(map(str.strip, JOB_TEMPLATE.split(sep='\n')))


def convert_bool_arg(arg_in):
    """
    Converts a boolean command line argument into a python boolean object.

    Parameters:
        arg_in:
            The given command line argument.

    Return:
        The interrupted boolean value of the command line argument.
    """

    false_options = ('f', 'false', '0', 'n', 'no', 'nan', 'none')

    if arg_in.lower() in false_options:
        return False

    return True


def strfdelta(tdelta, fmt='{H:02}:{M:02}:{S:02}', inputtype='timedelta'):
    """Convert a datetime.timedelta object or a regular number to a custom-
    formatted string, just like the stftime() method does for datetime.datetime
    objects.

    The fmt argument allows custom formatting to be specified.  Fields can 
    include seconds, minutes, hours, days, and weeks.  Each field is optional.

    Some examples:
        '{D:02}d {H:02}h {M:02}m {S:02}s' --> '05d 08h 04m 02s' (default)
        '{W}w {D}d {H}:{M:02}:{S:02}'     --> '4w 5d 8:04:02'
        '{D:2}d {H:2}:{M:02}:{S:02}'      --> ' 5d  8:04:02'
        '{H}h {S}s'                       --> '72h 800s'

    The inputtype argument allows tdelta to be a regular number instead of the  
    default, which is a datetime.timedelta object.  Valid inputtype strings: 
        's', 'seconds', 
        'm', 'minutes', 
        'h', 'hours', 
        'd', 'days', 
        'w', 'weeks'
    """

    # Convert tdelta to integer seconds.
    if inputtype == 'timedelta':
        remainder = int(tdelta.total_seconds())
    elif inputtype in ['s', 'seconds']:
        remainder = int(tdelta)
    elif inputtype in ['m', 'minutes']:
        remainder = int(tdelta)*60
    elif inputtype in ['h', 'hours']:
        remainder = int(tdelta)*3600
    elif inputtype in ['d', 'days']:
        remainder = int(tdelta)*86400
    elif inputtype in ['w', 'weeks']:
        remainder = int(tdelta)*604800

    f = Formatter()
    desired_fields = [field_tuple[1] for field_tuple in f.parse(fmt)]
    possible_fields = ('W', 'D', 'H', 'M', 'S')
    constants = {'W': 604800, 'D': 86400, 'H': 3600, 'M': 60, 'S': 1}
    values = {}
    for field in possible_fields:
        if field in desired_fields and field in constants:
            values[field], remainder = divmod(remainder, constants[field])
    return f.format(fmt, **values)


class JobCreator:
    "Creates (and possibly runs) job scripts for creating annotated images."

    def __init__(self, slurm_dir: str, data_input_path: str, data_output_path: str,
                 groups: int = 50, index: int = 0, submit: bool = False, temp: bool = False, dry_run: bool = False):
        """
        Initializes a job creator.

        Parameters:
            slurm_dir (str):
                A full path to a directory that will hold stdout and stderr of 
                the submitted jobs.

            submit (bool):
                If True the created jobs will be immediately submitted.

            temp (bool):
                An argument to specif if the batch files should only be 
                temporarily held.

            dry_run (bool):
                If True, creates the batch files for the jobs and simulates
                job submission.
        """

        self.slurm_dir = slurm_dir
        self.data_input_path = data_input_path
        self.data_output_path = data_output_path

        # How many distance calculations should be run in a single batch script
        self.groups: int = groups
        self.index: int = index

        if not os.path.exists(self.data_output_path):
            os.makedirs(self.data_output_path)

        # A full path to a directory that will hold the created job files.
        self.output_dir = os.path.join(self.slurm_dir, "batch_scripts")

        self.slurm_out = os.path.join(self.slurm_dir, "batch_out")
        self.slurm_err = os.path.join(self.slurm_dir, "batch_err")

        # Make sure new directories exist
        for folder in [self.output_dir, self.slurm_out, self.slurm_err]:
            if not os.path.exists(folder):
                os.makedirs(folder)

        self.submit = submit
        self.temp = temp
        self.dry_run = dry_run

        # Get all the different job argument combinations
        self.job_args = self.get_job_arg_combinations()

        self.begin_job_procession()

    def get_job_arg_combinations(self):
        """
        Infers argument job combnations from the data path.
        """

        job_args: List[Dict[str, str]] = []

        # Get all the fasta files from the data path
        fasta_files: List[str] = []

        # A list of all valid fasta file extensions
        fasta_ext = (".fasta", ".fna", ".ffn", ".faa", ".frn", ".fas")

        for ext in fasta_ext:
            fasta_files.extend(
                glob(os.path.join(self.data_input_path, '**' + ext)))

        # Create combinations for each different pair
        for kmerset1, kmerset2 in itertools.combinations(fasta_files, 2):

            new_arg_dict = dict()

            new_arg_dict["kmerset1"] = kmerset1 + ".21mer.nkc.gz"
            new_arg_dict["kmerset2"] = kmerset2 + ".21mer.nkc.gz"

            new_arg_dict["kmerset1_freq"] = kmerset1 + ".CharFreq"
            new_arg_dict["kmerset2_freq"] = kmerset2 + ".CharFreq"

            # Only added it to the job argument list if the kemer analysis
            # file and CharFreq file exist for both kemer groups
            if all(os.path.exists(path) for path in new_arg_dict.values()):

                # Create a name for the otuput file
                kmerset1_name, _ = os.path.basename(
                    kmerset1).rsplit('.', maxsplit=1)
                kmerset2_name, _ = os.path.basename(
                    kmerset2).rsplit('.', maxsplit=1)

                output_path = kmerset1_name + '-' + kmerset2_name + '.txt'
                output_path = os.path.join(self.data_output_path, output_path)

                new_arg_dict["D2S_out"] = output_path

                job_args.append(new_arg_dict)

        return job_args

    def begin_job_procession(self):
        """
        Creates, submits (if requested) and removes (if requested) jobs for
        all the tiles within the project directory.
        """

        # Begin by creating all the required job files
        self.create_job_files()

        # Submit all the created jobs
        if self.submit:
            self.submit_jobs()

        # Delete all the created jobs
        if self.temp:
            self.delete_jobs()

        return

    def create_job_files(self):
        """
        Creates a job files for each tile within the project folder.
        """

        # Create a list of all the commands that need to be run to compute the
        # distances
        d2s_cmds = [' '.join(f'--{param_name} {param_value}' for param_name, param_value in job_arg.items())
                    for job_arg in self.job_args]
        d2s_cmds = [
            "python{py_ver_short} -W ignore {python_filepath!r} {param_str}".format(
                py_ver_short=PYTHON_VERSION.rsplit(".", maxsplit=1)[0],
                python_filepath=os.path.join(
                    ROOT_DIR, "calc_d2s", "Calculate_D2S.py"),
                param_str=param_str) for param_str in d2s_cmds
        ]

        d2s_cmds = [d2s_cmds[i:i + self.groups]
                    for i in range(0, len(d2s_cmds), self.groups)]

        for param_id, d2s_cmd in enumerate(d2s_cmds, 0):

            # Create a timedelta object to format the amount of time we
            # for this job
            job_time: timedelta = timedelta(minutes=JOB_TIME * len(d2s_cmd))

            # Create a string of all the parameter names with their
            # corresponding parameter values.
            d2s_cmd: List[str] = [d2s_cmd[i:i + NCPUS]
                                  for i in range(0, len(d2s_cmd), NCPUS)]
            d2s_cmd = map(' & \n'.join, d2s_cmd)
            d2s_cmd: str = ' & \nwait\n'.join(d2s_cmd)
            d2s_cmd += ' & \nwait'

            file_name: str = f"d2s_{self.index}_{param_id}"

            stdout_path = os.path.join(
                self.slurm_out, f"{file_name}_out.txt")
            stderr_path = os.path.join(
                self.slurm_err, f"{file_name}_err.txt")

            FORMATTED_TEMPLATE = JOB_TEMPLATE.format(
                file_name=file_name,
                stdout_file=stdout_path,
                stderr_file=stderr_path,
                d2s_cmd=d2s_cmd,
                job_time=strfdelta(job_time),
                job_mem=JOB_MEM,
                job_nodes=JOB_NODES,
                ncpus=NCPUS
            )

            job_filename = f"{file_name}_job.sh"
            job_path = os.path.join(self.output_dir, job_filename)

            with open(job_path, "w+") as job_file:
                print(FORMATTED_TEMPLATE, end='',
                      file=job_file, flush=True)

        return

    def submit_jobs(self):
        """
        Submit the jobs via slurm on the current machine.
        """

        # Make a list of all the job files within the output directory
        job_file_list = [filename for filename in os.listdir(
            self.output_dir) if filename.endswith(".sh")]
        # Get the full path to each of the filenames
        job_file_list = list(map(lambda filename: os.path.join(
            self.output_dir, filename), job_file_list))

        job_file_list = sorted(job_file_list)

        # Submit each of the jobs within the job file list
        for job_file in job_file_list:

            if self.dry_run:
                print(f"[DRY RUN] qsub {job_file}")
            else:
                os.system(f"qsub {job_file}")
                time.sleep(0.5)

        return

    def delete_jobs(self):
        """
        Deletes all the job files within the output directory.
        """

        shutil.rmtree(self.output_dir)

        return


def main():

    parser = argparse.ArgumentParser(description="Creates (and possibly runs) "
                                     "job scripts for the jackknife workflow.")

    parser.add_argument('--slurm_dir', type=str, default=os.path.join(ROOT_DIR, "batch", "d2s_jobs"),
                        help='A full path to a directory to create slurm and batch files.')
    parser.add_argument('--data_input_path', type=str, required=True,
                        help='A full path to the nkc.gz and CharFreq files.')
    parser.add_argument('--data_output_path', type=str, required=True,
                        help='An output folder for the d2s script.')
    parser.add_argument('--group', type=int, required=False, default=500,
                        help='Indicates how many distance calculations are run in a single batch script.')
    parser.add_argument('--index', type=int, required=False, default=1,
                        help='Indicates the index sample.')

    parser.add_argument('-s', '--submit', type=convert_bool_arg, default=False, const=True, nargs='?',
                        help='If True the created jobs will be immediately submitted.')
    parser.add_argument('-t', '--temp', type=convert_bool_arg, default=False, const=True, nargs='?',
                        help='If True the created job folder will be deleted immediately after submitting the jobs.')
    parser.add_argument('-d', '--dry_run', type=convert_bool_arg, default=False, const=False, nargs='?',
                        help='If True the program will simulate job submission output text but will not submit the jobs.')

    args = parser.parse_args()

    JobCreator(args.slurm_dir, args.data_input_path, args.data_output_path,
               index=args.index, groups=args.group, submit=args.submit, temp=args.temp, dry_run=args.dry_run)


if __name__ == "__main__":

    main()

    exit(0)

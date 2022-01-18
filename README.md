# Jackknifing Workflow

## Jackknifing

To build our jackknife tree we first need to produce some jackknife samples from our genome data (surprise surprise!). The jackknifing step mostly makes use of `jackknife.py` found in the top level directory. In short `jackknife.py` reads in a fasta file a spits out a reduced version. `jackknife.py --help` does a pretty good job of explaining what command line arguments it's expecting, so I've taken the liberty of copying and pasting the output of running `jackknife.py` with the `--help` flag here
```
usage: jackknife.py [-h] --input_paths INPUT_PATHS [INPUT_PATHS ...] [--output_path OUTPUT_PATH] [-v] [--portion PORTION] [--chunk_size CHUNK_SIZE] [--threads THREADS]

Randomly removes a portion of data from a fasta file.

optional arguments:
  -h, --help            show this help message and exit
  --input_paths INPUT_PATHS [INPUT_PATHS ...]
                        A list of fasta files to perform the jacknifing process
  --output_path OUTPUT_PATH
                        A folder to write the reduced fasta files.
  -v, --verbose         Include to run the script in verbose mode. Useful for checking progress. WARN: This may cause the program to run slower.
  --portion PORTION     The portion of data to be removed. The default is a 40 percent reduction.
  --chunk_size CHUNK_SIZE
                        The size of the chunks that get randomly removed from sequences. Default is a chunk size of 100.
  --threads THREADS     The number of threads used to run the jackknife algorithm. If 0 threads are specified then it will default to os.cpu_count().
```
Here's an example use of `jackknife.py`
```
python3 jackknife.py --input_path ~/AEH.fasta --output_path ~/jn_yeast --portion=40 --chunk_size=100 --threads=4
```
If you're going through this workflow by yourself, chances are you've got hundereds of genomes to jackknife across many different samples. To see how this process could be automated check `batch\array_jn\jn_array.sh` and `batch\array_jn\submit_jn_wait.sh`. The `batch\array_jn\jn_array.sh` file is a PBS job which produces a jackknife sample for each fasta file found in `ARRAY_TARGET` and saves the outputs to `OUTPUT_DIR`. You might want to change around the number of cpus, memory and walltime for the job file depending on how big/how many fasta files you have. The `batch\array_jn\submit_jn_wait.sh` will automatically submit multiple `batch\array_jn\jn_array.sh` jobs with different sample values.

## Jellyfish

The next step in this workflow is to use jellyfish as well as a few custom purpose python scripts to perform some analysis on our jackknifed samples. This means you will either have to install the `jellyfish` program on your computer, or on you're working on a hpc, you could try running `module load jellyfish/*version*` or just `module load jellyfish` to load jellyfish (to check if they have jellyfish at all run `module avail` to view all availed modules). Once jellyfish has been setup you will need to run the following string of commands on each jacknified fasta file
```

k=*kmer-count* # eg k=21
s=10000000000

file=*jackknifed fasta file* # eg file=~/AEH_red_40.fasta
jellyfish count -m $k -s $s -t $CPUS_PER_TASK -o $file.$k.jf $file && \
jellyfish dump -ct $file.$k.jf | sort -k1,1 | python2 jf_scripts/Kmers_2_NumbericRepresentation.py -o $file.${k}mer.nkc.gz && \
python2 jf_scripts/Composition_of_InputSeqs.py --fasta $file --freq $file.CharFreq && \
touch $file.done &
```
Note that these commands usually require a lot of memory and are therefore best run as batch jobs. Again doing this for each individual jackknifed file can be very tedious. To see how this process can be automated, take a look at `batch\jellyfish_jobs\submit_jellyfish.sh` and `batch\jellyfish_jobs\run_jellyfish.sh`. The `batch\jellyfish_jobs\run_jellyfish.sh` is a PBS job file which runs the above string of commands on every fasta file found within `ARRAY_TARGET`. The `batch\jellyfish_jobs\submit_jellyfish.sh` bash script simply submits the `batch\jellyfish_jobs\run_jellyfish.sh` batch script for different samples. A successful completion of these commands for a fasta file should have produced a `.*kmer-count*.jf` `.*kmer-count*mer.nkc.gz`, `.CharFreq` and a `.done` file. For example here's what a directory containing just `~/AEH_red_40.fasta` with `*kmer-count*=21` would look like
```
AEH_red_40.fasta               
AEH_red_40.fasta.21.jf         
AEH_red_40.fasta.21mer.nkc.gz  
AEH_red_40.fasta.CharFreq      
AEH_red_40.fasta.done
```

## Distance Calculations

The end goal is to create a distance tree, so first we need to compute all the pair-wise distances between each jackknifed sample. Computing the distance can be done using `jf_scripts\Calculate_D2S_mod.py`. Again, `jf_scripts\Calculate_D2S.py --help` does a pretty good job of explaining what command line arguments it's expecting, so here's it's output
```
usage: Calculate_D2S.py [-h] --kmerset1 KmerSet1.21mers.gz --kmerset1_freq
                        KmerSet1.21mers.charFreq --kmerset2 KmerSet2.21mers.gz
                        --kmerset2_freq KmerSet2.21mers.charFreq
                        [--D2S_out D2S.txt] [--debug]

Calculate the D2S score between two Kmer sets.

optional arguments:
  -h, --help            show this help message and exit
  --kmerset1 KmerSet1.21mers.gz
                        Kmers for dataset 1, can be gziped
  --kmerset1_freq KmerSet1.21mers.charFreq
                        Character frequency for dataset 1, can be gziped
  --kmerset2 KmerSet2.21mers.gz
                        Kmer for dataset 2, can be gziped
  --kmerset2_freq KmerSet2.21mers.charFreq
                        Character frequency for dataset 2, can be gziped
  --D2S_out D2S.txt     Output for D2S score (default: <open file '<stdout>',
                        mode 'w' at 0x2b8adbdca150>)
  --debug               Print DEBUG info (default: False)
```
Here is an example of using `Calculate_D2S.py` between the outputs of `AEH_red_40.fasta` and `AEG_red_40.fasta`
```
python2 ./calc_d2s/Calculate_D2S.py --kmerset1 AEG_red_40.fasta.21mer.nkc.gz --kmerset1_freq AEG_red_40.fasta.CharFreq --kmerset2 AEH_red_40.fasta.21mer.nkc.gz --kmerset2_freq AEH_red_40.fasta.CharFreq --D2S_out AEG-AEH.txt
```
If run successfully, this  commands should produce a single text file `AEG-AEH.txt`, containing some preamble and then the distance value preceded by a semi-colon. Again, doing this for each individual combination will be extremely time consuming. For this reason `.\calc_d2s\create_d2s_jobs.py` has been created to

1. Automatically create job files for every combination (and creates output distance file names that are compatible with `phylip_amalg.py`)
2. Submit the created job files
3. Remove the created job files once submitted, if desired (as many hundreds of job files can be created accross multiple samples)

The output of `python3 .\calc_d2s\create_d2s_jobs.py --help` is given below
```
usage: create_d2s_jobs.py [-h] [--slurm_dir SLURM_DIR] --data_input_path DATA_INPUT_PATH --data_output_path DATA_OUTPUT_PATH [--group GROUP] [--index INDEX] [-s [SUBMIT]] [-t [TEMP]] [-d [DRY_RUN]]

Creates (and possibly runs) job scripts for the jackknife workflow.

optional arguments:
  -h, --help            show this help message and exit
  --slurm_dir SLURM_DIR
                        A full path to a directory to create slurm and batch files.
  --data_input_path DATA_INPUT_PATH
                        A full path to the nkc.gz and CharFreq files.
  --data_output_path DATA_OUTPUT_PATH
                        An output folder for the d2s script.
  --group GROUP         Indicates how many distance calculations are run in a single batch script.
  --index INDEX         Indicates the index sample.
  -s [SUBMIT], --submit [SUBMIT]
                        If True the created jobs will be immediately submitted.
  -t [TEMP], --temp [TEMP]
                        If True the created job folder will be deleted immediately after submitting the jobs.
  -d [DRY_RUN], --dry_run [DRY_RUN]
                        If True the program will simulate job submission output text but will not submit the jobs.
```
Note that `.\calc_d2s\create_d2s_jobs.py` assumes that all of the jellyfish outputs are stored in the same directory. Here's an example of `.\calc_d2s\create_d2s_jobs.py` in action
```
python3 calc_d2s/create_d2s_jobs.py --data_input_path ~/sample_1 --data_output_path ~/sample_1_D2S --temp T --submit T --dry_run F --index=1
```
For our above example, once all the jobs have completed after running the above example we should find the file `AEG-AEH.txt` in the directory `~/sample_1_D2S`.

## Distance Tree Creation

Now for the part we've all been waiting for ... creating the distance tree! First however, we're going to need to make a distance matrix. Of course, you could manually to this yourself but this can be time consuming and is very prone to error. Instead if you have all of your distance files in the same directory with the file name format `[Gene name 1]-[Gene name 1].txt` (make sure that none of your gene names are more than 10 characters long!!) you can run `PHYLIP/phylip_amalg.py` on the folder to automatically generate the distance matrix for you. Here's the output of running `python3 PHYLIP/phylip_amalg.py --help`
```
usage: phylip_amalg.py [-h] --data DATA --matrix MATRIX

Creates a distance matrix from individual distance files.

optional arguments:
  -h, --help       show this help message and exit
  --data DATA      A path to a directory or tarball that has the individual distances.
  --matrix MATRIX  A path to a text file to dump the contents of the matrix.
```
Pretty self explanatory. Once you have your distance matrix you will need to convert this into a distance tree. You can do this using the `neighbour` program found in the PHYLIP suite (see: https://evolution.genetics.washington.edu/phylip/getme-new1.html). This can be a bit painful to use since the program will prompt you for arguments. Instead I've created a modified version of the PHYLIP's neighbor program where you only need to specify the important arguments as command line arguments, see the git hub page for my modified version: https://github.com/Michae1CC/PHYLIP-neighbor

The modified PHYLIP's neighbor program has the following invocation (and won't prompt you for any other arguments once running)
```
neighbor [input-matrix-path] [outtree-path] [outfile-path]
```

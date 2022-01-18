#!/bin/bash
#-#PBS -l nodes=1:ppn=2,mem=60GB,vmem=60GB,walltime=24:00:00
#PBS -k oe
#PBS -j oe
# NOTE: This Job array line needs to change!
#PBS -J 0-196
#PBS -N D2S_jellyfish
#PBS -l select=1:ncpus=2:mpiprocs=2:mem=60GB:vmem=60GB,walltime=24:00:00

#PBS -j oe

#PBS -A NCMAS-d85

set -o errexit
cd $PBS_O_WORKDIR
#source ~/run_cmd.sh

module load python
module load jellyfish/2.2.10

k=21

# Load file.
# n=0
# while read fline
# do
#         n=$((n + 1))
#         filearray[$n]=$fline
# done < /home/s4430291/chanlab-genomics/jackknifing/jf_scripts/files2run/files2runH3H4.txt

filearray=(/30days/s4430291/Genomes_for_AFphylogeny_red_40_2/**.fna)

s=10000000000

# Run a given line in the array
file=${filearray[$PBS_ARRAY_INDEX]}
echo "## File:" $file
jellyfish count -m $k -s $s -t $NCPUS -o $file.$k.jf $file && \
jellyfish dump -ct $file.$k.jf | sort -k1,1 | python2 /home/s4430291/chanlab-genomics/jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o $file.${k}mer.nkc.gz && \
python2 /home/s4430291/chanlab-genomics/jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta $file --freq $file.CharFreq && \
touch $file.done



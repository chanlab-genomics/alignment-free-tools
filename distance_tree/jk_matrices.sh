#!/bin/bash

if [ $# -eq 2 ]
then

  for INDEX in `seq $1 1 $2`;
do
    echo "Generating matrix ${INDEX}"
    python3 ~/chanlab-genomics/jackknifing/PHYLIP/phylip_amalg.py --data /scratch/d85/mc7636/Yeast/D2S_archive/Genomes_for_AFphylogeny_red_40_${INDEX}_D2S.tz.gz --matrix /scratch/d85/mc7636/Yeast/jk_matrices/mat_${INDEX}.txt &
done

else
  echo "USAGE: jk_matrices.sh [START] [END]"
fi

wait

echo "Completed matrix generation ($1-$2)"
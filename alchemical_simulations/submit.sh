#!/bin/bash -x
#SBATCH -A nmdar
#SBATCH --nodes=1
#SBATCH -t 02:00:00

for file in $(ls | grep .json); do
  mkdir "${file%.*}_job"

  for repeat in {0..2}; do
    srun -N 1 -c 8 --cpu-bind=threads --distribution=block:cyclic:cyclic --threads-per-core=1 openfe quickrun $file -o ${file%.*}_job/results_${file%.*}/results_${repeat} -d ${file%.*}_job/results_${file%.*}/results_${repeat} &
  done

done

wait

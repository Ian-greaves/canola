#!/usr/bin/env bash

#SBATCH --job-name=Canola_FASTQC
#SBATCH --time=1:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mem=32GB

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load fastqc

SAMPLES=( $(cut -d " " -f 1 fastqcInputList.csv) );
INFILES=( $(cut -d " " -f 2 fastqcInputList.csv) );
 
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    fastqc /OSM/CBR/AF_HETEROSIS/work/2018_dataschool/data/${INFILES[$i]} -o /flush1/che249/Canola_data/fastqc/
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi





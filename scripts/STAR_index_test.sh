#!/usr/bin/env bash

#SBATCH --job-name=CanolaIndex
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --exclusive
#SBATCH --mem=32GB

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load star

#Working_directories /OSM/CBR/AF_HETEROSIS/work
LOGDIR=/flush1/che249/Canola_data/log
INPDIR=/OSM/CBR/AF_HETEROSIS/work/2018_dataschool/data
REFDIR=/OSM/CBR/AF_HETEROSIS/work/2018_dataschool/reference
ANODIR=/OSM/CBR/AF_HETEROSIS/work/2018_dataschool/annotation
OUTDIR=/flush1/che249/Canola_data/processed

STAR \
--runMode genomeGenerate \
--genomeDir ${REFDIR}/ \
--genomeFastaFiles ${REFDIR}/Brassica_napus_v4.1.chromosomes.fa \
--sjdbGTFfile ${ANODIR}/Brassica_napus.annotation_v5_exon.gtf \
--sjdbOverhang 149

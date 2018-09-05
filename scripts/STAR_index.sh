#!/usr/bin/env bash

#SBATCH --job-name=CanolaIndex
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mem=32GB

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load star

#Working_directories /flush1/gre486/Canola_data/work
LOGDIR=/flush1/gre486/Canola_data/log
INPDIR=/flush1/gre486/Canola_data/raw
REFDIR=/flush1/gre486/Canola_data/reference_genome
ANODIR=/flush1/gre486/Canola_data/annotation
OUTDIR=/flush1/gre486/Canola_data/processed

STAR \
--runMode genomeGenerate \
--genomeDir ${REFDIR}/ \
--genomeFastaFiles ${REFDIR}/Brassica_napus_v4.1.chromosomes.fa \
--sjdbGTFfile ${ANODIR}/Brassica_napus.annotation_v5_exon.gtf \
--sjdbOverhang 149
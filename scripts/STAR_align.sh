#!/usr/bin/env bash

#SBATCH --job-name=CanolaAlign
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mem=40GB

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load star

#Working_directories /flush1/gre486/Canola_data/work
INPDIR=/flush1/gre486/Canola_data/raw
REFDIR=/flush1/gre486/Canola_data/reference_genome
ANODIR=/flush1/gre486/Canola_data/annotation
OUTDIR=/flush1/gre486/Canola_data/processed/Alignment

#STAR \
#--runThreadN 8 \
#--genomeDir ${REFDIR}/ \
#--readFilesCommand zcat \
#--readFilesIn ${INPDIR}/C_4_1_1_1.fq.gz ${INPDIR}/C_4_1_1_2.fq.gz \
#--outFileNamePrefix ${OUTDIR}/C_4_1_1/ \   #gives standard prefix to all output files
#--outSAMstrandField intronMotif \
#--outFilterMismatchNmax 10 \  #maximum number of mismatches per pair, large number switches oﬀ this ﬁlter
#--outSAMtype BAM SortedByCoordinate \    #outputs a bam file which is sorted  
#--quantMode TranscriptomeSAM GeneCounts \    #produces gene counts 
#--outFilterMultimapNmax 10 \   #max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
#--outSAMattrIHstart 0 \   #required to be used with Stringtie which is new version of cufflinks?
#--outSAMmapqUnique 255 \  # value used to select uniqe read.  If value is the same for 10 reads it will randomly select the best
#--outSAMmultNmax -1 \   # only outputs one alignment??
#--chimSegmentMin 40  #allows chimeric reads to be mapped.  Minimum size is 40

SAMPLES=( $(cut -d " " -f 1 STARInputList.csv) );
INFILES_R1=( $(cut -d " " -f 2 STARInputList.csv) );
INFILES_R2=( $(cut -d " " -f 3 STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    STAR \
	--runThreadN 8 \
	--genomeDir ${REFDIR}/ \
	--readFilesCommand zcat \
	--readFilesIn ${INPDIR}/${INFILES_R1[$i]} ${INPDIR}/${INFILES_R2[$i]} \
	--outFileNamePrefix ${OUTDIR}/${SAMPLES[$i]} \
	--outSAMstrandField intronMotif \
	--outFilterMismatchNmax 10 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFilterMultimapNmax 10 \
	--outSAMattrIHstart 0 \
	--outSAMmapqUnique 255 \
	--outSAMmultNmax -1 \
	--chimSegmentMin 40
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi


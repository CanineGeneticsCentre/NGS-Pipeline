#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake-himem

#SBATCH -o logs/mergeBam_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}

source ${SAMPLE}.config

module load ${GATK}
module load ${SAMTOOLS}

PG_ID=$(samtools view -H ${DIR}/${SAMPLE}.L${LANE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 2 | sed 's/ID://')
PG_PN=$(samtools view -H ${DIR}/${SAMPLE}.L${LANE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 3 | sed 's/PN://')
PG_VN=$(samtools view -H ${DIR}/${SAMPLE}.L${LANE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 4 | sed 's/VN://')
PG_CL=$(samtools view -H ${DIR}/${SAMPLE}.L${LANE}.aligned.bam | grep '@PG' | grep 'ID:bwa' | cut -f 5 | sed 's/CL://')

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" MergeBamAlignment \
  --REFERENCE_SEQUENCE ${FASTA}/${GENOME}.fasta \
  --UNMAPPED_BAM ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
  --ALIGNED_BAM ${DIR}/${SAMPLE}.L${LANE}.aligned.bam \
  --OUTPUT ${DIR}/${SAMPLE}.L${LANE}.merged.bam \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --SORT_ORDER "unsorted" \
  --PROGRAM_RECORD_ID "${PG_ID}" \
  --PROGRAM_GROUP_VERSION "${PG_VN}" \
  --PROGRAM_GROUP_COMMAND_LINE "${PG_CL}" \
  --PROGRAM_GROUP_NAME "${PG_PN}" \
  --CLIP_ADAPTERS false \
  --UNMAP_CONTAMINANT_READS true \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ATTRIBUTES_TO_REMOVE NM \
  --ATTRIBUTES_TO_REMOVE MD \
  --CREATE_INDEX false \
  --ADD_MATE_CIGAR true \
  --IS_BISULFITE_SEQUENCE false \
  --ALIGNED_READS_ONLY false \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --TMP_DIR ${HOME}/hpc-work/tmp/


input_size=$(stat -c%s "${DIR}/${SAMPLE}.L${LANE}.unaligned.bam")
output_size=$(stat -c%s "${DIR}/${SAMPLE}.L${LANE}.merged.bam")

if [ output_size > input_size ]; then
  mv ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam ${DIR}/${SAMPLE}.L${LANE}.aligned.bam ${DIR}/${SAMPLE}.L${LANE}.fastq.gz tmp_files/
else
  exit 1;
fi
#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time 08:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/alignFastq_%A-%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
#LANE=$SLURM_ARRAY_TASK_ID
DIR=rg${SLURM_ARRAY_TASK_ID}

source ${SAMPLE}.config

#module load ${BWA}
#module load ${SAMTOOLS}

${BWA} mem -K 100000000 -p -v 3 -t 16 -Y ${FASTA}/${GENOME}.fasta ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.fastq.gz | ${SAMTOOLS} view -h -b -o ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.aligned.bam

#bwa-mem2 mem -K 100000000 -p -v 3 -t 16 -Y ${FASTA}/${GENOME}.fasta ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.fastq.gz | samtools view -h -b -o ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.aligned.bam

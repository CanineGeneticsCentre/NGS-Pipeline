#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake-himem

#SBATCH -o logs/mark-adapters_%A-%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
#LANE=$SLURM_ARRAY_TASK_ID
DIR=rg${SLURM_ARRAY_TASK_ID}
source ${SAMPLE}.config

module load ${GATK}


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" MarkIlluminaAdapters \
  --INPUT ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.unaligned.bam \
  --OUTPUT ${DIR}/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.adaptMarked.bam \
  --METRICS metrics/${SAMPLE}.RG${SLURM_ARRAY_TASK_ID}.adaptMarked.metrics.txt

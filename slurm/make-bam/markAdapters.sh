#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake-himem

#SBATCH -o logs/mark-adapters_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load ${GATK}


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" MarkIlluminaAdapters \
  --INPUT ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
  --OUTPUT ${DIR}/${SAMPLE}.L${LANE}.adaptMarked.bam \
  --METRICS metrics/${SAMPLE}.L${LANE}.adaptMarked.metrics.txt

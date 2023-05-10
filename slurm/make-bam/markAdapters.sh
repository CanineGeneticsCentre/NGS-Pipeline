#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=BEGIN,FAIL,INVALID_DEPEND,END
#SBATCH -p skylake

#SBATCH -o logs/mark-adapters_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load ${GATK}


gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" MarkIlluminaAdapters \
  --INPUT ${DIR}/${SAMPLE}.L${LANE}.unaligned.bam \
  --OUTPUT ${DIR}/${SAMPLE}.L${LANE}.adaptMarked.bam \
  --METRICS metrics/${SAMPLE}.L${LANE}.adaptMarked.metrics.txt

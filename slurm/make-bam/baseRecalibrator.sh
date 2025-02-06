#!/usr/bin/env bash

#! RUN : sbatch sortSam.sh <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/baseRecal-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
source ${SAMPLE}.config
n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx3G" BaseRecalibrator \
  -R ${FASTA}/${GENOME}.fasta \
  -I ${SAMPLE}.sorted.bam \
  --use-original-qualities \
  -O base_recal/${SAMPLE}.${n}.bsqr.txt \
  --known-sites ${BQSR} \
  -L intervals/${n}-scattered.interval_list 

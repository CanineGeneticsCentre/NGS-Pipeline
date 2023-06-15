#!/usr/bin/env bash

#! RUN : sbatch applyBSQR.sh <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o logs/applyBQSR-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
source ${SAMPLE}.config
n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" ApplyBQSR \
  -R ${FASTA}/${GENOME}.fasta \
  -I ${SAMPLE}.sorted.bam \
  -O base_recal/${SAMPLE}.${n}.bam \
  -L intervals/${n}-scattered.interval_list \
  -bqsr ${SAMPLE}.bsqr.out \
  --static-quantized-quals 10 \
  --static-quantized-quals 20 \
  --static-quantized-quals 30 \
  --add-output-sam-program-record \
  --create-output-bam-md5 \
  --use-original-qualities 

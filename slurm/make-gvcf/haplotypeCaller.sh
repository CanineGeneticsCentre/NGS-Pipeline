#!/usr/bin/env bash

#! RUN : sbatch haplotypeCaller.sh <SAMPLE>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 06:00:00
#SBATCH --mail-type=FAIL
##SBATCH --no-requeue
#SBATCH -p cclake

#SBATCH -o logs/haplotypeCaller-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
REF=$2
PCR_MODEL=$3
source ${SAMPLE}.config

module load ${GATK}

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" HaplotypeCaller \
  -R ${FASTA}/${GENOME}.fasta \
  -I ${SAMPLE}-${REF}.bam \
  -L intervals/${n}-scattered.interval_list \
  -O gvcf/${SAMPLE}-${REF}.${n}.g.vcf \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  -ERC GVCF \
  -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
  --pcr-indel-model ${PCR_MODEL}

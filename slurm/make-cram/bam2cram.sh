#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 08:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake

#SBATCH -o /home/%u/hpc-work/logs/bam2cram_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1

#source ${SAMPLE}.config
source /home/es904/scripts/NGS-Pipeline/ngs-pipeline-cf4.config

rm -rf ${SAMPLE}-cf4.cram

${SAMTOOLS} view -T ${FASTA}/${GENOME}.fasta -C -o ${SAMPLE}-cf4.cram ${SAMPLE}-cf4.bam
${SAMTOOLS} index ${SAMPLE}-cf4.cram

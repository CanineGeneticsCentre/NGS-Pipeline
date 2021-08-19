#!/usr/bin/env bash

#! RUN : sbatch fastq2sam.sh <SAMPLE>>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=16
#! How much wallclock time will be required?
#SBATCH --time 12:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o ../logs/job-%A_%a.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load bwa/0.7.12                      # bwa
module load samtools/1.10                   # samtools

SAMPLE=$1
source ../${SAMPLE}.config


FILE1=$(ls ../*s_${SLURM_ARRAY_TASK_ID}.r_1.fq.gz)
FILE2=$(ls ../*s_${SLURM_ARRAY_TASK_ID}.r_2.fq.gz)

#bwa mem -M -t 16 ${FASTA}/${GENOME}.fasta ${FILE1} ${FILE2} | samtools view -1 - > ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam
bwa mem -M -t 16 ${FASTA}/${GENOME}.fasta ${FILE1} ${FILE2} | samtools sort -o ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam

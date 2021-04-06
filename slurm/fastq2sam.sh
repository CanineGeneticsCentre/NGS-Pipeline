#!/usr/bin/env bash

#! RUN : sbatch fastq2sam.sh <TAG> <CHR>

#! sbatch directives begin here ###############################
#! Which project should be charged:
#SBATCH -A MELLERSH-SL3-CPU
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

#SBATCH -o logs/job-%j.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load bwa/0.7.12                      # bwa
module load samtools/1.10                   # samtools

SAMPLE=$1
source ${SAMPLE}.config

bwa mem -M -t 16 ${FASTA}/${GENOME}.fasta ${SAMPLE}\_R1.fastq.gz ${SAMPLE}\_R2.fastq.gz | samtools view -1 - > ${SAMPLE}.aligned.bam
#!/usr/bin/env bash

#! RUN : sbatch -A MELLERSH-SL3-CPU miseqStats.sh <FILE1> <FILE2> <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=4
#! How much wallclock time will be required?
#SBATCH --time 00:05:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=BEGIN,FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o /home/%u/hpc-work/logs/job-%j.out



. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load bwa-0.7.17-gcc-5.4.0-42mry2g                      # bwa
module load samtools-1.15-gcc-5.4.0-elpvwsy                   # samtools

FILE1=$1
FILE2=$2
SAMPLE=$3

source /home/es904/scripts/NGS-Pipeline/ngs-pipeline-cf4.config

#bwa mem -M -t 16 ${FASTA}/${GENOME}.fasta ${FILE1} ${FILE2} | samtools view -1 - > ${SAMPLE}.s_${SLURM_ARRAY_TASK_ID}.aligned.bam
bwa mem -M -t 4 ${FASTA}/${GENOME}.fasta ${FILE1} ${FILE2} | samtools sort -o ${SAMPLE}.aligned.bam

samtools flagstat ${SAMPLE}.aligned.bam > ${SAMPLE}.flagstat.out

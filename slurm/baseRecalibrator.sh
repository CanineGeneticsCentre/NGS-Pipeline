#!/usr/bin/env bash

#! RUN : sbatch sortSam.sh <SAMPLE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 00:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake
#SBATCH --mem=5gb

#SBATCH -o ../logs/job-%A_%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load jdk-8u141-b15-gcc-5.4.0-p4aaopt 
module load gatk/4.1.0.0                    # GATK 4.1

SAMPLE=$1
source ${SAMPLE}.config


intervals=`head -${SLURM_ARRAY_TASK_ID} sequence_grouping.txt | tail -1 | sed s/"\t"/" -L "/g`

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx2G" BaseRecalibrator -R ${FASTA}/${GENOME}.fasta -I ${SAMPLE}.sorted.bam --use-original-qualities -O ${SAMPLE}.${SLURM_ARRAY_TASK_ID}.bsqr.txt --known-sites ${BQSR} -L ${intervals}

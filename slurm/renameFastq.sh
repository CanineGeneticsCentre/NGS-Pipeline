#!/usr/bin/env bash

#! RUN : sbatch fastq2sam.sh <DIR>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 02:00:00
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


source ../${SAMPLE}.config

DIR=$1
FILE=`head -${SLURM_ARRAY_TASK_ID} files.list | tail -1`

mv ${FILE} ${DIR}_${SLURM_ARRAY_TASK_ID}.fq
gzip ${DIR}_${SLURM_ARRAY_TASK_ID}.fq

mv ${DIR}_${SLURM_ARRAY_TASK_ID}.fq.gz ${WGS}/${SAMPLE}/
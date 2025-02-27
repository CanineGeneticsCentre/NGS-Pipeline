#!/usr/bin/env bash

#! RUN : sbatch fastq2sam.sh <DIR>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 03:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL,INVALID_DEPEND,END
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p cclake

#SBATCH -o ../logs/job-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

#source ../${SAMPLE}.config

DIR=$1
FILE=`head -${SLURM_ARRAY_TASK_ID} files.list | tail -1`

echo $FILE

mv ${FILE} ${DIR}.s_${SLURM_ARRAY_TASK_ID}.fq
gzip ${DIR}.s_${SLURM_ARRAY_TASK_ID}.fq

mv ${DIR}.s_${SLURM_ARRAY_TASK_ID}.fq.gz ../

echo ${DIR}.s_${SLURM_ARRAY_TASK_ID}.fq.gz

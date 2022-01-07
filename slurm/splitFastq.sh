#!/usr/bin/env bash

#! RUN : sbatch splitFastq.sh <FILE>

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

#SBATCH -o logs/job-%j.out

module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

SAMPLE=$1
FILE=$2

DIR="${FILE%%.*}"
SCRIPTS=`dirname $0`

source ${SAMPLE}.config

mkdir $DIR
cd $DIR

#split $f into chunks of approx 5000,000,000 lines
zcat ../$FILE | split -l500000000 --additional-suffix=.fq

FQ=`ls | wc -l`
sbatch -A ${ACCOUNT} -J renameFastq --array=1-${FQ} ${SCRIPTS}/slurm/renameFastq.sh ${DIR}

cd ../
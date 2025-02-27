#!/usr/bin/env bash

#! RUN : sbatch splitFastq.sh <FILE>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=16
#! How much wallclock time will be required?
#SBATCH --time 02:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=BEGIN,FAIL,INVALID_DEPEND
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p cclake

#SBATCH -o logs/job-%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE=$1
FILE=$2

DIR="${FILE%%.*}"

#source ${SAMPLE}.config

mkdir $DIR
cd $DIR

#split $f into chunks of approx 5000,000,000 lines
zcat ../$FILE | split -l500000000 --additional-suffix=.fq

ls -1 x*.fq > files.list
FQ=`wc -l files.list | cut -f1 -d ' '`

cat files.list

#echo sbatch -A ${ACCOUNT} -J renameFastq --array=1-${FQ} --export=SAMPLE=${SAMPLE} ${SCRIPTS}/slurm/renameFastq.sh ${DIR}
sbatch -J ${SAMPLE}.renameFastq --array=1-${FQ} --export=SAMPLE=${SAMPLE} ${SCRIPTS}/slurm/renameFastq.sh ${DIR} 

cd ../

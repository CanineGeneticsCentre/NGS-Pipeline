#!/usr/bin/env bash

#! RUN : sbatch createGenomicsDB.sh <SAMPLE_MAP> <REF>

#! sbatch directives begin here ###############################
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4   ## hopefully gives me 4 cores on the same node, therefore 4 x 6840MB of RAM... (maybe!)
#SBATCH --time 08:00:00
#SBATCH --mail-type=ALL
#SBATCH -p cclake-himem

#SBATCH --output=logs/job-%A_%a.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE_MAP=$1
REF=$2
source ${REF}.config

module load ${GATK}

CHR_COUNT=`wc -l ${FASTA}/intervals/gdb-chromsomes.intervals | awk '{print $1}'`
standard_opts=" --batch-size 50 --reader-threads 4 --genomicsdb-shared-posixfs-optimizations "

if [[ ${SLURM_ARRAY_TASK_ID} -le ${CHR_COUNT} ]]; then 
  # Use gdb-chromsomes.intervals
  INTERVALS=`head -${SLURM_ARRAY_TASK_ID} ${FASTA}/intervals/gdb-chromsomes.intervals | tail -1 | sed s/" "/" -L "/g`
  DIR=`echo ${INTERVALS} | cut -f 1 -d' ' | cut -d'_' -f 1 | cut -f 1 -d':'`
  if [[ ${#DIR} -lt 4 ]] ; then
    DIR="chr"${DIR}
  fi
else 
  #Use gdb-scaffolds.intervals
  LINE=`expr ${SLURM_ARRAY_TASK_ID} - ${CHR_COUNT}`
  DIR=`head -${LINE} ${FASTA}/intervals/gdb-scaffolds.intervals | tail -1`
  INTERVALS="${FASTA}/intervals/${DIR}.list"
  standard_opts+="--merge-contigs-into-num-partitions 1 "
fi

rm -rf ${GDB}/${DIR}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx6G" GenomicsDBImport \
    ${standard_opts} \
    --sample-name-map ${SAMPLE_MAP} \
    --tmp-dir ${HOME}/hpc-work/tmp/ \
    --genomicsdb-workspace-path ${GDB}/${DIR} \
    -L ${INTERVALS}

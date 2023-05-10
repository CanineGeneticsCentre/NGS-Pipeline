#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake-himem

#SBATCH -o logs/sort-bam_%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

SAMPLE=$1
source  ${SAMPLE}.config

module load $GATK

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" SortSam  \
  --INPUT ${SAMPLE}.aligned.unsorted.dedup.bam \
  --OUTPUT ${SAMPLE}.sorted.bam \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE false \
  --SORT_ORDER "coordinate" \
  --VALIDATION_STRINGENCY SILENT
  --TMP_DIR ${HOME}/hpc-work/tmp/

if [ ! -f ${SAMPLE}.sorted.bam ]; then
	exit 1;
else
	bam_size=$(wc -c < ${SAMPLE}.sorted.bam)
	if [ $bam_size -ge 50000000 ]; then
		rm -rf ${SAMPLE}.aligned.unsorted.dedup.bam
	fi
fi

#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time 12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH -p cclake-himem
#SBATCH --mem=25000

#SBATCH -o logs/sort-bam_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment

SAMPLE=$1
source  ${SAMPLE}.config

module load $GATK

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx24G" SortSam  \
  --INPUT ${SAMPLE}.aligned.unsorted.dedup.bam \
  --OUTPUT ${SAMPLE}.sorted.bam \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE false \
  --SORT_ORDER "coordinate" \
  --VALIDATION_STRINGENCY SILENT \
  --TMP_DIR ${HOME}/hpc-work/tmp/

if [ ! -f ${SAMPLE}.sorted.bam ]; then
	exit 1;
else
  input_size=$(stat -c%s "${SAMPLE}.aligned.unsorted.dedup.bam")
  output_size=$(stat -c%s "${SAMPLE}.sorted.bam")
  if [ output_size > input_size ]; then
    #mv ${SAMPLE}.aligned.unsorted.dedup.bam tmp_files/
    rm -rf ${SAMPLE}.aligned.unsorted.dedup.bam
  else
    exit 1;
  fi
#	bam_size=$(wc -c < ${SAMPLE}.sorted.bam)
#	if [ $bam_size -ge 50000000 ]; then
#		rm -rf ${SAMPLE}.aligned.unsorted.dedup.bam
#	fi
fi

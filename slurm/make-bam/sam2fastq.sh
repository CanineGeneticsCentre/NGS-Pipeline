#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake

#SBATCH -o logs/sam2fastq_%A-%a.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment


SAMPLE=$1
LANE=$SLURM_ARRAY_TASK_ID
DIR=lane${LANE}
source ${SAMPLE}.config

module load ${GATK}

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx4G" SamToFastq \
  --INPUT ${DIR}/${SAMPLE}.L${LANE}.adaptMarked.bam \
  --FASTQ ${DIR}/${SAMPLE}.L${LANE}.fastq.gz \
  --CLIPPING_ATTRIBUTE XT \
  --CLIPPING_ACTION 2 \
  --INTERLEAVE true \
  --NON_PF true \
  --TMP_DIR ${HOME}/hpc-work/tmp/


if [ ! -f ${DIR}/${SAMPLE}.L${LANE}.fastq.gz ]; then
	exit 1;
fi

FILES=(../*.s_${LANE}.*fq.gz)
((input_size = $(stat -c%s ${FILES[0]}) + $(stat -c%s ${FILES[1]})))
output_size=$(stat -c%s "${DIR}/${SAMPLE}.L${LANE}.fastq.gz")
if [ output_size > input_size ]; then
  mv ${DIR}/${SAMPLE}.L${LANE}.adaptMarked.bam tmp_files/
fi
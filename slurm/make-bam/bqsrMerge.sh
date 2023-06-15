#!/usr/bin/env bash

#! RUN : sbatch mergeBam.sh <SAMPLE> <INTERVALS>

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH -p cclake

#SBATCH -o logs/bqsrMerge-%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge
module load rhel7/default-ccl


SAMPLE=$1
INTERVALS=$2
REF=$3
source ${SAMPLE}.config

module load ${GATK}

BAMS=""
for i in `seq 0 $(($INTERVALS-1))`; do n=$(printf "%04d" $i); BAMS+="-I base_recal/${SAMPLE}.$n.bam "; done

gatk --java-options  "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx3G" GatherBamFiles ${BAMS} -O /dev/stdout | \
gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx3G" SortSam \
  --INPUT /dev/stdin  \
  --OUTPUT ${SAMPLE}-${REF}.bam \
  --SORT_ORDER "coordinate" \
  --TMP_DIR ${HOME}/hpc-work/tmp/ \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE true


if [ ! -f ${SAMPLE}-${REF}.bam ]; then
  echo "ERROR - ${SAMPLE}-${REF}.bam does not exist";
	exit 1;
fi

input_size=$(stat -c%s "${SAMPLE}.sorted.bam")
output_size=$(stat -c%s "${SAMPLE}-${REF}.bam")
if [ output_size > input_size ]; then
  rm -rf base_recal;
  rm -rf ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.bai ${SAMPLE}.bsqr.out
  #mv ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.bai ${SAMPLE}.bsqr.out tmp_files/
fi

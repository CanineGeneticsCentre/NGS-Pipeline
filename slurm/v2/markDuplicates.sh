#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time 06:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,INVALID_DEPEND
#SBATCH -p skylake-himem

#SBATCH -o logs/mark-duplicates_%j.out

. /etc/profile.d/modules.sh                 # Leave this line (enables the module command)
module purge                                # Removes all modules still loaded
module load rhel7/default-peta4             # REQUIRED - loads the basic environment

module load gatk-4.2.5.0-gcc-5.4.0-hzdcjga

SAMPLE=$1
LANES=$2

INPUT=''
for LANE in `seq 1 ${LANES}`; do 
  INPUT+=" --INPUT lane${LANE}/${SAMPLE}.L${LANE}.merged.bam"
done

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" MarkDuplicates ${INPUT} \
  --OUTPUT ${SAMPLE}.aligned.unsorted.dedup.bam \
  --METRICS_FILE metrics/${SAMPLE}.duplicate_metrics \
  --VALIDATION_STRINGENCY SILENT \
  --ASSUME_SORT_ORDER "queryname" \
  --CLEAR_DT "false" \
  --READ_NAME_REGEX null \
  --TMP_DIR ${HOME}/hpc-work/tmp/

#!/usr/bin/env bash

#! sbatch directives begin here ###############################
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time 06:00:00
#SBATCH --mail-type=FAIL
#SBATCH -p cclake-himem
#SBATCH --mem=12000

#SBATCH -o logs/mark-duplicates_%j.out

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-ccl              # REQUIRED - loads the basic environment


SAMPLE=$1
READ_GROUPS=$2
source ${SAMPLE}.config

module load ${GATK}


INPUT=''
for RG in `seq 1 ${READ_GROUPS}`; do 
  INPUT+=" --INPUT rg${RG}/${SAMPLE}.RG${RG}.merged.bam"
  ((input_size+=$(stat -c%s "rg${RG}/${SAMPLE}.RG${RG}.merged.bam")))
done

gatk --java-options "-Djava.io.tmpdir=${HOME}/hpc-work/tmp/ -Xmx10G" MarkDuplicates ${INPUT} \
  --OUTPUT ${SAMPLE}.aligned.unsorted.dedup.bam \
  --METRICS_FILE metrics/${SAMPLE}.duplicate_metrics \
  --VALIDATION_STRINGENCY SILENT \
  --ASSUME_SORT_ORDER "queryname" \
  --CLEAR_DT "false" \
  --READ_NAME_REGEX null \
  --TMP_DIR ${HOME}/hpc-work/tmp/

# If output file from MarkDuplicates is larger than the sum of the input files, DELETE input files
output_size=$(stat -c%s "${SAMPLE}.aligned.unsorted.dedup.bam")
if [ output_size > input_size ]; then
  for RG in `seq 1 ${READ_GROUPS}`; do 
    rm -rf rg${RG}/${SAMPLE}.RG${RG}.merged.bam
    rm -rf rg${RG}
  done
fi
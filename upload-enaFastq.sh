#!/usr/bin/env bash

#! RUN : bash upload-enaFastq.sh <SAMPLE>
#! Eg. : bash upload-enaFastq.sh BoT_37630

SAMPLE=$1
SCRIPTS=`dirname $0`

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided to upload"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE

# Count how many fastq files we have
COUNT=`ls ${WGS}/${SAMPLE}/*.fastq.gz | wc -l`
if [ $COUNT -le 1 ]; then
  rsync --progress -auvh ${WGS}/${SAMPLE}/*.fq.gz ./
  cat *.r_1.fq.gz > ${SAMPLE}_R1.fastq.gz
  cat *.r_2.fq.gz > ${SAMPLE}_R2.fastq.gz
else
  rsync --progress -auvh ${WGS}/${SAMPLE}/*.fastq.gz ./
fi


sbatch -J ${SAMPLE}.ascp --export=SAMPLE=${SAMPLE},CONDA_PREFIX,ASPERA_SCP_PASS ${SCRIPTS}/slurm/ena-ascp.sh


if [ $COUNT -le 1 ]; then
  md5sum *.fastq.gz > ${SAMPLE}.fastq.md5
  for f in `ls *.fastq.gz ${SAMPLE}.fastq.md5`; do 
    sync --progress -auvh0. $f ${WGS}/${SAMPLE}/;
  done
fi
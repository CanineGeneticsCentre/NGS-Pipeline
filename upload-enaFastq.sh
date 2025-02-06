#!/usr/bin/env bash

#! RUN : bash upload-enaFastq.sh <SAMPLE>
#! Eg. : bash upload-enaFastq.sh BoT_37630

SAMPLE_LIST=$1
SCRIPTS=`dirname $0`
export CONDA_PREFIX='/rds/project/rds-Qr3fy2NTCy0/Software/local/miniconda3'

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: No SAMPLE_LIST provided to upload"; exit 1; }

for SAMPLE in `cat ${SAMPLE_LIST}`; do
	if [[ -f ${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz && -f ${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz ]]; then
		#ls -lh ${SAMPLE}/fastq*
		continue;
	fi
  mkdir -p $SAMPLE/fastq; cd $SAMPLE
  # Count how many fastq files we have
  COUNT=`ls ${WGS}/${SAMPLE}/*.fastq.gz | wc -l`
  if [ $COUNT -le 1 ]; then
    rsync --progress -auvh --no-links ${WGS}/${SAMPLE}/*.fq.gz ./
    cat *.r_1.fq.gz > fastq/${SAMPLE}_R1.fastq.gz
    cat *.r_2.fq.gz > fastq/${SAMPLE}_R2.fastq.gz
#    rm -rf *.r_1.fq.gz *.r_2.fq.gz
  else
    rsync --progress -auvh --no-links ${WGS}/${SAMPLE}/*.fastq.gz fastq
  fi
  if [[ ! -f fastq/${SAMPLE}.fastq.md5 ]]; then
    cd fastq; md5sum *.fastq.gz > ${SAMPLE}.fastq.md5; 
    cat ${SAMPLE}.fastq.md5; echo
    scp ${SAMPLE}.fastq.md5 kcgc:/home/kcgc/tmp/ena/fastq/
    cd ../../
  fi
done;

SAMPLES=$(wc -l ${SAMPLE_LIST} | cut -f1 -d' ')

JOB_NAME=$(basename $SAMPLE_LIST .list)
#sbatch -J enaFastq.${JOB_NAME} --array=1-${SAMPLES}%4 --export=SAMPLE_LIST=${SAMPLE_LIST},CONDA_PREFIX,ASPERA_SCP_PASS ${SCRIPTS}/slurm/ena-ftp.sh fastq
sbatch -J enaFastq.${JOB_NAME} --array=1-${SAMPLES}%4 --export=SAMPLE_LIST=${SAMPLE_LIST},CONDA_PREFIX,ASPERA_SCP_PASS ${SCRIPTS}/slurm/ena-ascp.sh fastq



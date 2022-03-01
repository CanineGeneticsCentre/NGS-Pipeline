#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE> <REF>
#! Eg. : bash ngs-pipeline.sh CS_35365 cf4

SAMPLE=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $SAMPLE.config; source $SAMPLE.config

# Copy files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/*.fastq.gz ./

# Count how many fastq files we have
COUNT=`ls *.fastq.gz | wc -l`
# Work out number of lanes used... assuming pair-end, therefore divide total number of files by 2 - forward and reverse files for each lane
LANES=$((COUNT / 2))

if [ $LANES -le 2 ];then
  echo "Need to split FASTQ files"
  for f in `ls *.fastq.gz`; do 
    echo sbatch -A ${ACCOUNT} -J splitFastq --export=SCRIPTS=${SCRIPTS} ${SCRIPTS}/slurm/splitFastq.sh ${SAMPLE} ${f}
    sbatch -A ${ACCOUNT} -J splitFastq --export=SCRIPTS=${SCRIPTS} ${SCRIPTS}/slurm/splitFastq.sh ${SAMPLE} ${f}
  done
fi
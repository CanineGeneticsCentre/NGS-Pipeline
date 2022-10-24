#!/usr/bin/env bash

#! RUN : bash upload-enaFastq.sh <SAMPLE>
#! Eg. : bash upload-enaFastq.sh BoT_37630

SAMPLE=$1
SCRIPTS=`dirname $0`

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided to upload"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE
rsync --progress -auvh ${WGS}/${SAMPLE}/*.fq.gz ./

echo sbatch -J ${SAMPLE}.ascp --export=SAMPLE=${SAMPLE},CONDA_PREFIX,ASPERA_SCP_PASS ${SCRIPTS}/slurm/ena-ascp.sh
sbatch -J ${SAMPLE}.ascp --export=SAMPLE=${SAMPLE},CONDA_PREFIX,ASPERA_SCP_PASS ${SCRIPTS}/slurm/ena-ascp.sh

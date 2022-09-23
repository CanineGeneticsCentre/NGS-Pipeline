#!/usr/bin/env bash

#! RUN : bash upload-enaFastq.sh <SAMPLE>
#! Eg. : bash upload-enaFastq.sh BoT_37630

SAMPLE=$1

[[ -z "$SAMPLE" ]] && { echo "ERROR: No list of SAMPLES provided for this run"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE
rsync --progress -auvh ${WGS}/${SAMPLE}/*.fq.gz ./

sbatch -J ${SAMPLE}.ascp --export=SAMPLE=${SAMPLE},CONDA_PREFIX=${CONDA_PREFIX},ASPERA_SCP_PASS=${ASPERA_SCP_PASS} slurm/ena-ascp.sh

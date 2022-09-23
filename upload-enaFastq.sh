#!/usr/bin/env bash

SAMPLE=$1

mkdir -p $SAMPLE/logs; cd $SAMPLE
rsync --progress -auvh ${WGS}/${SAMPLE}/*.fq.gz ./

sbatch -J ${SAMPLE}.ascp --export=SAMPLE=${SAMPLE},CONDA_PREFIX=${CONDA_PREFIX},ASPERA_SCP_PASS=${ASPERA_SCP_PASS} slurm/ena-ascp.sh

#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE>
#! Eg. : bash ngs-pipeline.sh CS_35365

SAMPLE=$1
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $SAMPLE.config; source $SAMPLE.config

jid1=$(sbatch -J ${SAMPLE}.fastq2sam ${SCRIPTS}/slurm/fastq2sam.sh ${SAMPLE})
echo $jid1
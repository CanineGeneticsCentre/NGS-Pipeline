#!/usr/bin/env bash

#! RUN : bash genomics-db.sh <SAMPLE> <REF>
#! Eg. : bash genomics-db.sh sanples.list cf4

SAMPLE_LIST=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: No list of SAMPLES provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date
DIR=`echo $RANDOM | md5sum | head -c 10`

echo "Creating directory ${DIR}"
mkdir -p $DIR/logs; cd $DIR
mv ../${SAMPLE_LIST} .

cp $CFG ${REF}.config; source ${REF}.config


for s in `cat ${SAMPLE_LIST}`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* ./;
done

INTERVALS=`wc -l ${FASTA}/genomicsDB.intervals | awk '{print $1}'`

jid1=$(sbatch -A ${ACCOUNT} -J GenomicsDB --array=1-${INTERVALS} ${SCRIPTS}/slurm/gvcf2GenomicsDB.sh ${SAMPLES} ${REF})


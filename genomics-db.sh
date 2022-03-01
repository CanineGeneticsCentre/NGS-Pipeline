#!/usr/bin/env bash

#! RUN : bash genomics-db.sh <SAMPLE> <REF>
#! Eg. : bash genomics-db.sh sanples.list cf4

SAMPLE_LIST=$1
REF=$2
CHR=$3
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: No list of SAMPLES provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

#DIR='67c7afff2f'; cd $DIR; source ${REF}.config

DIR=`echo $RANDOM | md5sum | head -c 10`
echo "Creating directory ${DIR}"
mkdir -p $DIR/logs; cd $DIR
mv ../${SAMPLE_LIST} .
cp $CFG ${REF}.config; source ${REF}.config


for s in `cat ${SAMPLE_LIST}`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* ./;
done

if [ -z "$CHR" ]; then
  INTERVALS=`wc -l ${FASTA}/genomicsDB.intervals | awk '{print $1}'`
  ARRAY="1-${INTERVALS}"
else
  ARRAY="${CHR}-${CHR}"
fi

echo sbatch -A ${ACCOUNT} -J GenomicsDB --array=${ARRAY} ${SCRIPTS}/slurm/gvcf2GenomicsDB.sh ${SAMPLE_LIST} ${REF}
sbatch -A ${ACCOUNT} -J GenomicsDB --array=${ARRAY} ${SCRIPTS}/slurm/gvcf2GenomicsDB.sh ${SAMPLE_LIST} ${REF}


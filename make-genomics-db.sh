#!/usr/bin/env bash

#! RUN : bash make-genomics-db.sh <SAMPLE> <REF> [<CHR>]
#! Eg. : bash make-genomics-db.sh samples.list cf4 [38]

SAMPLE_LIST=$1
REF=$2
CHR=$3
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: No list of SAMPLES provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

#DIR='67c7afff2f'; cd $DIR; source ${REF}.config
#DIR=`echo $RANDOM | md5sum | head -c 10`
DIR=$(basename $SAMPLE_LIST .list)
echo "Creating directory ${DIR}"
mkdir -p $DIR/logs; cd $DIR
mv ../${SAMPLE_LIST} .
cp $CFG ${REF}.config; source ${REF}.config


for s in `cat ${SAMPLE_LIST}`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* ./;
done


# If chromosome is NOT set then use all lines in intervals list to scatter the update job
if [ -z "$CHR" ]; then
  INTERVALS=`wc -l ${FASTA}/${REF}-genomicsDB.intervals | awk '{print $1}'`
  ARRAY="1-${INTERVALS}"
else
  MIN=`grep -n "^$CHR:" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | head -1`
  MAX=`grep -n "^$CHR:" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | tail -1`
  ARRAY="${MIN}-${MAX}"
fi

echo sbatch -A ${ACCOUNT} -J GenomicsDB --array=${ARRAY} ${SCRIPTS}/slurm/createGenomicsDB.sh ${SAMPLE_LIST} ${REF}
sbatch -A ${ACCOUNT} -J GenomicsDB --array=${ARRAY} ${SCRIPTS}/slurm/createGenomicsDB.sh ${SAMPLE_LIST} ${REF}


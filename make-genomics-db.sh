#!/usr/bin/env bash

#! RUN : bash make-genomics-db.sh <SAMPLES MAP> <REF>
#! Eg. : bash make-genomics-db.sh samples.map cf4

SAMPLE_MAP=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE_MAP" ]] && { echo "ERROR: No SAMPLE map provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }


DIR=$(basename $SAMPLE_MAP .map)
echo "Creating directory ${DIR}"
mkdir $DIR; cd $DIR
mkdir logs gvcf

cp ../${SAMPLE_MAP} .
cp $CFG ${REF}.config; source ${REF}.config

#INTERVALS=`wc -l ${FASTA}/${REF}-genomicsDB.intervals | awk '{print $1}'`
#Get total number of folder we will be making - one per CHR (gdb-chromsomes.intervals) plus a number for the scaffolds (gdb-scaffolds.intervals)
INTERVALS=`cat ${FASTA}/intervals/gdb*.intervals | wc -l`

for s in `cat ${SAMPLE_MAP}| cut -f 1`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* gvcf/;
done

sbatch -A ${ACCOUNT} -J GenomicsDB --array=1-$INTERVALS ${SCRIPTS}/slurm/createGenomicsDB.sh ${SAMPLE_MAP} ${REF}
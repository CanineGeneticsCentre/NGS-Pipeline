#!/usr/bin/env bash

#! RUN : bash update-genomics-db.sh <SAMPLE> <REF> [<CHR>]
#! Eg. : bash update-genomics-db.sh samples.list cf4 [38]

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

#Get total number of folder we will be making - one per CHR (gdb-chromsomes.intervals) plus a number for the scaffolds (gdb-scaffolds.intervals)
####INTERVALS=`cat ${FASTA}/intervals/gdb*.intervals | wc -l`
INTERVALS=1

for s in `cat ${SAMPLE_MAP}| cut -f 1`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* gvcf/;
done

sbatch -A ${ACCOUNT} -J GenomicsDB.update --array=1-$INTERVALS ${SCRIPTS}/slurm/updateGenomicsDB.sh ${SAMPLE_MAP} ${REF}
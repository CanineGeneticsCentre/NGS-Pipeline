#!/usr/bin/env bash

#! RUN : bash export-genomics-db.sh <REF> [<CHR>]
#! Eg. : bash export-genomics-db.sh cf4 [38]

REF=$1
CHR=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

source ${CFG};
mkdir -p ${GENOME}/logs; cd $GENOME
mkdir -p snpEff
cp $CFG ${REF}.config; source ${REF}.config


# If chromosome is NOT set then use all lines in intervals list to scatter the update job
if [ -z "$CHR" ]; then
  INTERVALS=`wc -l ${FASTA}/${REF}-genomicsDB.intervals | awk '{print $1}'`
  ARRAY="1-${INTERVALS}"
else
  MIN=`grep -n "^$CHR:" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | head -1`
  MAX=`grep -n "^$CHR:" $FASTA/${REF}-genomicsDB.intervals | awk  -F':' ' { print $1 } ' | tail -1`
  ARRAY="${MIN}-${MAX}"
fi

#echo sbatch -A ${ACCOUNT} -J ${REF}.VCF --array=${ARRAY} --export=SCRIPTS=${SCRIPTS} ${SCRIPTS}/slurm/GenomicsDB2vcf.sh ${REF}
jid=$(sbatch -A ${ACCOUNT} -J ${REF}.VCF --array=${ARRAY} --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/GenomicsDB2vcf.sh)


if [ -z "$CHR" ]; then
  for chr in `cut -f 1 -d':' ${FASTA}/${REF}-genomicsDB.intervals | cut -f 1 -d' ' | cut -d'_' -f 1 | sort -n | uniq`; do
    sbatch -A ${ACCOUNT} -J ${REF}-$chr.snpEff --dependency=afterok:${jid1##* } --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/annotateVcf.sh $chr;
  }
else
  sbatch -A ${ACCOUNT} -J ${REF}-${CHR}.snpEff --dependency=afterok:${jid1##* } --export=SCRIPTS=${SCRIPTS},REF=${REF} ${SCRIPTS}/slurm/annotateVcf.sh $CHR;
fi


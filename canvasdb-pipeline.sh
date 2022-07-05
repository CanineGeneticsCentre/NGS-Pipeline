#!/usr/bin/env bash

#! RUN : bash canvasdb-pipeline.sh <SAMPLE> <REF>
#! Eg. : bash canvasdb-pipeline.sh CS_35365 cf4

SAMPLE=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }

date

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $REF.config; source $REF.config

# Copy files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/${SAMPLE}-${REF}.g.vcf.gz* ./

# WORKING ON CHR 29 only...
# bcftools view CS_10842-cf3.g.vcf.gz --regions 29 > CS_10842-cf3-29.g.vcf
# bgzip CS_10842-cf3-29.g.vcf
# tabix -p vcf CS_10842-cf3-29.g.vcf.gz



jid=$(sbatch -A ${ACCOUNT} -J ${REF}.VCF --export=SAMPLE=${SAMPLE},REF=${REF} ${SCRIPTS}/slurm/gvcf2vcf.sh)
jid2=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.filterVCF --dependency=afterok:${jid1##* } --export=SAMPLE=${SAMPLE},REF=${REF} ${SCRIPTS}/slurm/filterVCF.sh ${SAMPLE})
jid3=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.vcfStats --dependency=afterok:${jid2##* } --export=SAMPLE=${SAMPLE},REF=${REF} ${SCRIPTS}/slurm/vcfStats.sh ${SAMPLE})
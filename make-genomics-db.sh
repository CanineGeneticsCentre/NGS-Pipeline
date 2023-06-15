#!/usr/bin/env bash

#! RUN : bash make-genomics-db.sh <SAMPLE> <REF> [<CHR>]
#! Eg. : bash make-genomics-db.sh samples.list cf4 [38]

SAMPLE_MAP=$1
REF=$2
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline-${REF}.config"

[[ -z "$SAMPLE_MAP" ]] && { echo "ERROR: No SAMPLE map provided for this run"; exit 1; }
[[ -z "$REF" ]] && { echo "ERROR: No REFERENCE provided for this run"; exit 1; }


DIR=$(basename $SAMPLE_MAP .map)
echo "Creating directory ${DIR}"
mkdir -p $DIR/logs; cd $DIR
cp ../${SAMPLE_MAP} .
cp $CFG ${REF}.config; source ${REF}.config

module load ${GATK}
#gatk SplitIntervals -R ${FASTA}/${GENOME}.fasta -L ${INTERVAL_LIST} --scatter-count ${INTERVALS} -O intervals --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION


for s in `cat ${SAMPLE_MAP}| cut -f 1`; do
  # Copy g.vcf files from RCS... will only copy if files no present or rcs version is newer
  rsync --progress -av ${WGS}/${s}/${s}-${REF}.g.vcf* ./;
done

#sbatch -A ${ACCOUNT} -J GenomicsDB --array=0-$(($INTERVALS-1)) ${SCRIPTS}/slurm/createGenomicsDB.sh ${SAMPLE_MAP} ${REF}


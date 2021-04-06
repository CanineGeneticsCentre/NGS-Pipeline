#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE>
#! Eg. : bash ngs-pipeline.sh CS_35365

SAMPLE=$1
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $SAMPLE.config; source $SAMPLE.config

if [[ ! -e ${SAMPLE}_R1.fastq.gz && -e ${SAMPLE}_R2.fastq.gz ]]
then 
  cp ${WGS}/${SAMPLE}/${SAMPLE}_R1.fastq.gz ${WGS}/${SAMPLE}/${SAMPLE}_R2.fastq.gz ./
fi


# fastq2sam - convert FASTQ to aligned SAM files
jid1=$(sbatch -J ${SAMPLE}.fastq2sam ${SCRIPTS}/slurm/fastq2sam.sh ${SAMPLE})
#jid2=$(sbatch -J ${SAMPLE}.sam2bam --dependency=afterok:${jid1##* } ${SCRIPTS}/slurm/sam2bam.sh ${SAMPLE} ${RUN_NAME})
#jid3=$(sbatch -J ${SAMPLE}.validateSam --dependency=afterok:${jid2##* } ${SCRIPTS}/slurm/validateSam.sh ${SAMPLE})
#jid4=$(sbatch -J ${SAMPLE}.markDuplicates --dependency=afterok:${jid3##* } ${SCRIPTS}/slurm/markDuplicates.sh ${SAMPLE})
#jid5=$(sbatch -J ${SAMPLE}.realignReads --dependency=afterok:${jid4##* } ${SCRIPTS}/slurm/realignReads.sh ${SAMPLE})
#jid6=$(sbatch -J ${SAMPLE}.realignIndels --dependency=afterok:${jid5##* } ${SCRIPTS}/slurm/realignIndels.sh ${SAMPLE})
#jid7=$(sbatch -J ${SAMPLE}.baseRecalibrator --dependency=afterok:${jid6##* } ${SCRIPTS}/slurm/baseRecalibrator.sh ${SAMPLE})
#jid8=$(sbatch -J ${SAMPLE}.printReads --dependency=afterok:${jid7##* } ${SCRIPTS}/slurm/printReads.sh ${SAMPLE})
#jid9=$(sbatch -J ${SAMPLE}.CollectMetrics --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/collectMetrics.sh ${SAMPLE})
#jid10=$(sbatch -J ${SAMPLE}.validateSam2 --dependency=afterok:${jid8##* } ${SCRIPTS}/slurm/validateSam2.sh ${SAMPLE})
#jid11=$(sbatch -J ${SAMPLE}.FinalStep --dependency=afterok:${jid9##* }:${jid10##* } ${SCRIPTS}/slurm/finalStep.sh ${SAMPLE})

echo $jid2
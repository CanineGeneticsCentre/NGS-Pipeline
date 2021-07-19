#!/usr/bin/env bash

#! RUN : bash ngs-pipeline.sh <SAMPLE>
#! Eg. : bash ngs-pipeline.sh CS_35365

SAMPLE=$1
SCRIPTS=`dirname $0`
CFG="${SCRIPTS}/ngs-pipeline.config"

[[ -z "$SAMPLE" ]] && { echo "ERROR: No SAMPLE provided for this run"; exit 1; }

mkdir -p $SAMPLE/logs; cd $SAMPLE
cp $CFG $SAMPLE.config; source $SAMPLE.config
mkdir -p $GENOME; 

# Copy files from RCS... will only copy if files no present or rcs version is newer
rsync --progress -av ${WGS}/${SAMPLE}/*.fq.gz ./


#if [[ ! -e ${SAMPLE}\_R1.fastq.gz && ! -e ${SAMPLE}\_R2.fastq.gz ]]
#then 
#  echo "You need to copy the FASTQ files into ${SAMPLE} directory before running this pipeline."
#  echo "FASTQ files can be found here - ${WGS}/${SAMPLE}"
#  echo "\t cp ${WGS}/${SAMPLE}/${SAMPLE}\_R1.fastq.gz ${SAMPLE}"
#  echo "\t cp ${WGS}/${SAMPLE}/${SAMPLE}\_R2.fastq.gz ${SAMPLE}"
#  exit 1;
#fi



# Count how many fastq files we have
COUNT=`ls *.fq.gz | wc -l`
# Work out number of lanes used... assuming pair-end, therefore divide total number of files by 2 - forward and reverse files for each lane
LANES=$((COUNT / 2))

cd $GENOME
# fastq2sam - Submit job array to align samples to ref genome
jid1=$(sbatch -A ${ACCOUNT} -J ${SAMPLE}.fastq2sam --array=1-${LANES} ${SCRIPTS}/slurm/fastq2sam.sh ${SAMPLE})




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

echo $jid1
